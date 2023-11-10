using Test
using Pkg.Artifacts
using NCDatasets
import JET
import ClimaComms
import Infiltrator

using RRTMGP
using RRTMGP.Vmrs
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.AngularDiscretizations
using RRTMGP.RTE
using RRTMGP.RTESolver

import CLIMAParameters as CP

include(joinpath(pkgdir(RRTMGP), "parameters", "create_parameters.jl"))
# overriding some parameters to match with RRTMGP FORTRAN code

include("reference_files.jl")
include("read_rfmip_clear_sky.jl")
#---------------------------------------------------------------
function clear_sky(
    context,
    ::Type{OPC},
    ::Type{SRC},
    ::Type{VMR},
    ::Type{FT};
    exfiltrate = false,
) where {FT <: AbstractFloat, OPC, SRC, VMR}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = create_insolation_parameters(FT, overrides)
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    opc = Symbol(OPC)
    lw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :lw) # lw lookup tables for gas optics
    sw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :sw) # sw lookup tables for gas optics
    input_file = get_ref_filename(:atmos_state, :clearsky)         # clear-sky atmos state
    # reference data files for comparison
    flux_up_file_lw = get_ref_filename(:comparison, :clearsky, λ = :lw, flux_up_dn = :flux_up, opc = opc)
    flux_dn_file_lw = get_ref_filename(:comparison, :clearsky, λ = :lw, flux_up_dn = :flux_dn, opc = opc)
    flux_up_file_sw = get_ref_filename(:comparison, :clearsky, λ = :sw, flux_up_dn = :flux_up, opc = opc)
    flux_dn_file_sw = get_ref_filename(:comparison, :clearsky, λ = :sw, flux_up_dn = :flux_dn, opc = opc)

    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    max_threads = 256
    exp_no = 1
    n_gauss_angles = 1

    # reading longwave lookup data
    ds_lw = Dataset(lw_file, "r")
    lookup_lw, idx_gases = LookUpLW(ds_lw, FT, DA)
    close(ds_lw)

    # reading shortwave lookup data
    ds_sw = Dataset(sw_file, "r")
    lookup_sw, idx_gases = LookUpSW(ds_sw, FT, DA)
    close(ds_sw)

    # reading rfmip data to atmospheric state
    ds_lw_in = Dataset(input_file, "r")

    (as, sfc_emis, sfc_alb_direct, zenith, toa_flux, usecol) =
        setup_rfmip_as(context, ds_lw_in, idx_gases, exp_no, lookup_lw, FT, VMR, max_threads, param_set)
    close(ds_lw_in)

    ncol, nlay, ngpt_lw = as.ncol, as.nlay, lookup_lw.n_gpt
    nlev = nlay + 1
    op = OPC(FT, ncol, nlay, DA) # allocating optical properties object

    # setting up longwave problem
    ngpt_lw = lookup_lw.n_gpt
    src_lw = source_func_longwave(param_set, FT, ncol, nlay, opc, DA)         # allocating longwave source function object
    bcs_lw = LwBCs{FT, typeof(sfc_emis), Nothing}(sfc_emis, nothing) # setting up boundary conditions
    ang_disc = nothing#AngularDiscretization(opc, FT, n_gauss_angles, DA)  # initializing Angular discretization
    fluxb_lw = FluxLW(ncol, nlay, FT, DA)                          # flux storage for bandwise calculations
    flux_lw = FluxLW(ncol, nlay, FT, DA)                           # longwave fluxes
    # setting up shortwave problem
    ngpt_sw = lookup_sw.n_gpt
    src_sw = source_func_shortwave(FT, ncol, nlay, opc, DA)        # allocating shortwave source function object
    inc_flux_diffuse = nothing
    sfc_alb_diffuse = FTA2D(deepcopy(sfc_alb_direct))
    bcs_sw = SwBCs(zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

    fluxb_sw = FluxSW(ncol, nlay, FT, DA) # flux storage for bandwise calculations
    flux_sw = FluxSW(ncol, nlay, FT, DA)  # shortwave fluxes for band calculations
    #--------------------------------------------------
    # initializing RTE solver
    slv = Solver(context, as, op, src_lw, src_sw, bcs_lw, bcs_sw, fluxb_lw, fluxb_sw, flux_lw, flux_sw)
    exfiltrate && Infiltrator.@exfiltrate
    println("calling longwave solver; ncol = $ncol")
    @time solve_lw!(slv, max_threads, lookup_lw)
    if device isa ClimaComms.CPUSingleThreaded
        JET.@test_opt solve_lw!(slv, max_threads, lookup_lw)
        @test_broken (@allocated solve_lw!(slv, max_threads, lookup_lw)) == 0
        @test (@allocated solve_lw!(slv, max_threads, lookup_lw)) ≤ 448
    end

    println("calling shortwave solver; ncol = $ncol")
    @time solve_sw!(slv, max_threads, lookup_sw)
    if device isa ClimaComms.CPUSingleThreaded
        JET.@test_opt solve_sw!(slv, max_threads, lookup_sw)
        @test_broken (@allocated solve_sw!(slv, max_threads, lookup_sw)) == 0
        @test (@allocated solve_sw!(slv, max_threads, lookup_sw)) ≤ 448
    end

    # comparing longwave fluxes with data from RRTMGP FORTRAN code
    flip_ind = nlev:-1:1

    ds_flux_up_lw = Dataset(flux_up_file_lw, "r")
    comp_flux_up_lw = Array(ds_flux_up_lw["rlu"])[flip_ind, :, exp_no]
    close(ds_flux_up_lw)

    ds_flux_dn_lw = Dataset(flux_dn_file_lw, "r")
    comp_flux_dn_lw = Array(ds_flux_dn_lw["rld"])[flip_ind, :, exp_no]
    close(ds_flux_dn_lw)

    flux_up_lw = Array(slv.flux_lw.flux_up)
    flux_dn_lw = Array(slv.flux_lw.flux_dn)

    max_err_flux_up_lw = maximum(abs.(flux_up_lw .- comp_flux_up_lw))
    max_err_flux_dn_lw = maximum(abs.(flux_dn_lw .- comp_flux_dn_lw))

    println("=======================================")
    println("Clear-sky longwave test - $opc")
    println("max_err_flux_up_lw = $max_err_flux_up_lw")
    println("max_err_flux_dn_lw = $max_err_flux_dn_lw")

    toler_lw = 1e-4
    #--------------------------------------------------------------
    # comparing shortwave fluxes with data from RRTMGP FORTRAN code
    ds_flux_up_sw = Dataset(flux_up_file_sw, "r")
    comp_flux_up_sw = Array(ds_flux_up_sw["rsu"])[flip_ind, :, exp_no]
    close(ds_flux_up_sw)

    ds_flux_dn_sw = Dataset(flux_dn_file_sw, "r")
    comp_flux_dn_sw = Array(ds_flux_dn_sw["rsd"])[flip_ind, :, exp_no]
    close(ds_flux_dn_sw)

    flux_up_sw = Array(slv.flux_sw.flux_up)
    flux_dn_sw = Array(slv.flux_sw.flux_dn)

    for i in 1:ncol
        if usecol[i] == 0
            flux_up_sw[:, i] .= FT(0)
            flux_dn_sw[:, i] .= FT(0)
        end
    end

    max_err_flux_up_sw = maximum(abs.(flux_up_sw .- comp_flux_up_sw))
    max_err_flux_dn_sw = maximum(abs.(flux_dn_sw .- comp_flux_dn_sw))

    println("Clear-sky shortwave test, opc  = $opc")
    println("max_err_flux_up_sw = $max_err_flux_up_sw")
    println("max_err_flux_dn_sw = $max_err_flux_dn_sw")
    println("=======================================")


    toler_sw = FT(0.001)

    @test max_err_flux_up_lw ≤ toler_lw broken = (FT == Float32)
    @test max_err_flux_dn_lw ≤ toler_lw broken = (FT == Float32)
    @test max_err_flux_up_sw ≤ toler_sw broken = (FT == Float32)
    @test max_err_flux_dn_sw ≤ toler_sw broken = (FT == Float32)
end
