using Test
using Pkg.Artifacts
using NCDatasets
using ClimaComms

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
function lw_rfmip(context, ::Type{OPC}, ::Type{SRC}, ::Type{VMR}, ::Type{FT}) where {FT <: AbstractFloat, OPC, SRC, VMR}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = create_insolation_parameters(FT, overrides)
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    opc = Symbol(OPC)
    lw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :lw) # lw lookup tables
    lw_input_file = get_ref_filename(:atmos_state, :clearsky)      # clear-sky atmos state
    # reference data files for comparison
    flux_up_file = get_ref_filename(:comparison, :clearsky, λ = :lw, flux_up_dn = :flux_up, opc = opc)
    flux_dn_file = get_ref_filename(:comparison, :clearsky, λ = :lw, flux_up_dn = :flux_dn, opc = opc)

    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    max_threads = Int(256)
    exp_no = 1
    n_gauss_angles = 1

    # reading longwave lookup data
    ds_lw = Dataset(lw_file, "r")
    lookup_lw, idx_gases = LookUpLW(ds_lw, FT, DA)
    close(ds_lw)
    # reading rfmip data to atmospheric state
    ds_lw_in = Dataset(lw_input_file, "r")

    (as, sfc_emis, _, _, _) =
        setup_rfmip_as(context, ds_lw_in, idx_gases, exp_no, lookup_lw, FT, VMR, max_threads, param_set)
    close(ds_lw_in)

    ncol, nlay, ngpt = as.ncol, as.nlay, lookup_lw.n_gpt
    nlev = nlay + 1
    op = OPC(FT, ncol, nlay, DA)                                   # allocating optical properties object
    src_lw = source_func_longwave(param_set, FT, ncol, nlay, opc, DA) # allocating longwave source function object
    bcs_lw = LwBCs{FT, typeof(sfc_emis), Nothing}(sfc_emis, nothing)            # setting up boundary conditions
    fluxb_lw = FluxLW(ncol, nlay, FT, DA)                          # flux storage for bandwise calculations
    flux_lw = FluxLW(ncol, nlay, FT, DA)                           # longwave fluxes

    # initializing RTE solver
    slv = Solver(context, as, op, src_lw, nothing, bcs_lw, nothing, fluxb_lw, nothing, flux_lw, nothing)

    solve_lw!(slv, max_threads, lookup_lw)
    # comparing with data from RRTMGP FORTRAN code
    flip_ind = nlev:-1:1

    ds_flux_up = Dataset(flux_up_file, "r")
    comp_flux_up = Array{FT}(Array(ds_flux_up["rlu"])[flip_ind, :, exp_no])
    close(ds_flux_up)

    ds_flux_dn = Dataset(flux_dn_file, "r")
    comp_flux_dn = Array{FT}(Array(ds_flux_dn["rld"])[flip_ind, :, exp_no])
    close(ds_flux_dn)

    max_err_flux_up = FT(maximum(abs.(Array(slv.flux_lw.flux_up) .- comp_flux_up)))
    max_err_flux_dn = FT(maximum(abs.(Array(slv.flux_lw.flux_dn) .- comp_flux_dn)))

    println("=======================================")
    println("Clear-sky longwave test - $opc")
    println("max_err_flux_up = $max_err_flux_up")
    println("max_err_flux_dn = $max_err_flux_dn")
    println("=======================================")
    toler = FT(1e-4)

    @test max_err_flux_up ≤ toler
    @test max_err_flux_dn ≤ toler
    return nothing
end

lw_rfmip(ClimaComms.context(), OneScalar, SourceLWNoScat, VmrGM, Float64)
lw_rfmip(ClimaComms.context(), TwoStream, SourceLW2Str, VmrGM, Float64)
