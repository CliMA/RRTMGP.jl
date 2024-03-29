using Test
using Pkg.Artifacts
using NCDatasets

import JET
import Infiltrator
import ClimaComms
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
import RRTMGP.Parameters.RRTMGPParameters
import ClimaParams as CP
# overriding some parameters to match with RRTMGP FORTRAN code

include("reference_files.jl")
include("read_all_sky.jl")

function all_sky(
    context,
    ::Type{OPC},
    ::Type{FT};
    ncol = 128,# repeats col#1 ncol times per RRTMGP example 
    use_lut::Bool = true,
    cldfrac = FT(1),
    exfiltrate = false,
) where {FT <: AbstractFloat, OPC}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)

    opc = Symbol(OPC)
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    max_threads = 256
    n_gauss_angles = 1

    lw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :lw)             # lw lookup tables for gas optics
    lw_cld_file = get_ref_filename(:lookup_tables, :cloudysky, λ = :lw)        # lw cloud lookup tables
    sw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :sw)             # sw lookup tables for gas optics
    sw_cld_file = get_ref_filename(:lookup_tables, :cloudysky, λ = :sw)        # lw cloud lookup tables

    input_file = get_ref_filename(:atmos_state, :cloudysky)                    # all-sky atmos state

    #reading longwave gas optics lookup data
    ds_lw = Dataset(lw_file, "r")
    lookup_lw, idx_gases = LookUpLW(ds_lw, FT, DA)
    close(ds_lw)
    # reading longwave cloud lookup data
    ds_lw_cld = Dataset(lw_cld_file, "r")
    lookup_lw_cld = use_lut ? LookUpCld(ds_lw_cld, FT, DA) : PadeCld(ds_lw_cld, FT, DA)
    close(ds_lw_cld)
    #reading shortwave gas optics lookup data
    ds_sw = Dataset(sw_file, "r")
    lookup_sw, idx_gases = LookUpSW(ds_sw, FT, DA)
    close(ds_sw)
    # reading longwave cloud lookup data
    ds_sw_cld = Dataset(sw_cld_file, "r")
    lookup_sw_cld = use_lut ? LookUpCld(ds_sw_cld, FT, DA) : PadeCld(ds_sw_cld, FT, DA)
    close(ds_sw_cld)
    # reading input file 
    ds_in = Dataset(input_file, "r")
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, cos_zenith, toa_flux, bot_at_1 = setup_allsky_as(
        context,
        ds_in,
        idx_gases,
        lookup_lw,
        lookup_sw,
        lookup_lw_cld,
        lookup_sw_cld,
        cldfrac,
        use_lut,
        ncol,
        FT,
        max_threads,
        param_set,
    )
    close(ds_in)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    #---reading comparison files -----------------------------------
    if bot_at_1
        flip_ind = 1:nlev
    else
        flip_ind = nlev:-1:1
    end

    # Setting up longwave problem---------------------------------------
    op = OPC(FT, ncol, nlay, DA) # allocating optical properties object
    src_lw = source_func_longwave(param_set, FT, ncol, nlay, opc, DA)   # allocating longwave source function object
    bcs_lw = LwBCs{FT, typeof(sfc_emis), Nothing}(sfc_emis, nothing)    # setting up boundary conditions
    fluxb_lw = FluxLW(ncol, nlay, FT, DA)                             # flux storage for bandwise calculations
    flux_lw = FluxLW(ncol, nlay, FT, DA)                              # longwave fluxes

    # Setting up shortwave problem---------------------------------------
    src_sw = source_func_shortwave(FT, ncol, nlay, opc, DA)        # allocating longwave source function object

    # setting up boundary conditions
    inc_flux_diffuse = nothing

    bcs_sw = SwBCs(cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

    fluxb_sw = FluxSW(ncol, nlay, FT, DA) # flux storage for bandwise calculations
    flux_sw = FluxSW(ncol, nlay, FT, DA)  # shortwave fluxes for band calculations    
    #-------------------------------------------------------------------
    # initializing RTE solver
    slv = Solver(context, as, op, src_lw, src_sw, bcs_lw, bcs_sw, fluxb_lw, fluxb_sw, flux_lw, flux_sw)
    #------calling solvers
    solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld)
    if device isa ClimaComms.CPUSingleThreaded
        JET.@test_opt solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld)
        @test (@allocated solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld)) == 0
        @test (@allocated solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld)) ≤ 736
    end

    exfiltrate && Infiltrator.@exfiltrate
    solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld)
    if device isa ClimaComms.CPUSingleThreaded
        JET.@test_opt solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld)
        @test (@allocated solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld)) == 0
        @test (@allocated solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld)) ≤ 736
    end
    #-------------
    # comparison
    method = use_lut ? "Lookup Table Interpolation method" : "PADE method"
    comp_flux_up_lw, comp_flux_dn_lw, comp_flux_up_sw, comp_flux_dn_sw = load_comparison_data(use_lut, bot_at_1, ncol)

    comp_flux_net_lw = comp_flux_up_lw .- comp_flux_dn_lw
    comp_flux_net_sw = comp_flux_up_sw .- comp_flux_dn_sw

    flux_up_lw = Array(slv.flux_lw.flux_up)
    flux_dn_lw = Array(slv.flux_lw.flux_dn)
    flux_net_lw = Array(slv.flux_lw.flux_net)

    max_err_flux_up_lw = maximum(abs.(flux_up_lw .- comp_flux_up_lw))
    max_err_flux_dn_lw = maximum(abs.(flux_dn_lw .- comp_flux_dn_lw))
    max_err_flux_net_lw = maximum(abs.(flux_net_lw .- comp_flux_net_lw))

    rel_err_flux_net_lw = abs.(flux_net_lw .- comp_flux_net_lw)

    for gcol in 1:ncol, glev in 1:nlev
        den = abs(comp_flux_net_lw[glev, gcol])
        if den > 10 * eps(FT)
            rel_err_flux_net_lw[glev, gcol] /= den
        end
    end
    max_rel_err_flux_net_lw = maximum(rel_err_flux_net_lw)
    color2 = :cyan
    printstyled("Cloudy-sky longwave test with ncol = $ncol, nlev = $nlev, OPC = $opc, FT = $FT\n", color = color2)
    printstyled("device = $device\n", color = color2)
    printstyled("$method\n\n", color = color2)
    println("L∞ error in flux_up           = $max_err_flux_up_lw")
    println("L∞ error in flux_dn           = $max_err_flux_dn_lw")
    println("L∞ error in flux_net          = $max_err_flux_net_lw")
    println("L∞ relative error in flux_net = $(max_rel_err_flux_net_lw * 100) %\n")

    flux_up_sw = Array(slv.flux_sw.flux_up)
    flux_dn_sw = Array(slv.flux_sw.flux_dn)
    flux_dn_dir_sw = Array(slv.flux_sw.flux_dn_dir)
    flux_net_sw = Array(slv.flux_sw.flux_net)

    max_err_flux_up_sw = maximum(abs.(flux_up_sw .- comp_flux_up_sw))
    max_err_flux_dn_sw = maximum(abs.(flux_dn_sw .- comp_flux_dn_sw))
    max_err_flux_net_sw = maximum(abs.(flux_net_sw .- comp_flux_net_sw))

    rel_err_flux_net_sw = abs.(flux_net_sw .- comp_flux_net_sw)

    for gcol in 1:ncol, glev in 1:nlev
        den = abs(comp_flux_net_sw[glev, gcol])
        if den > 10 * eps(FT)
            rel_err_flux_net_sw[glev, gcol] /= den
        end
    end
    max_rel_err_flux_net_sw = maximum(rel_err_flux_net_sw)

    printstyled("Cloudy-sky shortwave test with ncol = $ncol, nlev = $nlev, OPC = $opc, FT = $FT\n", color = color2)
    printstyled("device = $device\n", color = color2)
    printstyled("$method\n\n", color = color2)
    println("L∞ error in flux_up           = $max_err_flux_up_sw")
    println("L∞ error in flux_dn           = $max_err_flux_dn_sw")
    println("L∞ error in flux_net          = $max_err_flux_net_sw")
    println("L∞ relative error in flux_net = $(max_rel_err_flux_net_sw * 100) %\n")

    toler = Dict(Float64 => Float64(1e-5), Float32 => Float32(0.05))

    @test max_err_flux_up_lw ≤ toler[FT]
    @test max_err_flux_dn_lw ≤ toler[FT]
    @test max_err_flux_net_lw ≤ toler[FT]

    @test max_err_flux_up_sw ≤ toler[FT]
    @test max_err_flux_dn_sw ≤ toler[FT]
    @test max_err_flux_net_sw ≤ toler[FT]
end
