using Test
using Pkg.Artifacts
using NCDatasets

import JET
import Infiltrator
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

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
include("read_all_sky_with_aerosols.jl")

function all_sky_with_aerosols(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{FT},
    toler_lw,
    toler_sw;
    ncol = 128,# repeats col#1 ncol times per RRTMGP example 
    use_lut::Bool = true,
    cldfrac = FT(1),
    exfiltrate = false,
) where {FT <: AbstractFloat, SLVLW, SLVSW}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)

    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    n_gauss_angles = 1

    lw_file = get_lookup_filename(:gas, :lw)          # lw lookup tables for gas optics
    lw_cld_file = get_lookup_filename(:cloud, :lw)    # lw cloud lookup tables
    lw_aero_file = get_lookup_filename(:aerosol, :lw) # lw aerosol lookup tables

    sw_file = get_lookup_filename(:gas, :sw)          # sw lookup tables for gas optics
    sw_cld_file = get_lookup_filename(:cloud, :sw)    # sw cloud lookup tables
    sw_aero_file = get_lookup_filename(:aerosol, :sw) # sw aerosol lookup tables

    input_file = get_input_filename(:gas_clouds_aerosols, :lw) # all-sky atmos state

    #reading longwave gas optics lookup data
    ds_lw = Dataset(lw_file, "r")
    lookup_lw, idx_gases = LookUpLW(ds_lw, FT, DA)
    close(ds_lw)
    # reading longwave cloud lookup data
    ds_lw_cld = Dataset(lw_cld_file, "r")
    lookup_lw_cld = use_lut ? LookUpCld(ds_lw_cld, FT, DA) : PadeCld(ds_lw_cld, FT, DA)
    close(ds_lw_cld)
    # reading longwave aerosol lookup data
    ds_lw_aero = Dataset(lw_aero_file, "r")
    lookup_lw_aero = LookUpAerosolMerra(ds_lw_aero, FT, DA)
    close(ds_lw_aero)

    #reading shortwave gas optics lookup data
    ds_sw = Dataset(sw_file, "r")
    lookup_sw, idx_gases = LookUpSW(ds_sw, FT, DA)
    close(ds_sw)
    # reading longwave cloud lookup data
    ds_sw_cld = Dataset(sw_cld_file, "r")
    lookup_sw_cld = use_lut ? LookUpCld(ds_sw_cld, FT, DA) : PadeCld(ds_sw_cld, FT, DA)
    close(ds_sw_cld)

    # reading shortwave aerosol lookup data
    ds_sw_aero = Dataset(sw_aero_file, "r")
    lookup_sw_aero = LookUpAerosolMerra(ds_sw_aero, FT, DA)
    close(ds_sw_aero)

    # reading input file 
    ds_in = Dataset(input_file, "r")
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, cos_zenith, toa_flux, bot_at_1 = setup_allsky_with_aerosols_as(
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
    inc_flux = nothing
    slv_lw = SLVLW(FT, DA, context, param_set, nlay, ncol, sfc_emis, inc_flux)
    # Setting up shortwave problem---------------------------------------
    inc_flux_diffuse = nothing
    swbcs = (cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(FT, DA, context, nlay, ncol, swbcs...)
    #------calling solvers
    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, lookup_lw_aero)
    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld, lookup_sw_aero)
    #------------------------------------------------------------------
    # comparison
    method = use_lut ? "Lookup Table Interpolation method" : "PADE method"
    comp_flux_up_lw, comp_flux_dn_lw, comp_flux_up_sw, comp_flux_dn_sw = load_comparison_data(use_lut, bot_at_1, ncol)

    comp_flux_net_lw = comp_flux_up_lw .- comp_flux_dn_lw
    comp_flux_net_sw = comp_flux_up_sw .- comp_flux_dn_sw

    flux_up_lw = Array(slv_lw.flux.flux_up)
    flux_dn_lw = Array(slv_lw.flux.flux_dn)
    flux_net_lw = Array(slv_lw.flux.flux_net)

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
    printstyled("Cloudy-sky longwave test with ncol = $ncol, nlev = $nlev, Solver = $SLVLW, FT = $FT\n", color = color2)
    printstyled("device = $device\n", color = color2)
    printstyled("$method\n\n", color = color2)
    println("L∞ error in flux_up           = $max_err_flux_up_lw")
    println("L∞ error in flux_dn           = $max_err_flux_dn_lw")
    println("L∞ error in flux_net          = $max_err_flux_net_lw")
    println("L∞ relative error in flux_net = $(max_rel_err_flux_net_lw * 100) %\n")

    flux_up_sw = Array(slv_sw.flux.flux_up)
    flux_dn_sw = Array(slv_sw.flux.flux_dn)
    flux_dn_dir_sw = Array(slv_sw.flux.flux_dn_dir)
    flux_net_sw = Array(slv_sw.flux.flux_net)

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

    printstyled(
        "Cloudy-sky shortwave test with ncol = $ncol, nlev = $nlev, Solver = $SLVSW, FT = $FT\n",
        color = color2,
    )
    printstyled("device = $device\n", color = color2)
    printstyled("$method\n\n", color = color2)
    println("L∞ error in flux_up           = $max_err_flux_up_sw")
    println("L∞ error in flux_dn           = $max_err_flux_dn_sw")
    println("L∞ error in flux_net          = $max_err_flux_net_sw")
    println("L∞ relative error in flux_net = $(max_rel_err_flux_net_sw * 100) %\n")

    # The reference results for the longwave solver are generated using a non-scattering solver,
    # which differ from the results generated by the TwoStream currently used.
    @test max_err_flux_up_lw ≤ toler_lw[FT]
    @test max_err_flux_dn_lw ≤ toler_lw[FT]
    @test max_err_flux_net_lw ≤ toler_lw[FT]

    @test max_err_flux_up_sw ≤ toler_sw[FT]
    @test max_err_flux_dn_sw ≤ toler_sw[FT]
    @test max_err_flux_net_sw ≤ toler_sw[FT]

    return nothing
end
