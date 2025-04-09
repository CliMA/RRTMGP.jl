using Test
using Pkg.Artifacts
using NCDatasets

import JET
import Infiltrator
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

using RRTMGP
using RRTMGP: RRTMGPGridParams, RRTMGPSolver
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
using RRTMGP.ArtifactPaths
# overriding some parameters to match with RRTMGP FORTRAN code

@isdefined(api_methods) || include("api_method_utils.jl")
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
    lookup_lw, idx_gases = Dataset(lw_file, "r") do ds
        LookUpLW(ds, FT, DA)
    end
    # reading longwave cloud lookup data
    lookup_lw_cld = Dataset(lw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end
    # reading longwave aerosol lookup data
    lookup_lw_aero, idx_aerosol, idx_aerosize = Dataset(lw_aero_file, "r") do ds
        LookUpAerosolMerra(ds, FT, DA)
    end

    #reading shortwave gas optics lookup data
    lookup_sw, idx_gases = Dataset(sw_file, "r") do ds
        LookUpSW(ds, FT, DA)
    end
    # reading longwave cloud lookup data
    lookup_sw_cld = Dataset(sw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end

    # reading shortwave aerosol lookup data
    lookup_sw_aero, _, _ = Dataset(sw_aero_file, "r") do ds
        LookUpAerosolMerra(ds, FT, DA)
    end

    # reading input file
    ds_in = Dataset(input_file, "r")
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, cos_zenith, toa_flux, bot_at_1 = setup_allsky_with_aerosols_as(
        context,
        ds_in,
        idx_gases,
        idx_aerosol,
        idx_aerosize,
        lookup_lw,
        lookup_sw,
        lookup_lw_cld,
        lookup_sw_cld,
        cldfrac,
        ncol,
        FT,
        param_set,
    )
    close(ds_in)

    nlay, _ = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)
    # Setting up longwave problem
    inc_flux = nothing
    slv_lw = SLVLW(grid_params; params = param_set, sfc_emis, inc_flux)
    # Setting up shortwave problem
    inc_flux_diffuse = nothing
    swbcs = (; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(grid_params; swbcs...)

    #---------------- Exercise new api (start)
    bcs_lw = BCs.LwBCs(sfc_emis, inc_flux)
    bcs_sw = BCs.SwBCs(cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    radiation_method = RRTMGP.AllSkyRadiationWithClearSkyDiagnostics(
        true, # aerosol_radiation
        true, # reset_rng_seed
    )
    solver = RRTMGPSolver(grid_params, radiation_method, param_set, bcs_lw, bcs_sw, as)
    RRTMGP.update_sw_fluxes!(solver)
    RRTMGP.update_lw_fluxes!(solver)
    for m in api_methods
        getproperty(RRTMGP, m)(solver)
    end
    for name in RRTMGP.aerosol_names()
        RRTMGP.aero_column_mass_density(solver, name)
        RRTMGP.aero_radius(solver, name)
    end
    #---------------- Exercise new api (end)

    # calling solvers

    # unity metric scaling - i.e. shallow atmosphere approximation (no column expansion with height)
    # test scaling factors (e.g. when applying corrections to metric terms for deep atmospheres)
    # first, test do-nothing op eg. shallow atmospheres
    metric_scaling = nothing
    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, lookup_lw_aero, metric_scaling)
    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld, lookup_sw_aero, metric_scaling)

    # comparison
    method = "Lookup Table Interpolation method"
    comp_flux_up_lw, comp_flux_dn_lw, comp_flux_up_sw, comp_flux_dn_sw = load_comparison_data(bot_at_1, ncol)

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
    printstyled(
        "Cloudy-sky with aerosols longwave test with ncol = $ncol, nlev = $nlev, Solver = $SLVLW, FT = $FT\n",
        color = color2,
    )
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
        "Cloudy-sky with aerosols shortwave test with ncol = $ncol, nlev = $nlev, Solver = $SLVSW, FT = $FT\n",
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

    @test minimum(as.aerosol_state.aod_sw_ext) >= 0
    @test minimum(as.aerosol_state.aod_sw_sca) >= 0
    @test minimum(as.aerosol_state.aod_sw_ext .- as.aerosol_state.aod_sw_sca) >= 0

    # New problem instance for metric scaling test
    # Setting up longwave problem

    # Set up test variables
    test_flux_up_sw = deepcopy(DA(slv_sw.flux.flux_up))
    test_flux_dn_sw = deepcopy(DA(slv_sw.flux.flux_dn))
    test_flux_dn_dir_sw = deepcopy(DA(slv_sw.flux.flux_dn_dir))
    test_flux_net_sw = deepcopy(DA(slv_sw.flux.flux_net))

    test_flux_up_lw = deepcopy(DA(slv_lw.flux.flux_up))
    test_flux_dn_lw = deepcopy(DA(slv_lw.flux.flux_dn))
    test_flux_net_lw = deepcopy(DA(slv_lw.flux.flux_net))
    # Set up problem 
    inc_flux = nothing
    slv_lw = SLVLW(FT, DA, context, param_set, nlay, ncol, sfc_emis, inc_flux)
    # Setting up shortwave problem
    inc_flux_diffuse = nothing
    swbcs = (cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(FT, DA, context, nlay, ncol, swbcs...)

    metric_scaling = DA(one.(slv_sw.flux.flux_up) * FT(2))
    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, lookup_lw_aero, metric_scaling)
    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld, lookup_sw_aero, metric_scaling)

    flux_up_sw = DA(slv_sw.flux.flux_up)
    flux_dn_sw = DA(slv_sw.flux.flux_dn)
    flux_net_sw = DA(slv_sw.flux.flux_net)
    flux_up_lw = DA(slv_lw.flux.flux_up)
    flux_dn_lw = DA(slv_lw.flux.flux_dn)
    flux_net_lw = DA(slv_lw.flux.flux_net)
    flux_dn_dir_sw = DA(slv_sw.flux.flux_dn_dir)

    @test all(test_flux_up_sw == flux_up_sw .* metric_scaling)
    @test all(test_flux_dn_sw == flux_dn_sw .* metric_scaling)
    @test all(test_flux_net_sw == flux_net_sw .* metric_scaling)
    @test all(test_flux_up_lw == flux_up_lw .* metric_scaling)
    @test all(test_flux_dn_lw == flux_dn_lw .* metric_scaling)
    @test all(test_flux_net_lw == flux_net_lw .* metric_scaling)
    @test all(test_flux_dn_dir_sw[1, :] == flux_dn_dir_sw[1, :] .* metric_scaling[1, :])

    return nothing
end
