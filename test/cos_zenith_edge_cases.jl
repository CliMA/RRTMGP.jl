#=
Test for numerical stability when cos_zenith is 0 or very small.

These edge cases could potentially cause issues with division by cos_zenith
in the shortwave solver equations (e.g., exp(-τ / cos_zenith)).
=#

using Test
using Pkg.Artifacts
using NCDatasets

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

include("reference_files.jl")
include("read_all_sky_with_aerosols.jl")

#=
    test_cos_zenith_edge_cases(context, ::Type{SLVLW}, ::Type{SLVSW}, ::Type{FT}; ncol=4)

Test that the solver handles edge cases for cos_zenith gracefully:
- cos_zenith = 0 (sun at horizon)
- cos_zenith = very small positive value (sun just above horizon)
- cos_zenith = negative value (nighttime)

The solver should not produce NaN or Inf values in any of these cases.

Note: flux_dn_dir only stores the direct beam flux at the surface (level 1).
Levels 2-nlev are not computed and may contain uninitialized values - this is 
expected behavior, not a bug (testing only level 1 for flux_dn_dir).
=#
function test_cos_zenith_edge_cases(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{FT};
    ncol = 4,
    cldfrac = FT(1),
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

    # Reading longwave gas optics lookup data
    lookup_lw, idx_gases = Dataset(lw_file, "r") do ds
        LookUpLW(ds, FT, DA)
    end
    # Reading longwave cloud lookup data
    lookup_lw_cld = Dataset(lw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end
    # Reading longwave aerosol lookup data
    lookup_lw_aero, idx_aerosol, idx_aerosize = Dataset(lw_aero_file, "r") do ds
        LookUpAerosolMerra(ds, FT, DA)
    end

    # Reading shortwave gas optics lookup data
    lookup_sw, idx_gases = Dataset(sw_file, "r") do ds
        LookUpSW(ds, FT, DA)
    end
    # Reading longwave cloud lookup data
    lookup_sw_cld = Dataset(sw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end

    # Reading shortwave aerosol lookup data
    lookup_sw_aero, _, _ = Dataset(sw_aero_file, "r") do ds
        LookUpAerosolMerra(ds, FT, DA)
    end

    # Reading input file
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

    # Copy aerosol data from column 1 to all other columns so all columns have identical data
    # except for the cos_zenith values being tested
    for icol in 2:ncol
        as.aerosol_state.aero_size[:, :, icol] .= as.aerosol_state.aero_size[:, :, 1]
        as.aerosol_state.aero_mass[:, :, icol] .= as.aerosol_state.aero_mass[:, :, 1]
    end

    # Test edge case cos_zenith values
    # Column 1: cos_zenith = 0.5 (normal positive, reference case)
    # Column 2: cos_zenith = 0 (sun exactly at horizon)
    # Column 3: cos_zenith = very small positive (sun just above horizon)  
    # Column 4: cos_zenith = negative (nighttime)
    cos_zenith_test = DA([FT(0.5), FT(0), FT(1e-10), FT(-0.5)])

    # Setting up longwave problem
    inc_flux = nothing
    slv_lw = SLVLW(grid_params; params = param_set, sfc_emis, inc_flux)

    # Setting up shortwave problem with edge case cos_zenith values
    inc_flux_diffuse = nothing
    swbcs = (; cos_zenith = cos_zenith_test, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(grid_params; swbcs...)

    # Calling solvers - this should not throw errors or produce NaN/Inf
    metric_scaling = nothing
    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, lookup_lw_aero, metric_scaling)
    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld, lookup_sw_aero, metric_scaling)

    # Retrieve fluxes
    flux_up_sw = Array(slv_sw.flux.flux_up)
    flux_dn_sw = Array(slv_sw.flux.flux_dn)
    flux_dn_dir_sw = Array(slv_sw.flux.flux_dn_dir)
    flux_net_sw = Array(slv_sw.flux.flux_net)

    flux_up_lw = Array(slv_lw.flux.flux_up)
    flux_dn_lw = Array(slv_lw.flux.flux_dn)
    flux_net_lw = Array(slv_lw.flux.flux_net)

    color2 = :cyan
    printstyled(
        "Edge case cos_zenith test with ncol = $ncol, nlev = $nlev, Solver = $SLVSW, FT = $FT\n",
        color = color2,
    )
    printstyled("device = $device\n", color = color2)
    printstyled("Testing cos_zenith = [0.5, 0, 1e-10, -0.5]\n\n", color = color2)

    # Test 1: No NaN/Inf in shortwave fluxes
    @test all(isfinite, flux_up_sw)
    @test all(isfinite, flux_dn_sw)
    @test all(isfinite, flux_net_sw)
    # Note: flux_dn_dir only has valid data at level 1 (surface) 
    # Levels 2-nlev may contain uninitialized values
    @test all(isfinite, flux_dn_dir_sw[1, :])

    # Test 2: No NaN/Inf in longwave fluxes
    @test all(isfinite, flux_up_lw)
    @test all(isfinite, flux_dn_lw)
    @test all(isfinite, flux_net_lw)

    # Test 3: Column with normal cos_zenith (0.5) should have positive fluxes
    @test all(flux_up_sw[:, 1] .> FT(0))
    @test all(flux_dn_sw[:, 1] .> FT(0))
    @test flux_dn_dir_sw[1, 1] > FT(0)

    # Test 4: Columns with cos_zenith ≤ 0 should have zero shortwave fluxes
    # (columns 2 and 4: cos_zenith = 0 and cos_zenith = -0.5)
    @test all(flux_up_sw[:, 2] .== FT(0))
    @test all(flux_dn_sw[:, 2] .== FT(0))
    @test all(flux_up_sw[:, 4] .== FT(0))
    @test all(flux_dn_sw[:, 4] .== FT(0))

    # Test 5: Column with very small cos_zenith (1e-10) should still produce valid results
    # The flux should be very small but finite
    @test all(isfinite.(flux_up_sw[:, 3]))
    @test all(isfinite.(flux_dn_sw[:, 3]))
    @test isfinite(flux_dn_dir_sw[1, 3])

    # Test 6: Aerosol optical depth should be non-negative and physically consistent
    # Check no NaN or Inf in aerosol optical depths
    @test all(isfinite, as.aerosol_state.aod_sw_ext)
    @test all(isfinite, as.aerosol_state.aod_sw_sca)
    @test all(as.aerosol_state.aod_sw_ext .> FT(0))
    @test all(as.aerosol_state.aod_sw_sca .> FT(0))

    # Test 7: Verify all columns have the same AOD values (since we copied data from column 1 to all columns)
    # Note: Use Array() to avoid GPU scalar indexing when extracting the first element
    aod_ext_first = Array(as.aerosol_state.aod_sw_ext)[1]
    aod_sca_first = Array(as.aerosol_state.aod_sw_sca)[1]
    @test all(Array(as.aerosol_state.aod_sw_ext) .== aod_ext_first)
    @test all(Array(as.aerosol_state.aod_sw_sca) .== aod_sca_first)

    return nothing
end

# Run tests if executed directly
context = ClimaComms.context()

@testset "cos_zenith edge cases with TwoStream SW solver (Float64)" begin
    @time test_cos_zenith_edge_cases(context, NoScatLWRTE, TwoStreamSWRTE, Float64; ncol = 4)
end

@testset "cos_zenith edge cases with TwoStream SW solver (Float32)" begin
    @time test_cos_zenith_edge_cases(context, NoScatLWRTE, TwoStreamSWRTE, Float32; ncol = 4)
end

