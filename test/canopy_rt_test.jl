#=
End-to-end test for coupled atmosphere-canopy radiative transfer.
Validates energy conservation for SW and LW using the new package API.
=#

using NCDatasets
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

using RRTMGP
using RRTMGP: RRTMGPGridParams
using RRTMGP.Vmrs
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.RTE
using RRTMGP.RTESolver
using RRTMGP.ArtifactPaths
using RRTMGP.Canopy
using RRTMGP.Canopy: fill_leaf_optics!, fill_G!, fill_beta_dir!
import RRTMGP.Parameters.RRTMGPParameters
import ClimaParams as CP

using CanopyOptics
using CanopyOptics: G, spherical_leaves
using Unitful

using Test

include(joinpath(@__DIR__, "reference_files.jl"))
include(joinpath(@__DIR__, "read_clear_sky.jl"))

@testset "Canopy RT energy conservation" begin
    FT = Float64
    ncol = 1
    ncanlay = 5
    LAI_total = FT(3)
    soil_alb_val = FT(0.2)
    sza_deg = FT(30)
    ε_leaf = FT(0.95)
    soil_emis_val = FT(0.96)

    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)

    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)

    # Load SW lookup
    sw_file = get_lookup_filename(:gas, :sw)
    lookup_sw, idx_gases = Dataset(sw_file, "r") do ds
        LookUpSW(ds, FT, DA)
    end
    nbnd_sw = LookUpTables.get_n_bnd(lookup_sw)
    bnd_lims_wn = Array(lookup_sw.band_data.bnd_lims_wn)

    # Load LW lookup
    lw_file = get_lookup_filename(:gas, :lw)
    lookup_lw, _ = Dataset(lw_file, "r") do ds
        LookUpLW(ds, FT, DA)
    end
    nbnd_lw = LookUpTables.get_n_bnd(lookup_lw)

    # Setup atmospheric state
    input_file = get_input_filename(:gas, :lw)
    ds_in = Dataset(input_file, "r")
    (as, _, _, _, toa_flux_orig, _) =
        setup_clear_sky_as(context, ds_in, idx_gases, 1, lookup_sw, ncol, FT, Vmrs.VmrGM, param_set)
    close(ds_in)

    nlay_atm, _ = AtmosphericStates.get_dims(as)
    grid_atm = RRTMGPGridParams(FT; context, nlay = nlay_atm, ncol)

    cos_zenith = DA{FT, 1}([cosd(sza_deg)])
    toa_flux = DA{FT, 1}(Array(toa_flux_orig)[1:ncol])

    # --- Create CanopyState and fill via extension ---
    canopy = CanopyState(grid_atm; ncanlay, nbnd_sw, nbnd_lw)

    # Fill LAI
    LAI_layer = FT(LAI_total / ncanlay)
    for ilay in 1:ncanlay, gcol in 1:ncol
        canopy.lai_layer[ilay, gcol] = LAI_layer
    end

    # Leaf optics via PROSPECT
    leaf = CanopyOptics.LeafProspectProProperties{FT}(
        N = FT(1.4), Ccab = FT(40), Ccar = FT(10), Canth = FT(0.5),
        Cbrown = FT(0), Cw = FT(0.012), Cm = FT(0), Cprot = FT(0.001), Ccbc = FT(0.009),
    )
    opti = createLeafOpticalStruct((400.0:2500.0) * u"nm")
    fill_leaf_optics!(canopy, leaf, opti, bnd_lims_wn)

    # G function and beta_dir
    ld = spherical_leaves()
    fill_G!(canopy, ld, Array(cos_zenith))
    fill_beta_dir!(canopy, ld, Array(cos_zenith))

    # LW leaf optics (spectrally flat)
    ω_lw = FT(1) - ε_leaf
    for ibnd in 1:nbnd_lw, gcol in 1:ncol
        canopy.leaf_omega_lw[ibnd, gcol] = ω_lw
    end

    # Soil properties
    for ibnd in 1:nbnd_sw, gcol in 1:ncol
        canopy.soil_alb_direct[ibnd, gcol] = soil_alb_val
        canopy.soil_alb_diffuse[ibnd, gcol] = soil_alb_val
    end
    for ibnd in 1:nbnd_lw, gcol in 1:ncol
        canopy.soil_emis[ibnd, gcol] = soil_emis_val
    end

    # Canopy temperature profile (linear from T_soil to T_air)
    T_sfc_val = FT(Array(as.t_sfc)[1])
    T_air_val = FT(Array(as.t_lev)[1, 1])
    for gcol in 1:ncol
        canopy.t_soil[gcol] = T_sfc_val
        for ilev in 1:(ncanlay + 1)
            frac = FT(ilev - 1) / FT(ncanlay)
            canopy.t_canopy_lev[ilev, gcol] = T_sfc_val + frac * (T_air_val - T_sfc_val)
        end
    end

    # --- SW solve ---
    rad_output = CanopyRadiationOutput(grid_atm; ncanlay, nbnd_sw, η_car = FT(0.35))
    sws = TwoStreamSWCanopyRTE(grid_atm; canopy, cos_zenith, toa_flux)
    solve_sw!(sws, as, lookup_sw; rad_output)

    nlev_total = nlay_atm + ncanlay + 1
    flux_up_arr = Array(sws.flux.flux_up)
    flux_dn_arr = Array(sws.flux.flux_dn)

    toa_incoming = flux_dn_arr[nlev_total, 1]
    toa_reflected = flux_up_arr[nlev_total, 1]
    soil_dn = flux_dn_arr[1, 1]
    soil_up = flux_up_arr[1, 1]
    soil_absorbed = soil_dn - soil_up

    # Canopy absorption (flux convergence)
    canopy_absorbed_sw = FT(0)
    for ilay in 1:ncanlay
        flux_in = flux_dn_arr[ilay + 1, 1] + flux_up_arr[ilay, 1]
        flux_out = flux_dn_arr[ilay, 1] + flux_up_arr[ilay + 1, 1]
        canopy_absorbed_sw += flux_in - flux_out
    end

    # Atmospheric absorption
    atm_absorbed_sw = FT(0)
    for ilay in (ncanlay + 1):(nlay_atm + ncanlay)
        flux_in = flux_dn_arr[ilay + 1, 1] + flux_up_arr[ilay, 1]
        flux_out = flux_dn_arr[ilay, 1] + flux_up_arr[ilay + 1, 1]
        atm_absorbed_sw += flux_in - flux_out
    end

    sw_residual = toa_incoming - toa_reflected - atm_absorbed_sw - canopy_absorbed_sw - soil_absorbed

    println("\n--- SW Energy Conservation ---")
    println("  TOA incoming:    $(round(toa_incoming, digits=2))")
    println("  TOA reflected:   $(round(toa_reflected, digits=2))")
    println("  Atm absorbed:    $(round(atm_absorbed_sw, digits=2))")
    println("  Canopy absorbed: $(round(canopy_absorbed_sw, digits=2))")
    println("  Soil absorbed:   $(round(soil_absorbed, digits=2))")
    println("  Residual:        $(round(sw_residual, digits=6))")

    @test abs(sw_residual) < 1e-4
    @test canopy_absorbed_sw > 0  # canopy must absorb some SW
    @test toa_incoming > 0

    # --- LW solve ---
    lws = TwoStreamLWCanopyRTE(grid_atm; canopy, params = param_set)
    solve_lw!(lws, as, lookup_lw)

    flux_up_lw_arr = Array(lws.flux.flux_up)
    flux_dn_lw_arr = Array(lws.flux.flux_dn)

    lw_toa_dn = flux_dn_lw_arr[nlev_total, 1]
    lw_toa_up = flux_up_lw_arr[nlev_total, 1]
    lw_soil_dn = flux_dn_lw_arr[1, 1]
    lw_soil_up = flux_up_lw_arr[1, 1]
    lw_soil_absorbed = lw_soil_dn - lw_soil_up

    canopy_absorbed_lw = FT(0)
    for ilay in 1:ncanlay
        flux_in = flux_dn_lw_arr[ilay + 1, 1] + flux_up_lw_arr[ilay, 1]
        flux_out = flux_dn_lw_arr[ilay, 1] + flux_up_lw_arr[ilay + 1, 1]
        canopy_absorbed_lw += flux_in - flux_out
    end

    atm_absorbed_lw = FT(0)
    for ilay in (ncanlay + 1):(nlay_atm + ncanlay)
        flux_in = flux_dn_lw_arr[ilay + 1, 1] + flux_up_lw_arr[ilay, 1]
        flux_out = flux_dn_lw_arr[ilay, 1] + flux_up_lw_arr[ilay + 1, 1]
        atm_absorbed_lw += flux_in - flux_out
    end

    lw_residual = lw_toa_dn - lw_toa_up - atm_absorbed_lw - canopy_absorbed_lw - lw_soil_absorbed

    println("\n--- LW Energy Conservation ---")
    println("  TOA ↓:           $(round(lw_toa_dn, digits=2))")
    println("  TOA ↑ (OLR):     $(round(lw_toa_up, digits=2))")
    println("  Atm absorbed:    $(round(atm_absorbed_lw, digits=2))")
    println("  Canopy absorbed: $(round(canopy_absorbed_lw, digits=2))")
    println("  Soil absorbed:   $(round(lw_soil_absorbed, digits=2))")
    println("  Residual:        $(round(lw_residual, digits=6))")

    @test abs(lw_residual) < 1e-4
    @test lw_toa_up > 0  # OLR must be positive

    # --- fAPAR and green-APAR diagnostics ---
    # rad_output was already passed to solve_sw! above — apar_by_band is accumulated
    compute_canopy_radiation!(rad_output, sws.flux, lws.flux, canopy, nlay_atm, bnd_lims_wn)

    fAPAR_val = Array(rad_output.fAPAR)[1]
    total_apar = sum(Array(rad_output.apar)[:, 1])
    total_green = sum(Array(rad_output.green_apar)[:, 1])

    println("\n--- fAPAR / Green-APAR ---")
    println("  fAPAR:           $(round(fAPAR_val, digits=4))")
    println("  Total APAR:      $(round(total_apar, digits=2)) W/m²")
    println("  Green-APAR:      $(round(total_green, digits=2)) W/m² (η_car=0.35)")
    println("  f_cab_band[1:5]: $(round.(Array(canopy.f_cab_band)[1:5, 1], digits=3))")
    println("  f_car_band[1:5]: $(round.(Array(canopy.f_car_band)[1:5, 1], digits=3))")

    @test 0 < fAPAR_val < 1
    @test total_apar > 0
    @test total_green ≤ total_apar  # green ≤ total always
    @test total_green ≥ 0
    @test all(Array(canopy.f_cab_band)[:, 1] .≥ 0)
    @test all(Array(canopy.f_cab_band)[:, 1] .+ Array(canopy.f_car_band)[:, 1] .≤ 1 + eps(FT))
end
