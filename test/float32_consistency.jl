#=
Float32 ↔ Float64 consistency harness.

Solves the same radiation problems at Float32 and Float64 (states built from the
same NetCDF inputs, downcast per precision) and compares broadband fluxes. This
pins the end-to-end Float32 accuracy of the pipeline — input rounding, gas/cloud
optics, and the RTE solve — directly against the Float64 result, which is a much
sharper probe than the per-precision reference comparisons.

The thresholds below are a RATCHET: they are set from measured values (with ~2×
headroom) and must be tightened, never loosened, as numerics fixes land. For
context, the Fortran reference (rte-rrtmgp) accepts 0.35 W/m² absolute flux error
in its single-precision CI.
=#

using Test
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
using RRTMGP.AngularDiscretizations
using RRTMGP.RTE
using RRTMGP.RTESolver
import RRTMGP.Parameters.RRTMGPParameters
import ClimaParams as CP
using RRTMGP.ArtifactPaths

@isdefined(get_reference_filename) || include("reference_files.jl")
@isdefined(setup_clear_sky_as) || include("read_clear_sky.jl")
@isdefined(setup_cloudy_sky_as) || include("read_cloudy_sky.jl")

# Maximum |Float32 − Float64| flux thresholds [W/m²] (gray heating rate [K/s]).
# Measured-with-headroom; tighten as numerics fixes land, never loosen.
# Measured after the f32 numerics program (τ_thresh at eps^(1/4), the γ1−γ2
# identities, expm1 thin-layer factors, the SW resonance nudge, exact
# delta-scale forms, the addition-built direct-beam prepass, and the
# restructured LW two-stream source): every LW path sits at its ~2–5e-4
# interpolation-noise floor, 25–90× below the pre-program values (LW noscat was
# 1.1e-2 clear / 3.3e-2 cloudy; LW 2-stream 3.2e-2 clear / 3.4e-2 cloudy).
# SW: ~1.1e-2/1.6e-2 clear, ~4.9e-2/5.5e-2 cloudy — verified to be per-g-point
# interpolation/coefficient noise, NOT accumulation: repeating the runs with
# full Float64 broadband accumulation left the SW differences unchanged, so
# compensated/f64 accumulators were deliberately not added.
const F32_THRESHOLDS = (;
    gray_flux = 1.0e-3,
    gray_heating_rate = 1.0e-8,
    clear_lw_noscat = 1.0e-3,
    clear_lw_2stream = 1.0e-3,
    clear_sw = 3.0e-2,
    cloudy_lw_noscat = 1.0e-3,
    cloudy_lw_2stream = 1.0e-3,
    cloudy_sw = 1.2e-1,
)

maxabsdiff(a32, a64) = maximum(abs.(Float64.(Array(a32)) .- Array(a64)))

#------------------------------------------------------------------
# Gray atmosphere (two-stream LW + SW): no lookup tables involved.
function gray_consistency(; nlay = 60, ncol = 8)
    out = map((Float32, Float64)) do FT
        RRTMGP.solve_gray(FT; nlay, ncol)
    end
    o32, o64 = out
    return (;
        lw_net = maxabsdiff(o32.lw_net, o64.lw_net),
        sw_net = maxabsdiff(o32.sw_net, o64.sw_net),
        heating_rate = maxabsdiff(o32.heating_rate, o64.heating_rate),
    )
end

#------------------------------------------------------------------
# Clear-sky (RFMIP state): gas optics + LW (noscat and 2-stream) + SW 2-stream.
function clear_sky_fluxes(context, ::Type{FT}, ::Type{SLVLW}; ncol) where {FT, SLVLW}
    overrides =
        (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)

    lookup_lw, idx_gases = Dataset(get_lookup_filename(:gas, :lw), "r") do ds
        LookUpLW(ds, FT, DA)
    end
    lookup_sw, _ = Dataset(get_lookup_filename(:gas, :sw), "r") do ds
        LookUpSW(ds, FT, DA)
    end

    ds_lw_in = Dataset(get_input_filename(:gas, :lw), "r")
    (as, sfc_emis, sfc_alb_direct, cos_zenith, toa_flux, _) =
        setup_clear_sky_as(
            context,
            ds_lw_in,
            idx_gases,
            1,
            lookup_lw,
            ncol,
            FT,
            VmrGM,
            param_set,
        )
    close(ds_lw_in)

    nlay, _ = AtmosphericStates.get_dims(as)
    grid_params = RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)

    slv_lw = SLVLW(grid_params; params = param_set, sfc_emis, inc_flux = nothing)
    sfc_alb_diffuse = DA{FT, 2}(deepcopy(sfc_alb_direct))
    slv_sw = TwoStreamSWRTE(
        grid_params;
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse = nothing,
        sfc_alb_diffuse,
    )

    solve_lw!(slv_lw, as, lookup_lw)
    solve_sw!(slv_sw, as, lookup_sw)
    return (;
        lw_up = Array(slv_lw.flux.flux_up),
        lw_dn = Array(slv_lw.flux.flux_dn),
        sw_up = Array(slv_sw.flux.flux_up),
        sw_dn = Array(slv_sw.flux.flux_dn),
    )
end

#------------------------------------------------------------------
# Cloudy sky (cldfrac = 1 → deterministic McICA mask): cloud optics + δ-scaling.
function cloudy_sky_fluxes(context, ::Type{FT}, ::Type{SLVLW}; ncol) where {FT, SLVLW}
    overrides =
        (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)

    lookup_lw, idx_gases = Dataset(get_lookup_filename(:gas, :lw), "r") do ds
        LookUpLW(ds, FT, DA)
    end
    lookup_lw_cld = Dataset(get_lookup_filename(:cloud, :lw), "r") do ds
        LookUpCld(ds, FT, DA)
    end
    lookup_sw, _ = Dataset(get_lookup_filename(:gas, :sw), "r") do ds
        LookUpSW(ds, FT, DA)
    end
    lookup_sw_cld = Dataset(get_lookup_filename(:cloud, :sw), "r") do ds
        LookUpCld(ds, FT, DA)
    end

    ds_in = Dataset(get_input_filename(:gas_clouds, :lw), "r")
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, cos_zenith, toa_flux, _ =
        setup_cloudy_sky_as(
            context,
            ds_in,
            idx_gases,
            lookup_lw,
            lookup_sw,
            lookup_lw_cld,
            lookup_sw_cld,
            FT(1), # cldfrac = 1 everywhere the input has cloud: deterministic mask
            ncol,
            FT,
            param_set,
        )
    close(ds_in)

    nlay, ncol = AtmosphericStates.get_dims(as)
    grid_params = RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)

    slv_lw = SLVLW(grid_params; params = param_set, sfc_emis, inc_flux = nothing)
    slv_sw = TwoStreamSWRTE(
        grid_params;
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse = nothing,
        sfc_alb_diffuse,
    )

    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld)
    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld)
    return (;
        lw_up = Array(slv_lw.flux.flux_up),
        lw_dn = Array(slv_lw.flux.flux_dn),
        sw_up = Array(slv_sw.flux.flux_up),
        sw_dn = Array(slv_sw.flux.flux_dn),
    )
end

fluxdiffs(f32, f64) = (;
    lw_up = maxabsdiff(f32.lw_up, f64.lw_up),
    lw_dn = maxabsdiff(f32.lw_dn, f64.lw_dn),
    sw_up = maxabsdiff(f32.sw_up, f64.sw_up),
    sw_dn = maxabsdiff(f32.sw_dn, f64.sw_dn),
)
lwmax(d) = max(d.lw_up, d.lw_dn)
swmax(d) = max(d.sw_up, d.sw_dn)

@testset "Float32 ↔ Float64 consistency" begin
    context = ClimaComms.context()

    gray = gray_consistency()
    @info "f32↔f64 gray (2-stream)" gray.lw_net gray.sw_net gray.heating_rate
    @test max(gray.lw_net, gray.sw_net) ≤ F32_THRESHOLDS.gray_flux
    @test gray.heating_rate ≤ F32_THRESHOLDS.gray_heating_rate

    for (SLVLW, tag, thresh) in (
        (NoScatLWRTE, "noscat", F32_THRESHOLDS.clear_lw_noscat),
        (TwoStreamLWRTE, "2stream", F32_THRESHOLDS.clear_lw_2stream),
    )
        d = fluxdiffs(
            clear_sky_fluxes(context, Float32, SLVLW; ncol = 100),
            clear_sky_fluxes(context, Float64, SLVLW; ncol = 100),
        )
        @info "f32↔f64 clear-sky (LW $tag)" d.lw_up d.lw_dn d.sw_up d.sw_dn
        @test lwmax(d) ≤ thresh
        @test swmax(d) ≤ F32_THRESHOLDS.clear_sw
    end

    for (SLVLW, tag, thresh) in (
        (NoScatLWRTE, "noscat", F32_THRESHOLDS.cloudy_lw_noscat),
        (TwoStreamLWRTE, "2stream", F32_THRESHOLDS.cloudy_lw_2stream),
    )
        d = fluxdiffs(
            cloudy_sky_fluxes(context, Float32, SLVLW; ncol = 128),
            cloudy_sky_fluxes(context, Float64, SLVLW; ncol = 128),
        )
        @info "f32↔f64 cloudy-sky (LW $tag)" d.lw_up d.lw_dn d.sw_up d.sw_dn
        @test lwmax(d) ≤ thresh
        @test swmax(d) ≤ F32_THRESHOLDS.cloudy_sw
    end
end
