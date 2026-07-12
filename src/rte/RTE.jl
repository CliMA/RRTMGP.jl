module RTE
using Adapt
import ClimaComms
using ..AngularDiscretizations
using ..AtmosphericStates
using ..Sources
using ..Fluxes
using ..Optics
import ..RRTMGPGridParams
using ..BCs

import ..Parameters as RP

export NoScatLWRTE, TwoStreamLWRTE, NoScatSWRTE, TwoStreamSWRTE

"""
    NoScatLWRTE(grid_params::RRTMGPGridParams; params, sfc_emis, inc_flux)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a non-scattering longwave simulation.

# Fields
- `context`: ClimaComms context.
- `op`: Optical properties.
- `src`: Longwave source functions.
- `bcs`: Longwave boundary conditions.
- `fluxb`: Temporary storage for bandwise calculations.
- `flux`: Longwave fluxes.
- `angle_disc`: Angular discretization.
- `state_cache`: [`TransposedStateCache`](@ref
  RRTMGP.AtmosphericStates.TransposedStateCache) with column-first copies of
  the hot state arrays (refreshed at each spectral solve), or `nothing`.
"""
struct NoScatLWRTE{
    C,
    OP <: OneScalar,
    SL <: SourceLWNoScat,
    BC <: LwBCs,
    FXBL,
    FXL <: FluxLW,
    AD,
    TSC,
}
    context::C
    op::OP
    src::SL
    bcs::BC
    fluxb::FXBL
    flux::FXL
    angle_disc::AD
    state_cache::TSC
end
Adapt.@adapt_structure NoScatLWRTE

function NoScatLWRTE(
    grid_params::RRTMGPGridParams;
    params,
    sfc_emis,
    inc_flux,
    state_cache = TransposedStateCache(grid_params),
)
    (; context) = grid_params
    op = OneScalar(grid_params)
    src = SourceLWNoScat(grid_params; params)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(grid_params)
    flux = FluxLW(grid_params)
    ad = AngularDiscretization(grid_params, 1)
    return NoScatLWRTE(context, op, src, bcs, fluxb, flux, ad, state_cache)
end

"""
    TwoStreamLWRTE(grid_params::RRTMGPGridParams; params, sfc_emis, inc_flux)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a `2-stream` longwave simulation.

# Fields
- `context`: ClimaComms context.
- `op`: Optical properties.
- `src`: Longwave source functions.
- `bcs`: Longwave boundary conditions.
- `fluxb`: Temporary storage for bandwise calculations.
- `flux`: Longwave fluxes.
- `band_flux`: Optional per-band longwave fluxes `(nlev, ncol, n_bnd)`, or `nothing`.
- `state_cache`: [`TransposedStateCache`](@ref
  RRTMGP.AtmosphericStates.TransposedStateCache) with column-first copies of
  the hot state arrays (refreshed at each spectral solve), or `nothing`.
"""
struct TwoStreamLWRTE{
    C,
    OP <: TwoStream,
    SL <: SourceLW2Str,
    BC <: LwBCs,
    FXBL,
    FXL <: FluxLW,
    FXBND,
    TSC,
}
    context::C
    op::OP
    src::SL
    bcs::BC
    fluxb::FXBL
    flux::FXL
    band_flux::FXBND
    state_cache::TSC
end
Adapt.@adapt_structure TwoStreamLWRTE
# By default no per-band (spectrally-resolved) fluxes are retained.
TwoStreamLWRTE(context, op, src, bcs, fluxb, flux, band_flux) =
    TwoStreamLWRTE(context, op, src, bcs, fluxb, flux, band_flux, nothing)

function TwoStreamLWRTE(
    grid_params::RRTMGPGridParams;
    params,
    sfc_emis,
    inc_flux,
    state_cache = TransposedStateCache(grid_params),
)
    (; context) = grid_params
    op = TwoStream(grid_params)
    src = SourceLW2Str(grid_params; params)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(grid_params)
    flux = FluxLW(grid_params)
    return TwoStreamLWRTE(
        context,
        op,
        src,
        bcs,
        fluxb,
        flux,
        nothing,
        state_cache,
    )
end

"""
    NoScatSWRTE(grid_params::RRTMGPGridParams; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a non-scattering shortwave simulation.

# Fields
- `context`: ClimaComms context.
- `op`: Optical properties.
- `bcs`: Shortwave boundary conditions.
- `fluxb`: Temporary storage for bandwise calculations.
- `flux`: Shortwave fluxes.
- `state_cache`: [`TransposedStateCache`](@ref
  RRTMGP.AtmosphericStates.TransposedStateCache) with column-first copies of
  the hot state arrays (refreshed at each spectral solve), or `nothing`.
"""
struct NoScatSWRTE{C, OP <: OneScalar, BC <: SwBCs, FXBS, FXS <: FluxSW, TSC}
    context::C
    op::OP
    bcs::BC
    fluxb::FXBS
    flux::FXS
    state_cache::TSC
end
Adapt.@adapt_structure NoScatSWRTE

function NoScatSWRTE(
    grid_params::RRTMGPGridParams;
    cos_zenith,
    toa_flux,
    sfc_alb_direct,
    inc_flux_diffuse,
    sfc_alb_diffuse,
    state_cache = TransposedStateCache(grid_params),
)
    (; context) = grid_params
    op = OneScalar(grid_params)
    bcs = SwBCs(
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
    fluxb = FluxSW(grid_params)
    flux = FluxSW(grid_params)
    return NoScatSWRTE(context, op, bcs, fluxb, flux, state_cache)
end

"""
    TwoStreamSWRTE(grid_params::RRTMGPGridParams; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a `2-stream` shortwave simulation.

# Fields
- `context`: ClimaComms context.
- `op`: Optical properties.
- `src`: Shortwave source functions.
- `bcs`: Shortwave boundary conditions.
- `fluxb`: Temporary storage for bandwise calculations.
- `flux`: Shortwave fluxes.
- `band_flux`: Optional per-band shortwave fluxes `(nlev, ncol, n_bnd)`, or `nothing`.
- `state_cache`: [`TransposedStateCache`](@ref
  RRTMGP.AtmosphericStates.TransposedStateCache) with column-first copies of
  the hot state arrays (refreshed at each spectral solve), or `nothing`.
"""
struct TwoStreamSWRTE{
    C,
    OP <: TwoStream,
    SS,
    BC <: SwBCs,
    FXBS,
    FXS <: FluxSW,
    FXBND,
    TSC,
}
    context::C
    op::OP
    src::SS
    bcs::BC
    fluxb::FXBS
    flux::FXS
    band_flux::FXBND
    state_cache::TSC
end
Adapt.@adapt_structure TwoStreamSWRTE
# By default no per-band (spectrally-resolved) fluxes are retained.
TwoStreamSWRTE(context, op, src, bcs, fluxb, flux, band_flux) =
    TwoStreamSWRTE(context, op, src, bcs, fluxb, flux, band_flux, nothing)

function TwoStreamSWRTE(
    grid_params::RRTMGPGridParams;
    cos_zenith,
    toa_flux,
    sfc_alb_direct,
    inc_flux_diffuse,
    sfc_alb_diffuse,
    state_cache = TransposedStateCache(grid_params),
)
    (; context) = grid_params
    op = TwoStream(grid_params)
    src = SourceSW2Str(grid_params)
    bcs = SwBCs(
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
    fluxb = FluxSW(grid_params)
    flux = FluxSW(grid_params)
    return TwoStreamSWRTE(
        context,
        op,
        src,
        bcs,
        fluxb,
        flux,
        nothing,
        state_cache,
    )
end


end
