module RTE
using Adapt
import ClimaComms
using ..AngularDiscretizations
using ..AtmosphericStates
using DocStringExtensions
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
$(DocStringExtensions.FIELDS)
"""
struct NoScatLWRTE{
    C,
    OP <: OneScalar,
    SL <: SourceLWNoScat,
    BC <: LwBCs,
    FXBL,
    FXL <: FluxLW,
    AD,
}
    "ClimaComms context"
    context::C
    "optical properties"
    op::OP
    "longwave source functions"
    src::SL
    "longwave boundary conditions"
    bcs::BC
    "temporary storage for bandwise calculations"
    fluxb::FXBL
    "longwave fluxes"
    flux::FXL
    "Angular discretization"
    angle_disc::AD
end
Adapt.@adapt_structure NoScatLWRTE

function NoScatLWRTE(grid_params::RRTMGPGridParams; params, sfc_emis, inc_flux)
    (; context) = grid_params
    op = OneScalar(grid_params)
    src = SourceLWNoScat(grid_params; params)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(grid_params)
    flux = FluxLW(grid_params)
    ad = AngularDiscretization(grid_params, 1)
    return NoScatLWRTE(context, op, src, bcs, fluxb, flux, ad)
end

"""
    TwoStreamLWRTE(grid_params::RRTMGPGridParams; params, sfc_emis, inc_flux)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a `2-stream` longwave simulation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamLWRTE{
    C,
    OP <: TwoStream,
    SL <: SourceLW2Str,
    BC <: LwBCs,
    FXBL,
    FXL <: FluxLW,
    FXBND,
}
    "ClimaComms context"
    context::C
    "optical properties"
    op::OP
    "longwave source functions"
    src::SL
    "longwave boundary conditions"
    bcs::BC
    "temporary storage for bandwise calculations"
    fluxb::FXBL
    "longwave fluxes"
    flux::FXL
    "optional per-band longwave fluxes `(nlev, ncol, n_bnd)`, or `nothing`"
    band_flux::FXBND
end
Adapt.@adapt_structure TwoStreamLWRTE
# By default no per-band (spectrally-resolved) fluxes are retained.
TwoStreamLWRTE(context, op, src, bcs, fluxb, flux) =
    TwoStreamLWRTE(context, op, src, bcs, fluxb, flux, nothing)

function TwoStreamLWRTE(
    grid_params::RRTMGPGridParams;
    params,
    sfc_emis,
    inc_flux,
)
    (; context) = grid_params
    op = TwoStream(grid_params)
    src = SourceLW2Str(grid_params; params)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(grid_params)
    flux = FluxLW(grid_params)
    return TwoStreamLWRTE(context, op, src, bcs, fluxb, flux)
end

"""
    NoScatSWRTE(grid_params::RRTMGPGridParams; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a non-scattering shortwave simulation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct NoScatSWRTE{C, OP <: OneScalar, BC <: SwBCs, FXBS, FXS <: FluxSW}
    "ClimaComms context"
    context::C
    "optical properties"
    op::OP
    "shortwave boundary conditions"
    bcs::BC
    "temporary storage for bandwise calculations"
    fluxb::FXBS
    "shortwave fluxes"
    flux::FXS
end
Adapt.@adapt_structure NoScatSWRTE

function NoScatSWRTE(
    grid_params::RRTMGPGridParams;
    cos_zenith,
    toa_flux,
    sfc_alb_direct,
    inc_flux_diffuse,
    sfc_alb_diffuse,
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
    return NoScatSWRTE(context, op, bcs, fluxb, flux)
end

"""
    TwoStreamSWRTE(grid_params::RRTMGPGridParams; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a `2-stream` shortwave simulation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamSWRTE{
    C,
    OP <: TwoStream,
    SS,
    BC <: SwBCs,
    FXBS,
    FXS <: FluxSW,
    FXBND,
}
    "ClimaComms context"
    context::C
    "optical properties"
    op::OP
    "shortwave source functions"
    src::SS
    "shortwave boundary conditions"
    bcs::BC
    "temporary storage for bandwise calculations"
    fluxb::FXBS
    "shortwave fluxes"
    flux::FXS
    "optional per-band shortwave fluxes `(nlev, ncol, n_bnd)`, or `nothing`"
    band_flux::FXBND
end
Adapt.@adapt_structure TwoStreamSWRTE
# By default no per-band (spectrally-resolved) fluxes are retained.
TwoStreamSWRTE(context, op, src, bcs, fluxb, flux) =
    TwoStreamSWRTE(context, op, src, bcs, fluxb, flux, nothing)

function TwoStreamSWRTE(
    grid_params::RRTMGPGridParams;
    cos_zenith,
    toa_flux,
    sfc_alb_direct,
    inc_flux_diffuse,
    sfc_alb_diffuse,
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
    return TwoStreamSWRTE(context, op, src, bcs, fluxb, flux)
end


end
