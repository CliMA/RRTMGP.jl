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
    NoScatLWRTE(
        ::Type{FT},
        ::Type{DA},
        ::Type{OP},
        context,
        param_set,
        nlay,
        ncol,
        sfc_emis,
        inc_flux,
    )

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a non-scattering longwave simulation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct NoScatLWRTE{C, OP <: OneScalar, SL <: SourceLWNoScat, BC <: LwBCs, FXBL, FXL <: FluxLW, AD}
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

function NoScatLWRTE(::Type{FT}, ::Type{DA}, context, params, nlay, ncol, sfc_emis, inc_flux) where {FT, DA}
    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)
    @warn "Please call NoScatLWRTE with RRTMGPGridParams instead"
    return NoScatLWRTE(grid_params; params, sfc_emis, inc_flux)
end

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
    TwoStreamLWRTE(::Type{FT}, ::Type{DA}, context, param_set, nlay, ncol, sfc_emis, inc_flux)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a `2-stream` longwave simulation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamLWRTE{C, OP <: TwoStream, SL <: SourceLW2Str, BC <: LwBCs, FXBL, FXL <: FluxLW}
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
end
Adapt.@adapt_structure TwoStreamLWRTE

function TwoStreamLWRTE(::Type{FT}, ::Type{DA}, context, params, nlay, ncol, sfc_emis, inc_flux) where {FT, DA}
    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)
    @warn "Please call TwoStreamLWRTE with RRTMGPGridParams instead."
    return TwoStreamLWRTE(grid_params; params, sfc_emis, inc_flux)
end

function TwoStreamLWRTE(grid_params::RRTMGPGridParams; params, sfc_emis, inc_flux)
    (; context) = grid_params
    op = TwoStream(grid_params)
    src = SourceLW2Str(grid_params; params)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(grid_params)
    flux = FluxLW(grid_params)
    return TwoStreamLWRTE(context, op, src, bcs, fluxb, flux)
end

"""
    NoScatSWRTE(::Type{FT}, ::Type{DA}, context, nlay, ncol, swbcs...)

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

function NoScatSWRTE(::Type{FT}, ::Type{DA}, context, nlay, ncol, swbcs...) where {FT, DA}
    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)
    @warn "Please call NoScatSWRTE with RRTMGPGridParams instead"
    (cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse) = swbcs
    return NoScatSWRTE(grid_params; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
end

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
    bcs = SwBCs(cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    fluxb = FluxSW(grid_params)
    flux = FluxSW(grid_params)
    return NoScatSWRTE(context, op, bcs, fluxb, flux)
end

"""
    TwoStreamSWRTE(::Type{FT}, ::Type{DA}, context, nlay, ncol, swbcs...)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a `2-stream` shortwave simulation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamSWRTE{C, OP <: TwoStream, SS, BC <: SwBCs, FXBS, FXS <: FluxSW}
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
end
Adapt.@adapt_structure TwoStreamSWRTE

function TwoStreamSWRTE(::Type{FT}, ::Type{DA}, context, nlay, ncol, swbcs...) where {FT, DA}
    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)
    (cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse) = swbcs
    @warn "Please call TwoStreamSWRTE with RRTMGPGridParams instead."
    return TwoStreamSWRTE(grid_params; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
end

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
    bcs = SwBCs(cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    fluxb = FluxSW(grid_params)
    flux = FluxSW(grid_params)
    return TwoStreamSWRTE(context, op, src, bcs, fluxb, flux)
end


end
