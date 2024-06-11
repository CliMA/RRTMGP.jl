module RTE
using Adapt
import ClimaComms
using ..AtmosphericStates
using DocStringExtensions
using ..Sources
using ..Fluxes
using ..Optics
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
struct NoScatLWRTE{C, OP, SL <: SourceLWNoScat, BC <: LwBCs, FXBL, FXL <: FluxLW}
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
Adapt.@adapt_structure NoScatLWRTE

function NoScatLWRTE(
    ::Type{FT},
    ::Type{DA},
    ::Type{OP},
    context,
    param_set,
    nlay,
    ncol,
    sfc_emis,
    inc_flux,
) where {FT, DA, OP}
    op = OP(FT, ncol, nlay, DA)
    src = SourceLWNoScat(param_set, FT, DA, nlay, ncol)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(ncol, nlay, FT, DA)
    flux = FluxLW(ncol, nlay, FT, DA)
    return NoScatLWRTE(context, op, src, bcs, fluxb, flux)
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

function TwoStreamLWRTE(::Type{FT}, ::Type{DA}, context, param_set, nlay, ncol, sfc_emis, inc_flux) where {FT, DA}
    op = TwoStream(FT, ncol, nlay, DA)
    src = SourceLW2Str(param_set, FT, DA, nlay, ncol)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(ncol, nlay, FT, DA)
    flux = FluxLW(ncol, nlay, FT, DA)
    return TwoStreamLWRTE(context, op, src, bcs, fluxb, flux)
end

TwoStreamLWRTE(::Type{FT}, ::Type{DA}, ::Type{OP}, args...) where {FT, DA, OP} = TwoStreamLWRTE(FT, DA, args...)

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
    op = OneScalar(FT, ncol, nlay, DA)
    bcs = SwBCs(swbcs...)
    fluxb = FluxSW(ncol, nlay, FT, DA)
    flux = FluxSW(ncol, nlay, FT, DA)
    return NoScatSWRTE(context, op, bcs, fluxb, flux)
end

NoScatSWRTE(::Type{FT}, ::Type{DA}, ::Type{OP}, args...) where {FT, DA, OP} = NoScatSWRTE(FT, DA, args...)

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
    op = TwoStream(FT, ncol, nlay, DA)
    src = SourceSW2Str(FT, DA, nlay, ncol)
    bcs = SwBCs(swbcs...)
    fluxb = FluxSW(ncol, nlay, FT, DA)
    flux = FluxSW(ncol, nlay, FT, DA)
    return TwoStreamSWRTE(context, op, src, bcs, fluxb, flux)
end

TwoStreamSWRTE(::Type{FT}, ::Type{DA}, ::Type{OP}, args...) where {FT, DA, OP} = TwoStreamSWRTE(FT, DA, args...)
end
