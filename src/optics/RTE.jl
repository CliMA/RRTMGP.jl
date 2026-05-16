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

using ..Canopy

export NoScatLWRTE, TwoStreamLWRTE, NoScatSWRTE, TwoStreamSWRTE,
    TwoStreamSWCanopyRTE, TwoStreamLWCanopyRTE

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


"""
    TwoStreamSWCanopyRTE

A high-level RRTMGP data structure for `2-stream` shortwave simulation
with coupled atmosphere-canopy radiative transfer.

Carries two sets of optical properties: `op_atm` (atmospheric dimensions for
`compute_optical_props!`) and `op_exp` (expanded column = atmosphere + canopy).

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamSWCanopyRTE{C, OPA <: TwoStream, OPE <: TwoStream, SSA, SSE, BC <: SwBCs, FXBS, FXS <: FluxSW, CS <: CanopyState}
    "ClimaComms context"
    context::C
    "atmospheric optical properties (nlay_atm layers)"
    op_atm::OPA
    "expanded column optical properties (nlay_atm + ncanlay layers)"
    op_exp::OPE
    "atmospheric shortwave source workspace"
    src_atm::SSA
    "expanded column shortwave source"
    src_exp::SSE
    "shortwave boundary conditions (soil albedos)"
    bcs::BC
    "temporary storage for per-gpoint fluxes (expanded dims)"
    fluxb::FXBS
    "accumulated shortwave fluxes (expanded dims)"
    flux::FXS
    "canopy state"
    canopy::CS
end
Adapt.@adapt_structure TwoStreamSWCanopyRTE

function TwoStreamSWCanopyRTE(
    grid_params::RRTMGPGridParams;
    canopy::CanopyState,
    cos_zenith,
    toa_flux,
    inc_flux_diffuse = nothing,
)
    (; context) = grid_params
    FT = eltype(grid_params)
    nlay_atm = grid_params.nlay
    ncol = grid_params.ncol
    ncanlay = canopy.ncanlay

    # Atmospheric-dimension workspace
    op_atm = TwoStream(grid_params)
    src_atm = SourceSW2Str(grid_params)

    # Expanded-dimension workspace (atm + canopy)
    grid_exp = RRTMGPGridParams(FT; context, nlay = nlay_atm + ncanlay, ncol)
    op_exp = TwoStream(grid_exp)
    src_exp = SourceSW2Str(grid_exp)
    fluxb = FluxSW(grid_exp)
    flux = FluxSW(grid_exp)

    # BCs use soil albedos from canopy state
    bcs = SwBCs(cos_zenith, toa_flux, canopy.soil_alb_direct, inc_flux_diffuse, canopy.soil_alb_diffuse)

    return TwoStreamSWCanopyRTE(context, op_atm, op_exp, src_atm, src_exp, bcs, fluxb, flux, canopy)
end


"""
    TwoStreamLWCanopyRTE

A high-level RRTMGP data structure for `2-stream` longwave simulation
with coupled atmosphere-canopy radiative transfer.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamLWCanopyRTE{C, OPA <: TwoStream, OPE <: TwoStream, SLA, SLE, BC <: LwBCs, FXBL, FXL <: FluxLW, CS <: CanopyState}
    "ClimaComms context"
    context::C
    "atmospheric optical properties (nlay_atm layers)"
    op_atm::OPA
    "expanded column optical properties (nlay_atm + ncanlay layers)"
    op_exp::OPE
    "atmospheric longwave source"
    src_atm::SLA
    "expanded column longwave source"
    src_exp::SLE
    "longwave boundary conditions (soil emissivity)"
    bcs::BC
    "temporary storage for per-gpoint fluxes (expanded dims)"
    fluxb::FXBL
    "accumulated longwave fluxes (expanded dims)"
    flux::FXL
    "canopy state"
    canopy::CS
end
Adapt.@adapt_structure TwoStreamLWCanopyRTE

function TwoStreamLWCanopyRTE(
    grid_params::RRTMGPGridParams;
    canopy::CanopyState,
    params,
)
    (; context) = grid_params
    FT = eltype(grid_params)
    nlay_atm = grid_params.nlay
    ncol = grid_params.ncol
    ncanlay = canopy.ncanlay

    # Atmospheric-dimension workspace
    op_atm = TwoStream(grid_params)
    src_atm = SourceLW2Str(grid_params; params)

    # Expanded-dimension workspace
    grid_exp = RRTMGPGridParams(FT; context, nlay = nlay_atm + ncanlay, ncol)
    op_exp = TwoStream(grid_exp)
    src_exp = SourceLW2Str(grid_exp; params)
    fluxb = FluxLW(grid_exp)
    flux = FluxLW(grid_exp)

    # BCs use soil emissivity from canopy state
    bcs = LwBCs(canopy.soil_emis, nothing)

    return TwoStreamLWCanopyRTE(context, op_atm, op_exp, src_atm, src_exp, bcs, fluxb, flux, canopy)
end


end
