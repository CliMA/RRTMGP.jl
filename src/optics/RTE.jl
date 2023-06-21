module RTE
using CUDA
using Adapt
using ..AtmosphericStates
using DocStringExtensions
import GaussQuadrature
using ..Sources
using ..Fluxes
using ..Optics
using ..BCs

import ..Parameters as RP

export Solver

"""
    Solver(
        as,
        op,
        src_lw,
        src_sw,
        bcs_lw,
        bcs_sw,
        fluxb_lw,
        fluxb_sw,
        flux_lw,
        flux_sw
    )

The high-level RRTMGP data structure storing
the atmospheric state, optical properties,
sources, boundary conditions and fluxes
configurations for a given simulation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Solver{AS, OP, SL, SS, BCL, BCS, FXBL, FXBS, FXL, FXS}
    "atmospheric state"
    as::AS
    "optical properties"
    op::OP
    "source functions"
    src_lw::SL
    "source functions"
    src_sw::SS
    "boundary conditions"
    bcs_lw::BCL
    "boundary conditions"
    bcs_sw::BCS
    "temporary storage for bandwise calculations"
    fluxb_lw::FXBL
    "temporary storage for bandwise calculations"
    fluxb_sw::FXBS
    "fluxes for longwave problem"
    flux_lw::FXL
    "fluxes for shortwave problem"
    flux_sw::FXS
    function Solver(as, op, src_lw, src_sw, bcs_lw, bcs_sw, fluxb_lw, fluxb_sw, flux_lw, flux_sw)
        args = (as, op, src_lw, src_sw, bcs_lw, bcs_sw, fluxb_lw, fluxb_sw, flux_lw, flux_sw)
        targs = typeof.(args)
        FT = eltype(as.p_lev)
        FTA1D = typeof(as.t_sfc)
        FTA2D = typeof(as.p_lev)
        @assert targs[1] <: AbstractAtmosphericState{FT, FTA1D}
        @assert targs[2] <: AbstractOpticalProps{FT, FTA2D}
        @assert targs[3] <: Union{AbstractSourceLW{FT, FTA1D, FTA2D}, Nothing}
        @assert targs[4] <: Union{SourceSW2Str{FT, FTA1D, FTA2D}, Nothing}
        @assert targs[5] <: Union{LwBCs{FT}, Nothing}
        @assert targs[6] <: Union{SwBCs{FT}, Nothing}
        @assert targs[7] <: Union{FluxLW{FT, FTA2D}, Nothing}
        @assert targs[8] <: Union{FluxSW{FT, FTA2D}, Nothing}
        @assert targs[9] <: Union{FluxLW{FT, FTA2D}, Nothing}
        @assert targs[10] <: Union{FluxSW{FT, FTA2D}, Nothing}
        return new{targs...}(args...)
    end
end

float_type(s::Solver) = eltype(s.as.p_lev)


Adapt.@adapt_structure Solver
end
