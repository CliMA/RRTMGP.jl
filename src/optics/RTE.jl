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

using CLIMAParameters
using CLIMAParameters.Planet: grav, R_d, cp_d

export Solver

"""
    struct Solver{
        FT,
        I,
        FTA1D,
        FTA2D,
        AS,
        OP,
        SL,
        SS,
        BCL,
        BCS,
        FXBL,
        FXBS,
        FXL,
        FXS,
    }

The high-level RRTMGP data structure storing the atmospheric state, 
optical properties, sources, boundary conditions and fluxes configurations
for a given simulation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Solver{
    FT <: AbstractFloat,
    I <: Int,
    FTA1D <: AbstractArray{FT, 1},
    FTA2D <: AbstractArray{FT, 2},
    AS <: AbstractAtmosphericState{FT, I, FTA1D},
    OP <: AbstractOpticalProps{FT, FTA2D},
    SL <: Union{AbstractSourceLW{FT, FTA1D, FTA2D}, Nothing},
    SS <: Union{SourceSW2Str{FT, FTA1D, FTA2D}, Nothing},
    BCL <: Union{LwBCs{FT}, Nothing},
    BCS <: Union{SwBCs{FT}, Nothing},
    FXBL <: Union{FluxLW{FT, FTA2D}, Nothing},
    FXBS <: Union{FluxSW{FT, FTA2D}, Nothing},
    FXL <: Union{FluxLW{FT, FTA2D}, Nothing},
    FXS <: Union{FluxSW{FT, FTA2D}, Nothing},
}
    as::AS         # atmospheric state
    op::OP         # optical properties
    src_lw::SL     # source functions
    src_sw::SS     # source functions
    bcs_lw::BCL    # boundary conditions
    bcs_sw::BCS    # boundary conditions
    fluxb_lw::FXBL # temporary storage for bandwise calculations
    fluxb_sw::FXBS # temporary storage for bandwise calculations
    flux_lw::FXL   # fluxes for longwave problem
    flux_sw::FXS   # fluxes for shortwave problem
end

Solver(as, op, src_lw, src_sw, bcs_lw, bcs_sw, fluxb_lw, fluxb_sw, flux_lw, flux_sw) = Solver{
    eltype(as.p_lev),
    typeof(as.ncol),
    typeof(as.t_sfc),
    typeof(as.p_lev),
    typeof(as),
    typeof(op),
    typeof(src_lw),
    typeof(src_sw),
    typeof(bcs_lw),
    typeof(bcs_sw),
    typeof(fluxb_lw),
    typeof(fluxb_sw),
    typeof(flux_lw),
    typeof(flux_sw),
}(
    as,
    op,
    src_lw,
    src_sw,
    bcs_lw,
    bcs_sw,
    fluxb_lw,
    fluxb_sw,
    flux_lw,
    flux_sw,
)

Adapt.@adapt_structure Solver
end
