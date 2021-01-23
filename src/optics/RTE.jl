module RTE
using KernelAbstractions
using CUDA
using ..Device: array_type, array_device
using Adapt
using ..AtmosphericStates
using DocStringExtensions
import GaussQuadrature
using ..Sources
using ..Fluxes
using ..AngularDiscretizations
using ..Optics
using ..BCs

using CLIMAParameters
using CLIMAParameters.Planet: grav, R_d, cp_d

export gas_optics_gray_atmos!,
    update_profile_lw_kernel!, Solver, compute_gray_heating_rate!

struct Solver{
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    AS<:AbstractAtmosphericState{FT,FTA1D,FTA2D,I},
    OP<:AbstractOpticalProps{FT,FTA2D},
    SL<:Union{AbstractSource{FT,FTA2D},Nothing},
    SS<:Union{AbstractSource{FT,FTA2D},Nothing},
    BCL<:Union{AbstractBCs{FT,FTA1D},Nothing},
    BCS<:Union{AbstractBCs{FT,FTA1D},Nothing},
    AD<:Union{AngularDiscretization{FT,FTA1D,I},Nothing},
    FXBL<:Union{AbstractFlux{FT,FTA2D},Nothing},
    FXBS<:Union{AbstractFlux{FT,FTA2D},Nothing},
    FXL<:Union{AbstractFlux{FT,FTA2D},Nothing},
    FXS<:Union{AbstractFlux{FT,FTA2D},Nothing},
}
    as::AS         # atmoshperic state
    op::OP         # optical properties
    src_lw::SL    # source functions
    src_sw::SS    # source functions
    bcs_lw::BCL    # boundary conditions
    bcs_sw::BCS    # boundary conditions
    angle_disc::AD # angular discretization
    fluxb_lw::FXBL  # temporay storage for bandwise calculations
    fluxb_sw::FXBS  # temporay storage for bandwise calculations
    flux_lw::FXL   # fluxes for longwave problem
    flux_sw::FXS   # fluxes for longwave problem
end
Solver(
    as,
    op,
    src_lw,
    src_sw,
    bcs_lw,
    bcs_sw,
    angle_disc,
    fluxb_lw,
    fluxb_sw,
    flux_lw,
    flux_sw,
) = Solver{
    typeof(as),
    typeof(op),
    typeof(src_lw),
    typeof(src_sw),
    typeof(bcs_lw),
    typeof(bcs_sw),
    typeof(angle_disc),
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
    angle_disc,
    fluxb_lw,
    fluxb_sw,
    flux_lw,
    flux_sw,
)
Adapt.@adapt_structure Solver
end
