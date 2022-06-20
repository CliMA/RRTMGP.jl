module AtmosphericStates
using CUDA
using ..Device: array_type, array_device, CUDADevice, CPU
using DocStringExtensions
using Adapt
using GaussQuadrature
import ..Parameters as RP

using ..Vmrs

export AbstractAtmosphericState,
    AtmosphericState, ClearAtmosphericState, GrayAtmosphericState, setup_gray_as_pr_grid, setup_gray_as_alt_grid

abstract type AbstractAtmosphericState{FT, I, FTA1D} end

include("GrayAtmosphericStates.jl")

"""
    AtmosphericState{FT,FTA1D,FTA1DN,FTA2D,CLDP,CLDM,VMR,I} <:
        AbstractAtmosphericState{FT,I,FTA1D}

Atmospheric conditions, used to compute optical properties. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AtmosphericState{
    FT <: AbstractFloat,
    FTA1D <: AbstractArray{FT, 1},
    FTA1DN <: Union{AbstractArray{FT, 1}, Nothing},
    FTA2D <: AbstractArray{FT, 2},
    CLDP <: Union{AbstractArray{FT, 2}, Nothing},
    CLDM <: Union{AbstractArray{Bool, 2}, Nothing},
    VMR <: AbstractVmr{FT},
    I <: Int,
} <: AbstractAtmosphericState{FT, I, FTA1D}
    "longitude (`ncol`), optional"
    lon::FTA1DN
    "latitude (`ncol`), optional"
    lat::FTA1DN
    "Layer pressures `[Pa, mb]`; `(nlay,ncol)`"
    p_lay::FTA2D
    "Level pressures `[Pa, mb]`; `(nlay+1,ncol)`"
    p_lev::FTA2D
    "Layer temperatures `[K]`; `(nlay,ncol)`"
    t_lay::FTA2D
    "Level temperatures `[K]`; `(nlay+1,ncol)`"
    t_lev::FTA2D
    "Surface temperatures `[K]`; `(ncol)`"
    t_sfc::FTA1D
    "Number of molecules per cm^2 of dry air `(nlay, ncol)`"
    col_dry::FTA2D
    "volume mixing ratio of all relevant gases"
    vmr::VMR
    "effective radius of cloud liquid particles"
    cld_r_eff_liq::CLDP
    "effective radius of cloud ice particles"
    cld_r_eff_ice::CLDP
    "cloud water path"
    cld_path_liq::CLDP
    "cloud ice path"
    cld_path_ice::CLDP
    "cloud mask"
    cld_mask::CLDM
    "ice roughness, 1 = none, 2 = medium, 3 = rough"
    ice_rgh::I
    "Number of layers."
    nlay::I
    "Number of columns."
    ncol::I
    "Number of gases."
    ngas::I
end
Adapt.@adapt_structure AtmosphericState
end
