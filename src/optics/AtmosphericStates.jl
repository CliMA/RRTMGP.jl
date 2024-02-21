module AtmosphericStates
using CUDA
import ClimaComms
using DocStringExtensions
using Adapt
using GaussQuadrature
import ..Parameters as RP

using ..Vmrs

export AbstractAtmosphericState,
    AtmosphericState,
    GrayAtmosphericState,
    CloudState,
    setup_gray_as_pr_grid,
    setup_gray_as_alt_grid,
    MaxRandomOverlap,
    AbstractCloudMask,
    GrayOpticalThicknessSchneider2004,
    GrayOpticalThicknessOGorman2008,
    AbstractGrayOpticalThickness

abstract type AbstractAtmosphericState end

include("GrayAtmosphericStates.jl")

abstract type AbstractCloudMask end

struct MaxRandomOverlap <: AbstractCloudMask end
Adapt.@adapt_structure MaxRandomOverlap

"""
    AtmosphericState{FTA1D,FTA1DN,FTA2D,CLDP,CLDM,VMR} <:
        AbstractAtmosphericState

Atmospheric conditions, used to compute optical properties. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AtmosphericState{FTA1D, FTA1DN, FTA2D, VMR, CLD} <: AbstractAtmosphericState
    "longitude, in degrees (`ncol`), optional"
    lon::FTA1DN
    "latitude, in degrees (`ncol`), optional"
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
    "cloud state"
    cloud_state::CLD
end
Adapt.@adapt_structure AtmosphericState

# Number of layers
@inline get_nlay(as::AtmosphericState) = size(as.p_lay, 1)
# Number of columns
@inline get_ncol(as::AtmosphericState) = size(as.p_lay, 2)
# Number of layers and columns
@inline get_dims(as::AtmosphericState) = size(as.p_lay)

"""
    CloudState{CD, CF, CM, CMT}

Cloud state, used to compute optical properties.
"""
struct CloudState{CD, CF, CM, CMT}
    "effective radius of cloud liquid particles"
    cld_r_eff_liq::CD
    "effective radius of cloud ice particles"
    cld_r_eff_ice::CD
    "cloud water path"
    cld_path_liq::CD
    "cloud ice path"
    cld_path_ice::CD
    "cloud fraction"
    cld_frac::CF
    "cloud mask (longwave), = true if clouds are present"
    mask_lw::CM
    "cloud mask (shortwave), = true if clouds are present"
    mask_sw::CM
    "cloud mask type"
    mask_type::CMT
    "ice roughness, 1 = none, 2 = medium, 3 = rough"
    ice_rgh::Int
end
Adapt.@adapt_structure CloudState
end
