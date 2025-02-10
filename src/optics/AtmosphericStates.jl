module AtmosphericStates
import ClimaComms
using DocStringExtensions
using Adapt
import ..Parameters as RP

using ..Vmrs

export AbstractAtmosphericState,
    AtmosphericState,
    GrayAtmosphericState,
    CloudState,
    AerosolState,
    setup_gray_as_pr_grid,
    MaxRandomOverlap,
    AbstractCloudMask,
    GrayOpticalThicknessSchneider2004,
    GrayOpticalThicknessOGorman2008,
    AbstractGrayOpticalThickness

abstract type AbstractAtmosphericState end

include("gray_atmospheric_states.jl")

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
struct AtmosphericState{FTA1D, FTA1DN, FTA2D, D, VMR, CLD, AER} <: AbstractAtmosphericState
    "longitude, in degrees (`ncol`), optional"
    lon::FTA1DN
    "latitude, in degrees (`ncol`), optional"
    lat::FTA1DN
    "storage for col_dry [`molecules per cm^2 of dry air`], play `[Pa, mb]`, tlay `[K]`, rel_hum; `(4, nlay, ncol)`"
    layerdata::D
    "Level pressures `[Pa, mb]`; `(nlay+1,ncol)`"
    p_lev::FTA2D
    "Level temperatures `[K]`; `(nlay+1,ncol)`"
    t_lev::FTA2D
    "Surface temperatures `[K]`; `(ncol)`"
    t_sfc::FTA1D
    "volume mixing ratio of all relevant gases"
    vmr::VMR
    "cloud state"
    cloud_state::CLD
    "aerosol state"
    aerosol_state::AER
end
Adapt.@adapt_structure AtmosphericState

# Number of layers
@inline get_nlay(as::AtmosphericState) = size(as.layerdata, 2)
# Number of columns
@inline get_ncol(as::AtmosphericState) = size(as.layerdata, 3)
# Number of layers and columns
@inline get_dims(as::AtmosphericState) = size(as.layerdata, 2), size(as.layerdata, 3)

@inline getview_layerdata(as::AtmosphericState, gcol) = @inbounds view(as.layerdata, :, :, gcol)

# view of column amounts of dry air [molecules per cm^2 of dry air]
@inline getview_col_dry(as::AtmosphericState) = @inbounds view(as.layerdata, 1, :, :)
@inline getview_col_dry(as::AtmosphericState, gcol) = @inbounds view(as.layerdata, 1, :, gcol)

# view of layer pressures [Pa, mb]
@inline getview_p_lay(as::AtmosphericState) = @inbounds view(as.layerdata, 2, :, :)
@inline getview_p_lay(as::AtmosphericState, gcol) = @inbounds view(as.layerdata, 2, :, gcol)

# view of layer temperatures [K]
@inline getview_t_lay(as::AtmosphericState) = @inbounds view(as.layerdata, 3, :, :)
@inline getview_t_lay(as::AtmosphericState, gcol) = @inbounds view(as.layerdata, 3, :, gcol)

# view of relative humidity
@inline getview_rel_hum(as::AtmosphericState) = @inbounds view(as.layerdata, 4, :, :)
@inline getview_rel_hum(as::AtmosphericState, gcol) = @inbounds view(as.layerdata, 4, :, gcol)

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

"""
    AerosolState{A, B, D}

Aerosol state, used to compute optical properties.
"""
struct AerosolState{A, B, D}
    "shortwave aerosol optical depth"
    aod_sw_ext::A
    "shortwave aerosol optical depth (scattering component)"
    aod_sw_sca::A
    "aerosol mask, = true if any aerosol is present"
    aero_mask::B
    "aerosol size (microns)"
    aero_size::D
    "aerosol mass column (kg/m2)"
    aero_mass::D
end
Adapt.@adapt_structure AerosolState

end
