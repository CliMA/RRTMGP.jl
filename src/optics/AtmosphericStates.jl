module AtmosphericStates
import ClimaComms
using Adapt
import ..Parameters as RP

using ..VolumeMixingRatios

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

"""
    AbstractAtmosphericState

Abstract type for atmospheric states: [`AtmosphericState`](@ref) for the
spectral (correlated-``k``) methods and [`GrayAtmosphericState`](@ref) for
gray radiation.
"""
abstract type AbstractAtmosphericState end

include("gray_atmospheric_states.jl")

"""
    AbstractCloudMask

Abstract type for cloud-overlap sampling masks; see [`MaxRandomOverlap`](@ref).
"""
abstract type AbstractCloudMask end

"""
    MaxRandomOverlap

Maximum-random cloud overlap for McICA sampling: clouds in adjacent layers
overlap maximally, and clouds separated by clear sky overlap randomly. Used
by `build_cloud_mask!` to sample the g-point cloud masks from the cloud
fraction.
"""
struct MaxRandomOverlap <: AbstractCloudMask end
Adapt.@adapt_structure MaxRandomOverlap

"""
    AtmosphericState{FTA1D, FTA1DN, FTA2D, D, VMR, CLD, AER} <:
        AbstractAtmosphericState

Atmospheric conditions, used to compute optical properties.

# Fields
- `lon`: Longitude in degrees `(ncol)`; optional.
- `lat`: Latitude in degrees `(ncol)`; optional.
- `layerdata`: Storage for `col_dry` (column amount of dry air [molecules/cm²]),
  layer pressures [Pa, mb], layer temperatures [K], and relative humidity;
  `(4, nlay, ncol)`.
- `p_lev`: Level pressures [Pa, mb] `(nlay+1, ncol)`.
- `t_lev`: Level temperatures [K] `(nlay+1, ncol)`.
- `t_sfc`: Surface temperatures [K] `(ncol)`.
- `vmr`: Volume mixing ratios of all relevant gases.
- `cloud_state`: Cloud state.
- `aerosol_state`: Aerosol state.
"""
struct AtmosphericState{FTA1D, FTA1DN, FTA2D, D, VMR, CLD, AER} <:
       AbstractAtmosphericState
    lon::FTA1DN
    lat::FTA1DN
    layerdata::D
    p_lev::FTA2D
    t_lev::FTA2D
    t_sfc::FTA1D
    vmr::VMR
    cloud_state::CLD
    aerosol_state::AER
end
Adapt.@adapt_structure AtmosphericState

# Number of layers
@inline get_nlay(as::AtmosphericState) = size(as.layerdata, 2)
# Number of columns
@inline get_ncol(as::AtmosphericState) = size(as.layerdata, 3)
# Number of layers and columns
@inline get_dims(as::AtmosphericState) =
    size(as.layerdata, 2), size(as.layerdata, 3)

@inline getview_layerdata(as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, :, :, gcol)

# view of column amounts of dry air [molecules per cm^2 of dry air]
@inline getview_col_dry(as::AtmosphericState) =
    @inbounds view(as.layerdata, 1, :, :)
@inline getview_col_dry(as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, 1, :, gcol)

# view of layer pressures [Pa, mb]
@inline getview_p_lay(as::AtmosphericState) =
    @inbounds view(as.layerdata, 2, :, :)
@inline getview_p_lay(as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, 2, :, gcol)

# view of layer temperatures [K]
@inline getview_t_lay(as::AtmosphericState) =
    @inbounds view(as.layerdata, 3, :, :)
@inline getview_t_lay(as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, 3, :, gcol)

# view of relative humidity
@inline getview_rel_hum(as::AtmosphericState) =
    @inbounds view(as.layerdata, 4, :, :)
@inline getview_rel_hum(as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, 4, :, gcol)

"""
    CloudState{CD, CF, CC, CM, CMT}

Cloud state, used to compute optical properties.

# Fields
- `cld_r_eff_liq`: Effective radius of cloud liquid particles.
- `cld_r_eff_ice`: Effective radius of cloud ice particles.
- `cld_path_liq`: Cloud water path.
- `cld_path_ice`: Cloud ice path.
- `cld_frac`: Cloud fraction.
- `cld_cover_sw`: McICA effective shortwave cloud cover in `[0, 1]`; `(ncol,)`, or
  `nothing` if unused.
- `cld_cover_lw`: McICA effective longwave cloud cover in `[0, 1]`; `(ncol,)`, or
  `nothing` if unused.
- `mask_lw`: Cloud mask (longwave); `true` if clouds are present.
- `mask_sw`: Cloud mask (shortwave); `true` if clouds are present.
- `mask_type`: Cloud mask type.
- `ice_rgh`: Ice roughness; 1 = none, 2 = medium, 3 = rough.
"""
struct CloudState{CD, CF, CC, CM, CMT}
    cld_r_eff_liq::CD
    cld_r_eff_ice::CD
    cld_path_liq::CD
    cld_path_ice::CD
    cld_frac::CF
    cld_cover_sw::CC
    cld_cover_lw::CC
    mask_lw::CM
    mask_sw::CM
    mask_type::CMT
    ice_rgh::Int
end
Adapt.@adapt_structure CloudState

# Convenience constructor: callers that don't allocate cld_cover arrays
# can use this 9-argument form where both covers default to `nothing`.
function CloudState(
    cld_r_eff_liq,
    cld_r_eff_ice,
    cld_path_liq,
    cld_path_ice,
    cld_frac,
    mask_lw,
    mask_sw,
    mask_type,
    ice_rgh,
)
    return CloudState(
        cld_r_eff_liq,
        cld_r_eff_ice,
        cld_path_liq,
        cld_path_ice,
        cld_frac,
        nothing,
        nothing,
        mask_lw,
        mask_sw,
        mask_type,
        Int(ice_rgh),
    )
end


"""
    AerosolState{A, B, D}

Aerosol state, used to compute optical properties.

# Fields
- `aod_sw_ext`: Shortwave aerosol optical depth.
- `aod_sw_sca`: Shortwave aerosol optical depth (scattering component).
- `aero_mask`: Aerosol mask; `true` if any aerosol is present.
- `aero_size`: Aerosol size [microns].
- `aero_mass`: Aerosol column mass [kg/m²].
"""
struct AerosolState{A, B, D}
    aod_sw_ext::A
    aod_sw_sca::A
    aero_mask::B
    aero_size::D
    aero_mass::D
end
Adapt.@adapt_structure AerosolState

end
