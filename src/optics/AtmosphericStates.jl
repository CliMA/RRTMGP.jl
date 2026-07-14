module AtmosphericStates
import ClimaComms
using Adapt
import ..Parameters as RP
import ..RRTMGPGridParams

using ..VolumeMixingRatios

export AbstractAtmosphericState,
    AtmosphericState,
    GrayAtmosphericState,
    CloudState,
    AerosolState,
    TransposedStateCache,
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

# Whole-column-set views of the packed layer data (each `(nlay, ncol)`); the
# per-column `(nlay,)` field slices used by the gas-optics kernels are the
# `col_dry_view`/`p_lay_view`/`t_lay_view`/`rel_hum_view` accessors below.
# view of column amounts of dry air [molecules per cm^2 of dry air]
@inline getview_col_dry(as::AtmosphericState) =
    @inbounds view(as.layerdata, 1, :, :)
# view of layer pressures [Pa, mb]
@inline getview_p_lay(as::AtmosphericState) =
    @inbounds view(as.layerdata, 2, :, :)
# view of layer temperatures [K]
@inline getview_t_lay(as::AtmosphericState) =
    @inbounds view(as.layerdata, 3, :, :)
# view of relative humidity
@inline getview_rel_hum(as::AtmosphericState) =
    @inbounds view(as.layerdata, 4, :, :)

"""
    TransposedStateCache{D2, D3}

Column-first copies of the hot [`AtmosphericState`](@ref) fields, refreshed
once per solve (by `refresh_transposed_state!`). The gas-optics kernels
re-read the layer pressures, temperatures, dry-air column amounts, and level
temperatures for every g-point; in the caller-owned `AtmosphericState` those
arrays are vertical-first, so neighboring GPU threads (one per column) read
memory `nlay` or `4 nlay` elements apart. Reading the transposed copies
instead makes those loads consecutive (coalesced) across the threads of a
warp.

The cache trades one `permutedims!` per solve plus `(4 nlay + nlev) ncol`
extra storage for coalesced reads in the g-point loop. Pass the same cache to
the longwave and shortwave workspaces to share the storage, or `nothing` to
opt out (the kernels then read the `AtmosphericState` directly, as before).
The workspace constructors default to a cache on the GPU and to `nothing` on
the CPU, where the state's vertical-first layout is already cache-friendly
for the per-column sweeps.

# Fields
- `layerdata`: `(ncol, nlay, 4)` copy of `AtmosphericState.layerdata`
  (dry-air column amount, layer pressure, layer temperature, relative
  humidity).
- `t_lev`: `(ncol, nlay + 1)` copy of `AtmosphericState.t_lev`.
"""
struct TransposedStateCache{D2, D3}
    layerdata::D3
    t_lev::D2
end
Adapt.@adapt_structure TransposedStateCache

function TransposedStateCache(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    return TransposedStateCache(
        DA{FT, 3}(undef, ncol, nlay, 4),
        DA{FT, 2}(undef, ncol, nlay + 1),
    )
end

# refresh_transposed_state!(cache, as): copy the current `AtmosphericState`
# layer data and level temperatures into the column-first `TransposedStateCache`
# `cache`. Called at the start of every spectral `solve_lw!`/`solve_sw!`; a
# no-op when `cache === nothing`.
@inline refresh_transposed_state!(::Nothing, as::AtmosphericState) = nothing
function refresh_transposed_state!(
    cache::TransposedStateCache,
    as::AtmosphericState,
)
    _transpose_into!(cache.layerdata, as.layerdata)
    _transpose_into!(cache.t_lev, as.t_lev)
    return nothing
end

# Device arrays: `permutedims!` maps to a single native transpose kernel.
_transpose_into!(dest::AbstractArray{T, 3}, src::AbstractArray{T, 3}) where {T} =
    (permutedims!(dest, src, (3, 2, 1)); nothing)
_transpose_into!(dest::AbstractArray{T, 2}, src::AbstractArray{T, 2}) where {T} =
    (permutedims!(dest, src, (2, 1)); nothing)
# Host arrays: explicit loops (`permutedims!` allocates a few bytes on
# Julia ≤ 1.11, tripping the zero-allocation tests).
function _transpose_into!(dest::Array{T, 3}, src::Array{T, 3}) where {T}
    nfld, nlay, ncol = size(src)
    @inbounds for f in 1:nfld, glay in 1:nlay, gcol in 1:ncol
        dest[gcol, glay, f] = src[f, glay, gcol]
    end
    return nothing
end
function _transpose_into!(dest::Array{T, 2}, src::Array{T, 2}) where {T}
    nlev, ncol = size(src)
    @inbounds for glev in 1:nlev, gcol in 1:ncol
        dest[gcol, glev] = src[glev, gcol]
    end
    return nothing
end

# Kernel-facing per-column accessors for the hot state fields (each a `(nlay,)`
# view of dry-air column amount / layer pressure / layer temperature / relative
# humidity, plus the `(nlev,)` level temperatures). With a cache, they read the
# coalesced column-first copies (`cache.layerdata[gcol, :, fld]`, stride-1 over
# columns); without one, they read the `AtmosphericState` arrays directly
# (uncoalesced — the opt-out for memory-constrained configurations). Both cases
# are plain `view`s, so they compile on the GPU (unlike an in-kernel
# `PermutedDimsArray`). The cache stores the fields in the `AtmosphericState`
# order: dry-air column amount 1, pressure 2, temperature 3, relative humidity 4.
@inline col_dry_view(cache::TransposedStateCache, as, gcol) =
    @inbounds view(cache.layerdata, gcol, :, 1)
@inline col_dry_view(::Nothing, as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, 1, :, gcol)
@inline p_lay_view(cache::TransposedStateCache, as, gcol) =
    @inbounds view(cache.layerdata, gcol, :, 2)
@inline p_lay_view(::Nothing, as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, 2, :, gcol)
@inline t_lay_view(cache::TransposedStateCache, as, gcol) =
    @inbounds view(cache.layerdata, gcol, :, 3)
@inline t_lay_view(::Nothing, as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, 3, :, gcol)
@inline rel_hum_view(cache::TransposedStateCache, as, gcol) =
    @inbounds view(cache.layerdata, gcol, :, 4)
@inline rel_hum_view(::Nothing, as::AtmosphericState, gcol) =
    @inbounds view(as.layerdata, 4, :, gcol)
@inline t_lev_view(cache::TransposedStateCache, as, gcol) =
    @inbounds view(cache.t_lev, gcol, :)
@inline t_lev_view(::Nothing, as::AtmosphericState, gcol) =
    @inbounds view(as.t_lev, :, gcol)

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
