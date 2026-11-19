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

The gas concentrations, cloud inputs and aerosol inputs are re-read for every
g-point too. `TransposedStateCache(grid_params, as)` also allocates
column-first copies of whichever of those `as` carries, and
[`kernel_state`](@ref) hands the kernels a state that reads them. Each copy is
wrapped in a `PermutedDimsArray`, so it is indexed exactly like the array it
copies and the kernels need no second code path.
`TransposedStateCache(grid_params)` leaves the three `nothing`, and the kernels
then read those inputs from the `AtmosphericState` as before. A cache built
from a state has to be used with a state of the same shape: the same grid, and
the same optional parts (clouds, snow, aerosols).

# Fields
- `layerdata`: `(ncol, nlay, 4)` copy of `AtmosphericState.layerdata`
  (dry-air column amount, layer pressure, layer temperature, relative
  humidity).
- `t_lev`: `(ncol, nlay + 1)` copy of `AtmosphericState.t_lev`.
- `vmr`: column-first copies of the spatially varying volume mixing ratios, as
  a `NamedTuple` named like the fields they copy, or `nothing`.
- `cloud`: column-first copies of the cloud effective radii, water paths and
  cloud fraction, or `nothing`. The cloud masks are not copied: the kernels
  write them.
- `aerosol`: column-first copies of `aero_size` and `aero_mass`, or `nothing`.
"""
struct TransposedStateCache{D2, D3, V, C, A}
    layerdata::D3
    t_lev::D2
    vmr::V
    cloud::C
    aerosol::A
end
Adapt.@adapt_structure TransposedStateCache
TransposedStateCache(layerdata, t_lev) =
    TransposedStateCache(layerdata, t_lev, nothing, nothing, nothing)

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
    _refresh_inputs!(cache.vmr, as.vmr)
    _refresh_inputs!(cache.cloud, as.cloud_state)
    _refresh_inputs!(cache.aerosol, as.aerosol_state)
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
    CloudState{CD, PD, CF, CC, CM, CMT, SD}

Cloud state, used to compute optical properties.

# Fields
- `cld_r_eff_liq`: Effective radius of cloud liquid particles.
- `cld_r_eff_ice`: Effective radius of cloud ice particles.
- `cld_path_liq`: Cloud water path.
- `cld_path_ice`: Cloud ice path.
- `cld_r_eff_sno`: Effective radius of snow particles, or `nothing`.
- `cld_path_sno`: Snow water path, or `nothing` to disable snow.
- `cld_frac`: Cloud fraction.
- `cld_cover_sw`: McICA effective shortwave cloud cover in `[0, 1]`; `(ncol,)`, or
  `nothing` if unused.
- `cld_cover_lw`: McICA effective longwave cloud cover in `[0, 1]`; `(ncol,)`, or
  `nothing` if unused.
- `mask_lw`: Cloud mask (longwave); `true` if clouds are present. Also gates the
  snow increment — snow radiates only in sampled-cloudy layers.
- `mask_sw`: Cloud mask (shortwave); `true` if clouds are present. Also gates the
  snow increment.
- `mask_type`: Cloud mask type.
- `ice_rgh`: Ice roughness; 1 = none, 2 = medium, 3 = rough.
- `mcica_key`: One-element array holding the McICA key; the cloud mask is a
  function of `(mcica_key, column, g-point, layer)`. Set by
  `update_fluxes!(s, seedval)`, and must change between radiation steps. A key
  rather than a seed: it indexes a stateless draw, so the result does not depend
  on how many draws came before.
"""
struct CloudState{CD, PD, CF, CC, CM, CMT, SD}
    cld_r_eff_liq::CD
    cld_r_eff_ice::CD
    cld_path_liq::CD
    cld_path_ice::CD
    cld_r_eff_sno::PD
    cld_path_sno::PD
    cld_frac::CF
    cld_cover_sw::CC
    cld_cover_lw::CC
    mask_lw::CM
    mask_sw::CM
    mask_type::CMT
    ice_rgh::Int
    mcica_key::SD
end
Adapt.@adapt_structure CloudState

# McICA key salts, XORed into the key by `update_fluxes!` (rank) and
# `_build_cloud_mask!` (band). They must not cancel. `MCICA_RANK_SALT` is odd
# and so invertible mod 2^32, which means ANY band salt equals some rank's salt
# at exactly one rank index; the constant is chosen so that index is
# unreachable. Using `0x9e3779b9` for both put it at rank 2, where rank 1's
# shortwave key equaled rank 2's longwave key.
const MCICA_RANK_SALT = 0x9e3779b9
const MCICA_BAND_SALT_LW = 0x00000000
const MCICA_BAND_SALT_SW = 0x85ebca6b

# Without a key, the buffer is derived from `cld_frac` so it lands on the same
# device.
function CloudState(
    cld_r_eff_liq,
    cld_r_eff_ice,
    cld_path_liq,
    cld_path_ice,
    cld_frac,
    cld_cover_sw,
    cld_cover_lw,
    mask_lw,
    mask_sw,
    mask_type,
    ice_rgh::Integer;
    cld_r_eff_sno = nothing,
    cld_path_sno = nothing,
)
    # Snow needs both its path and radius, or neither: the optics kernel indexes
    # cld_r_eff_sno per layer, so a missing one would fault mid-kernel on GPU.
    isnothing(cld_path_sno) == isnothing(cld_r_eff_sno) || error(
        "snow optics need both cld_path_sno and cld_r_eff_sno, or neither",
    )
    mcica_key = similar(cld_frac, UInt32, 1)
    fill!(mcica_key, UInt32(0))
    return CloudState(
        cld_r_eff_liq,
        cld_r_eff_ice,
        cld_path_liq,
        cld_path_ice,
        cld_r_eff_sno,
        cld_path_sno,
        cld_frac,
        cld_cover_sw,
        cld_cover_lw,
        mask_lw,
        mask_sw,
        mask_type,
        Int(ice_rgh),
        mcica_key,
    )
end

# Convenience constructor: callers that don't allocate cld_cover arrays
# can use this form where both covers default to `nothing`.
function CloudState(
    cld_r_eff_liq,
    cld_r_eff_ice,
    cld_path_liq,
    cld_path_ice,
    cld_frac,
    mask_lw,
    mask_sw,
    mask_type,
    ice_rgh;
    cld_r_eff_sno = nothing,
    cld_path_sno = nothing,
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
        Int(ice_rgh);
        cld_r_eff_sno,
        cld_path_sno,
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

#####
##### Column-first copies of the per-g-point inputs
#####

# A column-first buffer shaped for `x`, wrapped so that it is indexed like `x`:
# `(nlay, ncol)` is stored `(ncol, nlay)`, and `(n, nlay, ncol)` is stored
# `(ncol, nlay, n)`, the layout `_transpose_into!` fills.
_column_first(::Nothing) = nothing
_column_first(x::AbstractArray{T, 2}) where {T} =
    PermutedDimsArray(similar(x, size(x, 2), size(x, 1)), (2, 1))
_column_first(x::AbstractArray{T, 3}) where {T} = PermutedDimsArray(
    similar(x, size(x, 3), size(x, 2), size(x, 1)),
    (3, 2, 1),
)

_refresh_field!(::Nothing, ::Nothing) = nothing
_refresh_field!(dest::PermutedDimsArray, src::AbstractArray) =
    _transpose_into!(parent(dest), src)
# A copy without its source, or a source without its copy: the cache was built
# from a state that carried different optional inputs (snow, say).
_refresh_field!(dest, src) = error(
    "this TransposedStateCache was built from an AtmosphericState with \
     different optional inputs than the one being solved; rebuild it with \
     `TransposedStateCache(grid_params, as)`",
)

# What is copied for each part of the state. Only inputs the kernels read are
# here; everything they write (masks, cloud cover, optical depths) stays in the
# caller's arrays, which `kernel_state` passes through.
_column_first_inputs(::Nothing) = nothing
_column_first_inputs(vmr::VmrGM) = (;
    vmr_h2o = _column_first(vmr.vmr_h2o),
    vmr_o3 = _column_first(vmr.vmr_o3),
)
_column_first_inputs(vmr::Vmr) = (; vmr = _column_first(vmr.vmr))
_column_first_inputs(cs::CloudState) = (;
    cld_r_eff_liq = _column_first(cs.cld_r_eff_liq),
    cld_r_eff_ice = _column_first(cs.cld_r_eff_ice),
    cld_path_liq = _column_first(cs.cld_path_liq),
    cld_path_ice = _column_first(cs.cld_path_ice),
    cld_r_eff_sno = _column_first(cs.cld_r_eff_sno),
    cld_path_sno = _column_first(cs.cld_path_sno),
    cld_frac = _column_first(cs.cld_frac),
)
_column_first_inputs(aero::AerosolState) = (;
    aero_size = _column_first(aero.aero_size),
    aero_mass = _column_first(aero.aero_mass),
)

function TransposedStateCache(
    grid_params::RRTMGPGridParams,
    as::AtmosphericState,
)
    (; layerdata, t_lev) = TransposedStateCache(grid_params)
    return TransposedStateCache(
        layerdata,
        t_lev,
        _column_first_inputs(as.vmr),
        _column_first_inputs(as.cloud_state),
        _column_first_inputs(as.aerosol_state),
    )
end

# A cache without copies of some part, or a state without that part, has
# nothing to refresh.
_refresh_inputs!(::Nothing, ::Nothing) = nothing
_refresh_inputs!(::Nothing, _) = nothing
_refresh_inputs!(_, ::Nothing) = nothing
function _refresh_inputs!(copies::NamedTuple, vmr::VmrGM)
    _refresh_field!(copies.vmr_h2o, vmr.vmr_h2o)
    _refresh_field!(copies.vmr_o3, vmr.vmr_o3)
    return nothing
end
_refresh_inputs!(copies::NamedTuple, vmr::Vmr) =
    _refresh_field!(copies.vmr, vmr.vmr)
function _refresh_inputs!(copies::NamedTuple, cs::CloudState)
    _refresh_field!(copies.cld_r_eff_liq, cs.cld_r_eff_liq)
    _refresh_field!(copies.cld_r_eff_ice, cs.cld_r_eff_ice)
    _refresh_field!(copies.cld_path_liq, cs.cld_path_liq)
    _refresh_field!(copies.cld_path_ice, cs.cld_path_ice)
    _refresh_field!(copies.cld_r_eff_sno, cs.cld_r_eff_sno)
    _refresh_field!(copies.cld_path_sno, cs.cld_path_sno)
    _refresh_field!(copies.cld_frac, cs.cld_frac)
    return nothing
end
function _refresh_inputs!(copies::NamedTuple, aero::AerosolState)
    _refresh_field!(copies.aero_size, aero.aero_size)
    _refresh_field!(copies.aero_mass, aero.aero_mass)
    return nothing
end

# The part of the state the kernels should read: the caller's own where the
# cache holds no copies, otherwise the same struct with its inputs replaced by
# the copies. The well-mixed concentrations, the masks, the cloud cover, the
# optical depths and the McICA key are the caller's arrays in both cases, so
# what the kernels write lands where the caller reads it.
@inline _kernel_inputs(::Nothing, part) = part
@inline _kernel_inputs(::NamedTuple, ::Nothing) = nothing
@inline _kernel_inputs(copies::NamedTuple, vmr::VmrGM) =
    VmrGM(copies.vmr_h2o, copies.vmr_o3, vmr.vmr)
@inline _kernel_inputs(copies::NamedTuple, ::Vmr) = Vmr(copies.vmr)
@inline _kernel_inputs(copies::NamedTuple, cs::CloudState) = CloudState(
    copies.cld_r_eff_liq,
    copies.cld_r_eff_ice,
    copies.cld_path_liq,
    copies.cld_path_ice,
    copies.cld_r_eff_sno,
    copies.cld_path_sno,
    copies.cld_frac,
    cs.cld_cover_sw,
    cs.cld_cover_lw,
    cs.mask_lw,
    cs.mask_sw,
    cs.mask_type,
    cs.ice_rgh,
    cs.mcica_key,
)
@inline _kernel_inputs(copies::NamedTuple, aero::AerosolState) = AerosolState(
    aero.aod_sw_ext,
    aero.aod_sw_sca,
    aero.aero_mask,
    copies.aero_size,
    copies.aero_mass,
)

"""
    kernel_state(cache, as::AtmosphericState)

The state the solve kernels read: `as` itself when `cache` is `nothing`,
otherwise `as` with the inputs `cache` holds column-first copies of replaced by
those copies. Call it after `refresh_transposed_state!`, which fills them.

Every array the kernels write is shared with `as`, so the result is only a
different route to the same inputs and never needs copying back.
"""
@inline kernel_state(::Nothing, as) = as
@inline kernel_state(cache::TransposedStateCache, as::AtmosphericState) =
    AtmosphericState(
        as.lon,
        as.lat,
        as.layerdata,
        as.p_lev,
        as.t_lev,
        as.t_sfc,
        _kernel_inputs(cache.vmr, as.vmr),
        _kernel_inputs(cache.cloud, as.cloud_state),
        _kernel_inputs(cache.aerosol, as.aerosol_state),
    )

end
