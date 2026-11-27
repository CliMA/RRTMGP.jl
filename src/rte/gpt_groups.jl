#####
##### Wide workspaces for a g-point-parallel solve.
#####

# The groups of a column are lanes of one warp, so their count is a power of two
# up to 32.
@inline function _check_ngpt_groups(ngpt_groups::Int)
    ispow2(ngpt_groups) && ngpt_groups <= 32 ||
        error("ngpt_groups must be a power of two up to 32, got $ngpt_groups")
    return nothing
end

# Planes of the grid's float type that a wide column holds per level: three for
# the optics, three for the longwave source and three for the shortwave's, and
# three flux planes in the longwave and four in the shortwave. The fused sweeps
# keep in them what each lane's sweep up leaves for its sweep down (see
# `LWGptGroups`).
const LW_GPT_GROUP_PLANES = 3 + 3 + 3
const SW_GPT_GROUP_PLANES = 3 + 3 + 4

# Planes of the band workspace `band` (levels x planes x columns), which the lanes
# of a column fill for each band and then share: the aerosol and cloud optics
# (`τ`, `ssa` and `g` of each) and the aerosol optical depths of its layers, and
# the parts of its gas optics that do not depend on the g-point: the packed
# interpolation indices, the interpolation fractions and column mixing ratios,
# the Planck function at the levels (longwave) or the water vapor mixing ratio
# (shortwave), the Planck function at the surface, and the
# scaling of each minor absorber.
const BAND_AEROSOL = 0
const BAND_CLOUD = 3
const BAND_GAS_INDICES = 9
const BAND_GAS_FRACTIONS = 10
const BAND_PLANCK = 16
const BAND_PLANCK_SFC = 17
const BAND_MINOR = 18
const BAND_MAX_MINOR = 11
# Three planes that hold, for every band, the humidity interpolation location
# and factor and the packed size bins of the aerosols of each layer (see
# `Optics.aerosol_searches`), which do not depend on the band.
const BAND_AEROSOL_SEARCHES = BAND_MINOR + BAND_MAX_MINOR
const BAND_PLANES = BAND_AEROSOL_SEARCHES + 2

"""
    check_band_layout(lookup)

Check that the gas optics lookup table `lookup` (a `LookUpLW` or `LookUpSW`) fits
the band workspace of the fused solves: that the interpolation indices of a
layer pack into the bits that the fused kernels give them (4 for the reference
temperature, 6 for the pressure, 4 for each `η`, 7 for the first minor
contributor of a band and 4 for their number), and that no band has more minor
contributors than `BAND_MAX_MINOR` in either part of the atmosphere.
"""
function check_band_layout(lookup)
    n_η, n_press = size(lookup.kmajor, 1), size(lookup.kmajor, 2)
    n_temp = length(lookup.ref_points.t_ref)
    bounds =
        (Array(lookup.minor_lower.bnd_st), Array(lookup.minor_upper.bnd_st))
    n_minor = maximum(b -> maximum(diff(b)), bounds)
    first_minor = maximum(b -> maximum(b), bounds)
    n_temp < 16 &&
    n_press < 64 &&
    n_η < 16 &&
    first_minor < 128 &&
    n_minor < 16 &&
    n_minor <= BAND_MAX_MINOR ||
        error("the gas optics lookup table does not fit the band workspace of \
               the fused solves: $n_temp temperatures, $n_press pressures, \
               $n_η η values, $n_minor minor contributors in a band")
    return nothing
end

"""
    gpt_group_bytes(FT, ncol, nlay, ngroups)

Device bytes that `LWGptGroups` and `SWGptGroups` together allocate for
`ngroups` groups of `ncol` columns on `nlay` layers, which
[`default_ngpt_groups`](@ref) bounds: the wide planes and surface sources, and
the band workspace of each column. The wide optics are over layers but are
charged a level each.
"""
function gpt_group_bytes(::Type{FT}, ncol, nlay, ngroups) where {FT}
    planes = LW_GPT_GROUP_PLANES + SW_GPT_GROUP_PLANES
    per_level = planes * sizeof(FT)
    per_column = 2 * sizeof(FT)
    band_optics = 2 * ncol * (nlay + 1) * BAND_PLANES * sizeof(FT)
    return ncol * ngroups * ((nlay + 1) * per_level + per_column) + band_optics
end

# The widened grid parameters, and the band workspace of each column (see
# `BAND_PLANES`).
function _gpt_group_workspaces(grid_params::RRTMGPGridParams, ngroups::Int)
    _check_ngpt_groups(ngroups)
    (; context, ncol, nlay, isothermal_boundary_layer) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    ncolw = ncol * ngroups
    # The inner constructor, because `nlay` here is already the total and the
    # keyword one would add the isothermal layer a second time.
    wide_grid_params = RRTMGPGridParams{FT, typeof(context)}(
        context,
        nlay,
        ncolw,
        isothermal_boundary_layer,
    )
    band = DA{FT}(undef, nlay + 1, BAND_PLANES, ncol)
    return (; wide_grid_params, band)
end

"""
    LWGptGroups(grid_params::RRTMGPGridParams, ngroups; params)

Workspaces for solving a column's longwave g-point loop with `ngroups` threads
(lanes of one warp), rather than by one thread's serial loop. `ngroups` is a
power of two up to 32, and a type parameter, so that the lanes' reductions
unroll.

Every workspace the solve indexes by column is widened to `ncol * ngroups`, so a
thread's `(column, lane)` pair selects its slice by one index, except `band`,
which holds what the g-points of one band share in each column: its aerosol and
cloud optics and parts of its gas optics (see `BAND_PLANES`). The fused sweep up
keeps in the wide planes what each lane's sweep down needs: the two coefficients
of each layer and the albedo and source of each level, of the clear sky in
`fluxb.flux_up`, `fluxb.flux_dn`, `src.albedo` and `src.src`, and of the cloudy
sky in `op.τ`, `op.ssa`, `fluxb.flux_net` and `src.lev_source`. The lanes'
cloud masks are bits in shared memory. `ngroups > 1` is implemented for CUDA; the
CPU driver ignores these workspaces.
"""
struct LWGptGroups{N, OP, SL, FB, BA}
    op::OP
    src::SL
    fluxb::FB
    band::BA
end
LWGptGroups{N}(op::OP, src::SL, fluxb::FB, band::BA) where {N, OP, SL, FB, BA} =
    LWGptGroups{N, OP, SL, FB, BA}(op, src, fluxb, band)
Adapt.adapt_structure(to, g::LWGptGroups{N}) where {N} = LWGptGroups{N}(
    Adapt.adapt(to, g.op),
    Adapt.adapt(to, g.src),
    Adapt.adapt(to, g.fluxb),
    Adapt.adapt(to, g.band),
)

function LWGptGroups(grid_params::RRTMGPGridParams, ngroups::Int; params)
    (; wide_grid_params, band) = _gpt_group_workspaces(grid_params, ngroups)
    return LWGptGroups{ngroups}(
        TwoStream(wide_grid_params),
        SourceLW2Str(wide_grid_params; params),
        FluxLW(wide_grid_params),
        band,
    )
end

"""
    SWGptGroups(grid_params::RRTMGPGridParams, ngroups)

The shortwave counterpart of [`LWGptGroups`](@ref). Its fused sweep down keeps
the clear sky's optics in `op` and both skies' direct beams in
`fluxb.flux_dn_dir` and `fluxb.flux_net`; its sweep up then keeps the two
coefficients of each layer of the clear sky in `op.τ` and `op.ssa` and of the
cloudy sky in `op.g` and the third plane of `src.leveldata`, and the albedo and
source of each level of the clear sky in `src.albedo` and `src.src` and of the
cloudy sky in `fluxb.flux_up` and `fluxb.flux_dn`.
"""
struct SWGptGroups{N, OP, SS, FB, BA}
    op::OP
    src::SS
    fluxb::FB
    band::BA
end
SWGptGroups{N}(op::OP, src::SS, fluxb::FB, band::BA) where {N, OP, SS, FB, BA} =
    SWGptGroups{N, OP, SS, FB, BA}(op, src, fluxb, band)
Adapt.adapt_structure(to, g::SWGptGroups{N}) where {N} = SWGptGroups{N}(
    Adapt.adapt(to, g.op),
    Adapt.adapt(to, g.src),
    Adapt.adapt(to, g.fluxb),
    Adapt.adapt(to, g.band),
)

function SWGptGroups(grid_params::RRTMGPGridParams, ngroups::Int)
    (; wide_grid_params, band) = _gpt_group_workspaces(grid_params, ngroups)
    return SWGptGroups{ngroups}(
        TwoStream(wide_grid_params),
        sw_group_sources(wide_grid_params),
        FluxSW(wide_grid_params),
        band,
    )
end

n_gpt_groups(::Nothing) = 1
n_gpt_groups(::Union{LWGptGroups{N}, SWGptGroups{N}}) where {N} = N

# The shortwave sources over the wide columns, with a third level plane, which the
# fused solve uses as scratch.
function sw_group_sources(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    sfc_source = DA{FT, 1}(undef, ncol)
    leveldata = Fluxes._coalesced_3d(DA, FT, ncol, nlay + 1, 3)
    albedo = view(leveldata, :, :, 1)
    src = view(leveldata, :, :, 2)
    return SourceSW2Str{typeof(sfc_source), typeof(leveldata), typeof(albedo)}(
        sfc_source,
        leveldata,
        albedo,
        src,
    )
end
