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
# the optics, three for the longwave source and two for the shortwave's, and
# three flux planes in the longwave and four in the shortwave.
const LW_GPT_GROUP_PLANES = 3 + 3 + 3
const SW_GPT_GROUP_PLANES = 3 + 2 + 4

"""
    gpt_group_bytes(FT, ncol, nlay, ngroups)

Device bytes that `LWGptGroups` and `SWGptGroups` together allocate for
`ngroups` groups of `ncol` columns on `nlay` layers, which
[`default_ngpt_groups`](@ref) bounds. Each band adds a byte of cloud mask per
level and a surface source per wide column, and the aerosol and cloud optics of
one band per column; the wide optics and the mask are over layers but are
charged a level each.
"""
function gpt_group_bytes(::Type{FT}, ncol, nlay, ngroups) where {FT}
    planes = LW_GPT_GROUP_PLANES + SW_GPT_GROUP_PLANES
    per_level = planes * sizeof(FT) + 2 * sizeof(Bool)
    per_column = 2 * sizeof(FT)
    band_optics = 2 * ncol * nlay * 6 * sizeof(FT)
    return ncol * ngroups * ((nlay + 1) * per_level + per_column) + band_optics
end

# The widened grid parameters, the mask over wide columns, and the aerosol and
# cloud optics (`τ`, `ssa` and `g` of each) of one band per column.
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
    mask = DA{Bool}(undef, nlay, ncolw)
    return (; wide_grid_params, mask, band = DA{FT}(undef, nlay, 6, ncol))
end

"""
    LWGptGroups(grid_params::RRTMGPGridParams, ngroups; params)

Workspaces for solving a column's longwave g-point loop with `ngroups` threads
(lanes of one warp), rather than by one thread's serial loop. `ngroups` is a
power of two up to 32, and a type parameter, so that the lanes' reductions
unroll.

Every workspace the solve indexes by column is widened to `ncol * ngroups`, so a
thread's `(column, lane)` pair selects its slice by one index, except `band`,
which holds the aerosol and cloud optics of one band per column, shared by the
column's lanes. `ngroups > 1` is implemented for CUDA; the CPU driver ignores
these workspaces.
"""
struct LWGptGroups{N, OP, SL, FB, MK, BA}
    op::OP
    src::SL
    fluxb::FB
    mask::MK
    band::BA
end
LWGptGroups{N}(
    op::OP,
    src::SL,
    fluxb::FB,
    mask::MK,
    band::BA,
) where {N, OP, SL, FB, MK, BA} =
    LWGptGroups{N, OP, SL, FB, MK, BA}(op, src, fluxb, mask, band)
Adapt.adapt_structure(to, g::LWGptGroups{N}) where {N} = LWGptGroups{N}(
    Adapt.adapt(to, g.op),
    Adapt.adapt(to, g.src),
    Adapt.adapt(to, g.fluxb),
    Adapt.adapt(to, g.mask),
    Adapt.adapt(to, g.band),
)

function LWGptGroups(grid_params::RRTMGPGridParams, ngroups::Int; params)
    (; wide_grid_params, mask, band) =
        _gpt_group_workspaces(grid_params, ngroups)
    return LWGptGroups{ngroups}(
        TwoStream(wide_grid_params),
        SourceLW2Str(wide_grid_params; params),
        FluxLW(wide_grid_params),
        mask,
        band,
    )
end

"""
    SWGptGroups(grid_params::RRTMGPGridParams, ngroups)

The shortwave counterpart of [`LWGptGroups`](@ref).
"""
struct SWGptGroups{N, OP, SS, FB, MK, BA}
    op::OP
    src::SS
    fluxb::FB
    mask::MK
    band::BA
end
SWGptGroups{N}(
    op::OP,
    src::SS,
    fluxb::FB,
    mask::MK,
    band::BA,
) where {N, OP, SS, FB, MK, BA} =
    SWGptGroups{N, OP, SS, FB, MK, BA}(op, src, fluxb, mask, band)
Adapt.adapt_structure(to, g::SWGptGroups{N}) where {N} = SWGptGroups{N}(
    Adapt.adapt(to, g.op),
    Adapt.adapt(to, g.src),
    Adapt.adapt(to, g.fluxb),
    Adapt.adapt(to, g.mask),
    Adapt.adapt(to, g.band),
)

function SWGptGroups(grid_params::RRTMGPGridParams, ngroups::Int)
    (; wide_grid_params, mask, band) =
        _gpt_group_workspaces(grid_params, ngroups)
    return SWGptGroups{ngroups}(
        TwoStream(wide_grid_params),
        SourceSW2Str(wide_grid_params),
        FluxSW(wide_grid_params),
        mask,
        band,
    )
end

n_gpt_groups(::Nothing) = 1
n_gpt_groups(::Union{LWGptGroups{N}, SWGptGroups{N}}) where {N} = N
