#####
##### Wide workspaces for a g-point-parallel solve.
#####

# Every group count a caller names goes through here. A count below one would
# pass every `ngpt_groups > 1` gate by failing it, leaving a serial solver where
# a grouped one was asked for, so it is an error rather than a clamp.
@inline function _check_ngpt_groups(ngpt_groups::Int)
    ngpt_groups >= 1 ||
        error("ngpt_groups must be at least 1, got $ngpt_groups")
    return nothing
end

# Planes of the grid's own float type one wide column holds per level, counted
# from what the constructors below build: three for the optics, three for the
# longwave source and two for the shortwave's, three flux planes in the longwave
# and four in the shortwave, and two partial accumulators of two planes each in
# the longwave and three each in the shortwave.
const LW_GPT_GROUP_PLANES = 3 + 3 + 3 + 2 * 2
const SW_GPT_GROUP_PLANES = 3 + 2 + 4 + 2 * 3
# Planes per band of the aerosol store: optical depth, single scattering albedo
# and asymmetry parameter.
const AEROSOL_BAND_PLANES = 3

"""
    gpt_group_bytes(FT, ncol, nlay, ngroups, naero_bands = 0)

Device bytes `LWGptGroups` and `SWGptGroups` together commit for `ngroups`
groups of `ncol` columns on `nlay` layers, which is the quantity
[`default_ngpt_groups`](@ref) bounds. A coupled run builds both bands, so this
prices both.

Each band adds a byte of cloud mask per level and, per wide column, a surface
source and a cloudy tally. The optics and the mask are over layers rather than
levels and are charged a level each, which keeps the count on the safe side of
the bound; `test/gpt_group_workspaces.jl` holds it to exactly that.

`naero_bands` is how many bands, over the longwave and the shortwave together,
the workspaces keep aerosol optics for. That store is three planes per band over
the layers of the `ncol` columns themselves: it does not widen with `ngroups`,
and it is counted exactly.
"""
function gpt_group_bytes(
    ::Type{FT},
    ncol,
    nlay,
    ngroups,
    naero_bands = 0,
) where {FT}
    planes = LW_GPT_GROUP_PLANES + SW_GPT_GROUP_PLANES
    per_level = planes * sizeof(FT) + 2 * sizeof(Bool)
    per_column = 2 * (sizeof(FT) + sizeof(Int32))
    aerosol = AEROSOL_BAND_PLANES * naero_bands * ncol * nlay * sizeof(FT)
    return ncol * ngroups * ((nlay + 1) * per_level + per_column) + aerosol
end

# An out-of-memory failure building these names a byte count and nothing about
# why the solve wanted it. Say what the group count costs and what to pass
# instead, then rethrow: a run that cannot afford the default wants an explicit
# `ngpt_groups`, not a quietly narrower solve.
function _try_gpt_group_alloc(
    f,
    grid_params::RRTMGPGridParams,
    ngroups::Int,
    naero_bands::Int = 0,
)
    try
        return f()
    catch
        (; ncol, nlay) = grid_params
        bytes = gpt_group_bytes(
            eltype(grid_params),
            ncol,
            nlay,
            ngroups,
            naero_bands,
        )
        hint = "pass `ngpt_groups` (1 solves the g-point loop serially)"
        @error "cannot allocate the g-point group workspaces" ngroups ncol nlay bytes hint
        rethrow()
    end
end

# Per-band aerosol optics for every column, or `nothing` when the solve keeps
# none. Aerosol optics do not depend on the g-point, so a grouped solve computes
# them once per band here and adds them to each of the band's g-points, instead
# of looking them up again for every one. Indexed `[gcol, glay, ibnd]`, so the
# threads of one group read it coalesced.
function _aerosol_band_store(grid_params::RRTMGPGridParams, naero_bands::Int)
    naero_bands > 0 || return nothing
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    plane() = DA{FT}(undef, ncol, nlay, naero_bands)
    # `AEROSOL_BAND_PLANES` of them, which is what `gpt_group_bytes` charges.
    return (; τ = plane(), ssa = plane(), g = plane())
end

# What the two bands' workspaces have in common: the widened grid parameters,
# the two partial accumulators, and the mask and cloud tally over wide columns.
# `planes` is how many flux planes a band accumulates, which is the only shape
# that differs between them.
function _gpt_group_workspaces(
    grid_params::RRTMGPGridParams,
    ngroups::Int,
    planes::Int,
)
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
    return (;
        wide_grid_params,
        acc_data = DA{FT}(undef, ncolw, nlay + 1, planes),
        acc_clear_data = DA{FT}(undef, ncolw, nlay + 1, planes),
        mask = DA{Bool}(undef, nlay, ncolw),
        cloudy = DA{Int32}(undef, ncolw),
    )
end

"""
    LWGptGroups(grid_params::RRTMGPGridParams, ngroups; params)

Workspaces for solving a column's longwave g-point loop in `ngroups`
independent pieces, covered by one grid rather than by one column's serial
loop.

At one thread per column the whole solve fits in a single wave of blocks, so
nothing is queued behind it: a block that finishes early frees its
multiprocessor to no one, and the launch lasts as long as its slowest column.
Columns are not equally expensive -- a column sweeps the cloudy increment only
where its mask is set, and in the shortwave a night column exits at once (see
[`SWGptGroups`](@ref)) -- so that difference is paid in idle warps instead of
being absorbed. Splitting the g-point loop multiplies the blocks, and the
scheduler refills a multiprocessor as each retires.

Every workspace the solve indexes by column is widened to `ncol * ngroups`, so a
thread's `(column, group)` pair selects its slice by one index and the kernel
body is the ungrouped one with that index substituted wherever it touches
scratch. Nothing is allocated or walked per group: the reduction recovers the
group from the wide column index.

`ngroups > 1` is implemented for CUDA. The CPU driver ignores these workspaces
and solves the whole g-point loop per column, which is correct either way.
"""
struct LWGptGroups{OP, SL, FB, MK, FX, CT, AE}
    ngroups::Int
    op::OP
    src::SL
    fluxb::FB
    mask::MK
    "all-sky partial accumulator over the wide columns, as views by flux name"
    acc::FX
    "clear-sky partial accumulator, as views by flux name"
    acc_clear::FX
    "`(ncol * ngroups)`: cloudy g-points seen by each group, summed for the cover"
    cloudy::CT
    "`(ncol, nlay, nbnd)` aerosol optics per band, or `nothing` to look them up per g-point"
    aero::AE
end
Adapt.@adapt_structure LWGptGroups

function LWGptGroups(
    grid_params::RRTMGPGridParams,
    ngroups::Int;
    params,
    naero_bands::Int = 0,
)
    _check_ngpt_groups(ngroups)
    _try_gpt_group_alloc(grid_params, ngroups, naero_bands) do
        # Up and down only: the reduction recomputes net from them, so a third
        # plane would be allocated per wide column and never written. The field
        # names match `FluxLW` because that is what these accumulate, but the
        # type does not, which is the point -- a `FluxLW` would carry the net
        # plane.
        (; wide_grid_params, acc_data, acc_clear_data, mask, cloudy) =
            _gpt_group_workspaces(grid_params, ngroups, 2)
        # `view(d, :, :, k)` is a contiguous block, so the accumulators index as
        # cheaply as the plain arrays the ungrouped path uses.
        fluxes(d) = (flux_up = view(d, :, :, 1), flux_dn = view(d, :, :, 2))
        LWGptGroups(
            ngroups,
            TwoStream(wide_grid_params),
            SourceLW2Str(wide_grid_params; params),
            FluxLW(wide_grid_params),
            mask,
            fluxes(acc_data),
            fluxes(acc_clear_data),
            cloudy,
            _aerosol_band_store(grid_params, naero_bands),
        )
    end
end

n_gpt_groups(::Nothing) = 1
n_gpt_groups(g::LWGptGroups) = g.ngroups

"""
    SWGptGroups(grid_params::RRTMGPGridParams, ngroups)

Workspaces for solving a column's shortwave g-point loop in `ngroups`
independent pieces, covered by one grid rather than by one column's serial
loop.

[`LWGptGroups`](@ref) carries the reasoning, which holds here unchanged.

Two things differ from the longwave. `FluxSW` carries a direct-beam profile, so
a group accumulates three planes rather than two. And the shortwave sweeps only
where the sun is up, so every group of a night column exits its loop at once,
which is the imbalance more blocks than one wave can absorb.
"""
struct SWGptGroups{OP, SS, FB, MK, FX, CT, AE}
    ngroups::Int
    op::OP
    src::SS
    fluxb::FB
    mask::MK
    "all-sky partial accumulator over the wide columns, as views by flux name"
    acc::FX
    "clear-sky partial accumulator, as views by flux name"
    acc_clear::FX
    "`(ncol * ngroups)`: cloudy g-points seen by each group, summed for the cover"
    cloudy::CT
    "`(ncol, nlay, nbnd)` aerosol optics per band, or `nothing` to look them up per g-point"
    aero::AE
end
Adapt.@adapt_structure SWGptGroups

function SWGptGroups(
    grid_params::RRTMGPGridParams,
    ngroups::Int;
    naero_bands::Int = 0,
)
    _check_ngpt_groups(ngroups)
    _try_gpt_group_alloc(grid_params, ngroups, naero_bands) do
        # Up, down and direct: the reduction recomputes net from up and down, so
        # a net plane would never be written (see `LWGptGroups`).
        (; wide_grid_params, acc_data, acc_clear_data, mask, cloudy) =
            _gpt_group_workspaces(grid_params, ngroups, 3)
        fluxes(d) = (
            flux_up = view(d, :, :, 1),
            flux_dn = view(d, :, :, 2),
            flux_dn_dir = view(d, :, :, 3),
        )
        SWGptGroups(
            ngroups,
            TwoStream(wide_grid_params),
            SourceSW2Str(wide_grid_params),
            FluxSW(wide_grid_params),
            mask,
            fluxes(acc_data),
            fluxes(acc_clear_data),
            cloudy,
            _aerosol_band_store(grid_params, naero_bands),
        )
    end
end

n_gpt_groups(g::SWGptGroups) = g.ngroups
