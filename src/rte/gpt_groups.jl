#####
##### Per-group workspaces for a g-point-parallel solve.
#####

"""
    LWGptGroups(grid_params::RRTMGPGridParams, ngroups; params)

Workspaces for solving a column's longwave g-point loop in `ngroups`
independent pieces, one kernel each.

At one thread per column the whole grid is a single wave -- 96 blocks of 256
threads against 108 SMs, 8 resident warps apiece -- so there is no second wave
to hide the latency behind and removing work from a thread buys almost nothing.
Splitting the g-point loop multiplies the resident parallelism instead, which
is the restructure priced in section 4l of the profiling study.

Each group needs its own `op`, `src`, `fluxb` and cloud mask, because those are
per-thread workspaces that the drivers index by column, and its own flux
accumulators for a reduction to sum afterwards. The scratch is allocated per
group rather than sliced out of one wide array so the solve kernel sees exactly
the types the ungrouped one does and pays no extra registers for the split; the
accumulators do carry a trailing group dimension, because the reduction has to
walk it.

`ngroups > 1` is implemented for CUDA. The CPU driver ignores these workspaces
and solves the whole g-point loop per column, which is correct either way.
"""
struct LWGptGroups{N, OP, SL, FB, FX, A4, MK, CT}
    op::NTuple{N, OP}
    src::NTuple{N, SL}
    fluxb::NTuple{N, FB}
    mask::NTuple{N, MK}
    "all-sky partial accumulators, views into `acc_data`"
    acc::NTuple{N, FX}
    "clear-sky partial accumulators, views into `acc_clear_data`"
    acc_clear::NTuple{N, FX}
    "`(ncol, nlev, 3, ngroups)`; the field axis is up, down, net as in `FluxLW`"
    acc_data::A4
    acc_clear_data::A4
    "`(ncol, ngroups)`: cloudy g-points seen by each group, summed for the cover"
    cloudy::CT
end
Adapt.@adapt_structure LWGptGroups

function LWGptGroups(grid_params::RRTMGPGridParams, ngroups::Int; params)
    ngroups ≥ 1 || error("ngroups must be at least 1, got $ngroups")
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    nlev = nlay + 1
    acc_data = DA{FT}(undef, ncol, nlev, 3, ngroups)
    acc_clear_data = DA{FT}(undef, ncol, nlev, 3, ngroups)
    # `view(d, :, :, k, g)` is a contiguous block, so the accumulators index as
    # cheaply as the plain arrays the ungrouped path uses.
    fluxes(d) = ntuple(ngroups) do g
        FluxLW(view(d, :, :, 1, g), view(d, :, :, 2, g), view(d, :, :, 3, g))
    end
    return LWGptGroups(
        ntuple(_ -> TwoStream(grid_params), ngroups),
        ntuple(_ -> SourceLW2Str(grid_params; params), ngroups),
        ntuple(_ -> FluxLW(grid_params), ngroups),
        ntuple(_ -> DA{Bool}(undef, nlay, ncol), ngroups),
        fluxes(acc_data),
        fluxes(acc_clear_data),
        acc_data,
        acc_clear_data,
        DA{Int32}(undef, ncol, ngroups),
    )
end

n_gpt_groups(::Nothing) = 1
n_gpt_groups(::LWGptGroups{N}) where {N} = N
