#####
##### Wide workspaces for a g-point-parallel solve.
#####

"""
    LWGptGroups(grid_params::RRTMGPGridParams, ngroups; params)

Workspaces for solving a column's longwave g-point loop in `ngroups`
independent pieces, covered by one grid rather than by one column's serial
loop.

At one thread per column the solve fills the device in a single wave of blocks,
so its tail and the imbalance between columns cannot be absorbed. Splitting the
g-point loop multiplies the blocks instead, and the scheduler refills a
multiprocessor as each retires.

Every workspace the solve indexes by column is widened to `ncol * ngroups`, so a
thread's `(column, group)` pair selects its slice by one index and the kernel
body is the ungrouped one with that index substituted wherever it touches
scratch. Nothing is allocated or walked per group: the reduction recovers the
group from the wide column index.

`ngroups > 1` is implemented for CUDA. The CPU driver ignores these workspaces
and solves the whole g-point loop per column, which is correct either way.
"""
struct LWGptGroups{OP, SL, FB, MK, FX, CT}
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
end
Adapt.@adapt_structure LWGptGroups

function LWGptGroups(
    grid_params::RRTMGPGridParams,
    ngroups::Int;
    params,
)
    ngroups ≥ 1 || error("ngroups must be at least 1, got $ngroups")
    (; context, ncol, nlay, isothermal_boundary_layer) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    nlev = nlay + 1
    ncolw = ncol * ngroups
    # The inner constructor, because `nlay` here is already the total and the
    # keyword one would add the isothermal layer a second time.
    wide = RRTMGPGridParams{FT, typeof(context)}(
        context,
        nlay,
        ncolw,
        isothermal_boundary_layer,
    )
    # Up and down only: the reduction recomputes net from them, so a third
    # plane would be allocated per wide column and never written. The field
    # names match `FluxLW` because that is what these accumulate, but the type
    # does not, which is the point -- a `FluxLW` would carry the net plane.
    acc_data = DA{FT}(undef, ncolw, nlev, 2)
    acc_clear_data = DA{FT}(undef, ncolw, nlev, 2)
    # `view(d, :, :, k)` is a contiguous block, so the accumulators index as
    # cheaply as the plain arrays the ungrouped path uses.
    fluxes(d) = (flux_up = view(d, :, :, 1), flux_dn = view(d, :, :, 2))
    return LWGptGroups(
        ngroups,
        TwoStream(wide),
        SourceLW2Str(wide; params),
        FluxLW(wide),
        DA{Bool}(undef, nlay, ncolw),
        fluxes(acc_data),
        fluxes(acc_clear_data),
        DA{Int32}(undef, ncolw),
    )
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
struct SWGptGroups{OP, SS, FB, MK, FX, CT}
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
end
Adapt.@adapt_structure SWGptGroups

function SWGptGroups(
    grid_params::RRTMGPGridParams,
    ngroups::Int;
)
    ngroups >= 1 || error("ngroups must be at least 1, got $ngroups")
    (; context, ncol, nlay, isothermal_boundary_layer) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    nlev = nlay + 1
    ncolw = ncol * ngroups
    # The inner constructor, because `nlay` here is already the total and the
    # keyword one would add the isothermal layer a second time.
    wide = RRTMGPGridParams{FT, typeof(context)}(
        context,
        nlay,
        ncolw,
        isothermal_boundary_layer,
    )
    # Up, down and direct: the reduction recomputes net from up and down, so a
    # net plane would never be written (see `LWGptGroups`).
    acc_data = DA{FT}(undef, ncolw, nlev, 3)
    acc_clear_data = DA{FT}(undef, ncolw, nlev, 3)
    fluxes(d) = (
        flux_up = view(d, :, :, 1),
        flux_dn = view(d, :, :, 2),
        flux_dn_dir = view(d, :, :, 3),
    )
    return SWGptGroups(
        ngroups,
        TwoStream(wide),
        SourceSW2Str(wide),
        FluxSW(wide),
        DA{Bool}(undef, nlay, ncolw),
        fluxes(acc_data),
        fluxes(acc_clear_data),
        DA{Int32}(undef, ncolw),
    )
end

n_gpt_groups(g::SWGptGroups) = g.ngroups
