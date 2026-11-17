#####
##### Wide workspaces for a g-point-parallel solve.
#####

"""
    LWGptGroups(grid_params::RRTMGPGridParams, ngroups; params)

Workspaces for solving a column's longwave g-point loop in `ngroups`
independent pieces, covered by one grid rather than by one column's serial
loop.

At one thread per column the solve's grid is 96 blocks of 256 threads, and the
kernel's 255 registers pin it to one block per SM, so 108 SMs hold 0.89 of that
grid at once. The device is already full. Splitting the g-point loop multiplies
the number of blocks instead: `ncol * ngroups` threads is some 14 waves, and the
scheduler refills an SM as each block retires, which is what absorbs the tail
and the imbalance between columns. This is the restructure priced in section 4l
of the profiling study.

Every workspace the solve indexes by column is widened to `ncol * ngroups`
columns, so a thread's `(column, group)` pair selects its own slice by a single
index and the kernel body is the ungrouped one with that index substituted for
the column wherever it touches scratch. Nothing is per-group-allocated and
nothing has a group axis to walk: the reduction recovers the group from the
wide column index.

`ngroups > 1` is implemented for CUDA. The CPU driver ignores these workspaces
and solves the whole g-point loop per column, which is correct either way.
"""
struct LWGptGroups{OP, SL, FB, MK, FX, CT}
    ngroups::Int
    """
    Solve each round of the g-point loop as two kernels, the gas and aerosol
    optics and then the sweeps, instead of one fused kernel.

    The fused kernel needs 255 registers, which pins it to one block per SM.
    Compiled on its own the optics phase needs 156 in the longwave and 128 in
    the shortwave with 32 bytes of stack, so it reaches 12 and 16 warps per SM
    where the fused kernel reaches 8. `op` and `src` are already device arrays
    that both halves address by wide column, so the phases exchange nothing a
    kernel boundary cannot carry and the split result is bit-for-bit.
    """
    split_optics::Bool
    "columns per group, so `wcol = gcol + (grp - 1) * ncol`"
    ncol::Int
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
    split_optics::Bool = true,
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
        split_optics,
        ncol,
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

splits_optics(g) = g.split_optics

"""
    SWGptGroups(grid_params::RRTMGPGridParams, ngroups)

Workspaces for solving a column's shortwave g-point loop in `ngroups`
independent pieces, covered by one grid rather than by one column's serial
loop.

[`LWGptGroups`](@ref) carries the reasoning, and it holds here unchanged: the
shortwave solve's grid is also 96 blocks of 256 threads pinned to one block per
SM by 255 registers, so 108 SMs hold 0.89 of it and the parallelism has to come
from more blocks rather than from more occupancy.

Two things differ from the longwave. `FluxSW` carries a direct-beam profile, so
a group accumulates three of them rather than two. And the shortwave sweeps only
where the sun is up, so every group of a night column exits its loop at once --
at fourteen waves an SM that retires such a block is handed another, which is
the load imbalance section 4l of the profiling study measured as unobtainable
while the grid was one wave.
"""
struct SWGptGroups{OP, SS, FB, MK, FX, CT}
    ngroups::Int
    "split the optics out of the sweeps; see [`LWGptGroups`](@ref)"
    split_optics::Bool
    "columns per group, so `wcol = gcol + (grp - 1) * ncol`"
    ncol::Int
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
    split_optics::Bool = true,
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
        split_optics,
        ncol,
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
