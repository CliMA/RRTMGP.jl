#####
##### Wide workspaces for a g-point-parallel solve.
#####

"""
    LWGptGroups(grid_params::RRTMGPGridParams, ngroups; params)

Workspaces for solving a column's longwave g-point loop in `ngroups`
independent pieces, all in one kernel launch.

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
struct LWGptGroups{OP, SL, FB, MK, FX, A3, CT}
    ngroups::Int
    "columns per group, so `wcol = gcol + (grp - 1) * ncol`"
    ncol::Int
    op::OP
    src::SL
    fluxb::FB
    mask::MK
    "all-sky partial accumulator over the wide columns, a view into `acc_data`"
    acc::FX
    "clear-sky partial accumulator, a view into `acc_clear_data`"
    acc_clear::FX
    "`(ncol * ngroups, nlev, 3)`; the field axis is up, down, net as in `FluxLW`"
    acc_data::A3
    acc_clear_data::A3
    "`(ncol * ngroups)`: cloudy g-points seen by each group, summed for the cover"
    cloudy::CT
end
Adapt.@adapt_structure LWGptGroups

function LWGptGroups(grid_params::RRTMGPGridParams, ngroups::Int; params)
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
    acc_data = DA{FT}(undef, ncolw, nlev, 3)
    acc_clear_data = DA{FT}(undef, ncolw, nlev, 3)
    # `view(d, :, :, k)` is a contiguous block, so the accumulators index as
    # cheaply as the plain arrays the ungrouped path uses.
    fluxes(d) = FluxLW(view(d, :, :, 1), view(d, :, :, 2), view(d, :, :, 3))
    return LWGptGroups(
        ngroups,
        ncol,
        TwoStream(wide),
        SourceLW2Str(wide; params),
        FluxLW(wide),
        DA{Bool}(undef, nlay, ncolw),
        fluxes(acc_data),
        fluxes(acc_clear_data),
        acc_data,
        acc_clear_data,
        DA{Int32}(undef, ncolw),
    )
end

n_gpt_groups(::Nothing) = 1
n_gpt_groups(g::LWGptGroups) = g.ngroups
