#####
##### One grid over `ncol * ngroups` threads, for either band.
#####

# Both halves of a round launch narrower than the solve's usual 256 threads,
# so that a capped register count reaches more warps per multiprocessor than
# the same count would at 256.
const GROUP_BLOCK_THREADS = 128

# The register budget each half of a round is compiled to. Shedding
# registers off these kernels buys warps, and the warps pay for
# themselves: each count below is the one its own kernel measured
# fastest at, over an interleaved sweep four occupancy steps deep.
# The budgets are per kernel because the two halves compose -- a
# matched pair landed within a few hundredths of a point of the sum
# of its two single-kernel effects at every cap measured -- so a band
# whose halves bottom out at different counts takes both. The sweep
# is a profiling run against a coupled model and lives outside this
# repo; the counts it chose:
# lw optics: 156 registers naturally, -11.41% at 96.
# lw sweeps: 140 registers naturally, -8.02% at 72.
# sw optics: 124 registers naturally, -4.65% at 96.
# sw sweeps: 117 registers naturally, -5.51% at 64.
const LW_OPTICS_MAXREGS = 96
const LW_SWEEPS_MAXREGS = 72
const SW_OPTICS_MAXREGS = 96
const SW_SWEEPS_MAXREGS = 64

@inline _optics_maxregs(::LWGptGroups) = LW_OPTICS_MAXREGS
@inline _optics_maxregs(::SWGptGroups) = SW_OPTICS_MAXREGS
@inline _sweeps_maxregs(::LWGptGroups) = LW_SWEEPS_MAXREGS
@inline _sweeps_maxregs(::SWGptGroups) = SW_SWEEPS_MAXREGS

# What the two bands do not share: the kernels a round launches, the source the
# longwave optics writes and the shortwave has none of, and the cloud-cover
# array the reduction fills.
@inline _optics_kernel(::LWGptGroups) = _lw_group_optics_CUDA!
@inline _optics_kernel(::SWGptGroups) = _sw_group_optics_CUDA!
@inline _sweeps_kernel(::LWGptGroups) = _lw_group_sweeps_CUDA!
@inline _sweeps_kernel(::SWGptGroups) = _sw_group_sweeps_CUDA!
@inline _reduce_kernel(::LWGptGroups) = _reduce_lw_groups_CUDA!
@inline _reduce_kernel(::SWGptGroups) = _reduce_sw_groups_CUDA!
@inline _optics_workspaces(g::LWGptGroups) = (g.op, g.src)
@inline _optics_workspaces(g::SWGptGroups) = (g.op,)
# The shortwave optics takes the zenith angle so it can leave a night column
# alone. The longwave has no such thread: every column emits.
@inline _optics_zenith(::LWGptGroups, bcs) = ()
@inline _optics_zenith(::SWGptGroups, bcs) = (bcs.cos_zenith,)
@inline only_or_nothing(::Tuple{}) = nothing
@inline only_or_nothing(t::Tuple{Any}) = t[1]

# Aerosol mask for every column, hoisted out of the grouped solves of both
# bands, which need it per column rather than per wide column.
function _aero_mask_CUDA!(aerosol_state::AerosolState, ncol)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    gcol ≤ ncol && _compute_aero_mask!(aerosol_state, gcol)
    return nothing
end

# The planes recurse rather than iterate: a loop indexes the tuple with a
# runtime variable, which defeats SROA -- a tuple indexed that way has to exist
# somewhere addressable, so it lands in a stack slot and the sweeps kernels
# spill. Recursion over `Base.tail` indexes by type instead and unrolls, which
# is only a fixed depth because a band accumulates two flux planes or three.
@inline function _zero_planes!(planes::Tuple, wcol, ilev)
    plane = first(planes)
    @inbounds plane[wcol, ilev] = zero(eltype(plane))
    return _zero_planes!(Base.tail(planes), wcol, ilev)
end
@inline _zero_planes!(::Tuple{}, wcol, ilev) = nothing

# The store the optics kernels add aerosols from: the band's own, when the
# workspaces keep one and the solve has aerosols at all.
@inline _aero_store(groups, ::Nothing) = nothing
@inline _aero_store(groups, ::LookUpAerosolMerra) = groups.aero

# What the per-band aerosol optics differ in between the bands: the shortwave
# delta-scales, reports the 550 nm optical depth, and has dark columns.
@inline _aero_delta_scaling(::LWGptGroups) = false
@inline _aero_delta_scaling(::SWGptGroups) = true
@inline _aero_band_550nm(::LWGptGroups, lkp_aero) = nothing
@inline _aero_band_550nm(::SWGptGroups, lkp_aero) = lkp_aero.iband_550nm

# Aerosol optics for every (column, band) pair, into the store the optics
# kernels read. One thread per pair, columns fastest, so a block's writes are
# coalesced. `iband_550nm` is `nothing` in the longwave, which reports no optical
# depth; `cos_zenith` is `nothing` there too, since every column emits.
#
# A dark shortwave column is skipped, as nothing reads its optics, except in the
# 550 nm band: that optical depth is a diagnostic every column has. Writing it
# here, once per column, is also what frees the optics kernel from nominating
# one of a column's concurrent groups as its writer.
function _aerosol_band_optics_CUDA!(
    aero::NamedTuple,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_aero::LookUpAerosolMerra,
    ncol,
    nbnd,
    iband_550nm,
    cos_zenith,
    ::Val{delta_scaling},
) where {delta_scaling}
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i ≤ ncol * nbnd
        gcol = (i - 1) % ncol + 1
        ibnd = (i - 1) ÷ ncol + 1
        collect_aod = !isnothing(iband_550nm) && ibnd == iband_550nm
        @inbounds if isnothing(cos_zenith) ||
                     cos_zenith[gcol] > 0 ||
                     collect_aod
            (; aerosol_state) = as
            Optics.compute_aerosol_band_optics!(
                view(aero.τ, gcol, :, ibnd),
                view(aero.ssa, gcol, :, ibnd),
                view(aero.g, gcol, :, ibnd),
                collect_aod ? view(aerosol_state.aod_sw_ext, gcol) : nothing,
                collect_aod ? view(aerosol_state.aod_sw_sca, gcol) : nothing,
                view(aerosol_state.aero_mask, :, gcol),
                view(aerosol_state.aero_size, :, :, gcol),
                view(aerosol_state.aero_mass, :, :, gcol),
                AtmosphericStates.rel_hum_view(state_cache, as, gcol),
                lookup_aero,
                ibnd,
                iband_550nm;
                delta_scaling,
            )
        end
    end
    return nothing
end

# Zero a group's slice of the partial accumulators. The reduction reads every
# slice whether or not a group had work, so a group with none has to write
# zeros rather than leave whatever the workspace last held. One body for both
# bands: the accumulator is a named tuple of flux planes, two in the longwave
# and three in the shortwave, and this writes whichever it is handed.
@inline function _zero_partials!(acc, acc_clear, wcol, nlev)
    for ilev in 1:nlev
        _zero_planes!(values(acc), wcol, ilev)
        _zero_planes!(values(acc_clear), wcol, ilev)
    end
    return nothing
end

@inline _cld_cover(::Nothing, _) = nothing
@inline _cld_cover(cs::CloudState, ::LWGptGroups) = cs.cld_cover_lw
@inline _cld_cover(cs::CloudState, ::SWGptGroups) = cs.cld_cover_sw

"""
    _solve_both_skies_grouped!(groups, flux, flux_clear, bcs, as, state_cache,
                               lookup, lookup_cld, lookup_aero)

Run the fused both-skies solve over `ncol * n_gpt_groups(groups)` threads, each
taking every `ngroups`-th g-point of one column, then reduce their partial
fluxes. `groups` is the band: `LWGptGroups` or `SWGptGroups`. Each round of that
g-point loop is two kernels, the gas and aerosol optics and then the sweeps; see
`_grouped_rounds!`.

A thread's wide column index is its own id: `wcol = gcol + (grp - 1) * ncol`
selects its slice of every workspace, and consecutive threads are consecutive
columns of one group, so the state reads and the scratch writes stay coalesced.

The gain is waves rather than occupancy. A round's kernels are register-bound
whatever the grid, so what grouping changes is that blocks queue behind them
instead of a multiprocessor idling until the launch's slowest column is done.

The fluxes are **not** bit-for-bit with the ungrouped solve: summing a column's
g-points as `ngroups` interleaved partials is a different order of the same sum.
Which mask each g-point draws is unchanged, since the McICA key indexes a
stateless draw in `(key, column, g-point, layer)`.

Two things are specific to the shortwave. The partials carry a direct-beam
profile as well as up and down, so the reduction sums three. And a night column
is left out of nearly all of a round: the optics kernel returns without filling
`op`, and the sweeps kernel zeroes that thread's slice rather than leaving
whatever the workspace last held, so the reduction sums zeros to the same result.
The exception is the 550 nm aerosol optical depth, which every column has
whether or not the sun is up. It is a read-modify-write on one per-column
address that concurrent groups would interleave, so the kernel nominates one
g-point as its only writer -- and that one thread runs for a dark column too,
being the only thing such a column still owes. It is 1 g-point in 224.
"""
function _solve_both_skies_grouped!(
    groups::Union{LWGptGroups, SWGptGroups},
    flux,
    flux_clear,
    bcs,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup,
    lookup_cld::Union{LookUpCld, Nothing},
    lookup_aero::Union{LookUpAerosolMerra, Nothing},
    lookup_precip::Union{LookUpPrecip, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup.band_data.major_gpt2bnd)
    ngroups = n_gpt_groups(groups)
    (; aerosol_state, cloud_state) = as

    # The aerosol mask depends on the column alone, so it is computed once here
    # instead of identically by every group.
    if aerosol_state isa AerosolState
        atx, abx = _configure_threadblock(ncol)
        @cuda always_inline = true threads = (atx) blocks = (abx) _aero_mask_CUDA!(
            aerosol_state,
            ncol,
        )
    end

    # Aerosol optics do not depend on the g-point, so every band's are computed
    # once here and the rounds below add them to each g-point from the store.
    aero = _aero_store(groups, lookup_aero)
    if !isnothing(aero)
        nbnd = size(aero.τ, 3)
        btx, bbx = _configure_threadblock(ncol * nbnd)
        @cuda always_inline = true threads = (btx) blocks = (bbx) _aerosol_band_optics_CUDA!(
            aero,
            as,
            state_cache,
            lookup_aero,
            ncol,
            nbnd,
            _aero_band_550nm(groups, lookup_aero),
            only_or_nothing(_optics_zenith(groups, bcs)),
            Val(_aero_delta_scaling(groups)),
        )
    end

    _grouped_rounds!(
        groups,
        bcs,
        as,
        state_cache,
        lookup,
        lookup_cld,
        lookup_aero,
        lookup_precip,
        nlay,
        ncol,
        ngroups,
        n_gpt,
    )

    rtx, rbx = _configure_threadblock(ncol * nlev)
    reduce_args = (
        flux,
        flux_clear,
        groups.acc,
        groups.acc_clear,
        groups.cloudy,
        _cld_cover(cloud_state, groups),
        ngroups,
        n_gpt,
        ncol,
        nlev,
    )
    reduce! = _reduce_kernel(groups)
    @cuda always_inline = true threads = (rtx) blocks = (rbx) reduce!(
        reduce_args...,
    )
    return nothing
end

# One launch of a round's half, at the band's block shape and register
# budget. `maxregs` is a compiler keyword, so a round's two halves compile
# once each at their own budget and the loop below launches what it is given.
@inline function _launch_round_half!(f, args, tx, bx, budget)
    @cuda always_inline = true maxregs = budget threads = (tx) blocks = (bx) f(
        args...,
    )
    return nothing
end

# One round of the grouped g-point loop is two launches: every thread's gas and
# aerosol optics, then every thread's sweeps over the optics now in `op`. The
# round index replaces a thread's own `for igpt in grp:ngroups:n_gpt`, so round
# `r` is that loop's `r`th iteration for every group at once.
function _grouped_rounds!(
    groups::Union{LWGptGroups, SWGptGroups},
    bcs,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup,
    lookup_cld::Union{LookUpCld, Nothing},
    lookup_aero::Union{LookUpAerosolMerra, Nothing},
    lookup_precip::Union{LookUpPrecip, Nothing},
    nlay,
    ncol,
    ngroups,
    n_gpt,
)
    nwide = ncol * ngroups
    tx, bx = _configure_threadblock(GROUP_BLOCK_THREADS, nwide)
    optics! = _optics_kernel(groups)
    sweeps! = _sweeps_kernel(groups)
    optics_maxregs = _optics_maxregs(groups)
    sweeps_maxregs = _sweeps_maxregs(groups)
    for round in 1:cld(n_gpt, ngroups)
        optics_args = (
            _optics_workspaces(groups)...,
            as,
            state_cache,
            lookup,
            lookup_aero,
            _aero_store(groups, lookup_aero),
            ncol,
            ngroups,
            n_gpt,
            round,
            _optics_zenith(groups, bcs)...,
        )
        _launch_round_half!(optics!, optics_args, tx, bx, optics_maxregs)
        sweeps_args = (
            groups.fluxb,
            groups.acc,
            groups.acc_clear,
            groups.src,
            bcs,
            groups.op,
            groups.mask,
            groups.cloudy,
            nlay,
            ncol,
            as,
            lookup,
            lookup_cld,
            lookup_precip,
            ngroups,
            n_gpt,
            round,
        )
        _launch_round_half!(sweeps!, sweeps_args, tx, bx, sweeps_maxregs)
    end
    return nothing
end
