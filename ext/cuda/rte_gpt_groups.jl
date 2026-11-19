#####
##### One grid over `ncol * ngroups` threads, for either band.
#####

# Both halves of a round launch narrower than the solve's usual 256 threads, so
# that their register counts (156 in the optics, 140 in the longwave sweeps)
# reach 12 warps per multiprocessor rather than the 8 they get at 256.
const GROUP_BLOCK_THREADS = 128

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

# Aerosol mask for every column, hoisted out of the grouped solves of both
# bands, which need it per column rather than per wide column.
function _aero_mask_CUDA!(aerosol_state::AerosolState, ncol)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    gcol ≤ ncol && _compute_aero_mask!(aerosol_state, gcol)
    return nothing
end

# Zero a group's slice of the partial accumulators. The reduction reads every
# slice whether or not a group had work, so a group with none has to write
# zeros rather than leave whatever the workspace last held. One body for both
# bands: the accumulator is a named tuple of flux planes, two in the longwave
# and three in the shortwave, and this writes whichever it is handed.
@inline function _zero_partials!(acc, acc_clear, wcol, nlev)
    FT = eltype(acc.flux_up)
    @inbounds for ilev in 1:nlev
        for plane in acc
            plane[wcol, ilev] = zero(FT)
        end
        for plane in acc_clear
            plane[wcol, ilev] = zero(FT)
        end
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
never sweeps, so the sweeps kernel zeroes that thread's slice rather than leaving
whatever the workspace last held, and the reduction sums zeros to the same
result. The 550 nm aerosol optical depth **is** bit-for-bit, and needs care to
be: it is a read-modify-write on one per-column address that concurrent groups
would interleave, so the kernel nominates one g-point as its only writer.
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
    for round in 1:cld(n_gpt, ngroups)
        optics_args = (
            _optics_workspaces(groups)...,
            as,
            state_cache,
            lookup,
            lookup_aero,
            ncol,
            ngroups,
            n_gpt,
            round,
        )
        @cuda always_inline = true threads = (tx) blocks = (bx) optics!(
            optics_args...,
        )
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
        @cuda always_inline = true threads = (tx) blocks = (bx) sweeps!(
            sweeps_args...,
        )
    end
    return nothing
end
