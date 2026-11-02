function rte_lw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    args = (flux_lw, src_lw, bcs_lw, op, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_lw_2stream_solve_CUDA!(
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    nlay,
    ncol,
    as::GrayAtmosphericState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    igpt, ibnd = 1, 1
    if gcol ≤ ncol
        compute_optical_props!(op, as, src_lw, gcol)
        rte_lw_2stream!(
            op,
            flux_lw,
            src_lw,
            bcs_lw,
            gcol,
            igpt,
            ibnd,
            nlev,
            ncol,
        )
        compute_net_flux!(flux_lw, gcol, nlev)
    end
    return nothing
end

function rte_lw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    band_flux,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    set_band_flux_to_zero!(band_flux)
    tx, bx = _configure_threadblock(ncol)
    args = (
        flux,
        flux_lw,
        band_flux,
        src_lw,
        bcs_lw,
        op,
        nlay,
        ncol,
        as,
        state_cache,
        lookup_lw,
        lookup_lw_cld,
        lookup_lw_aero,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_lw_2stream_solve_CUDA!(
    flux::FluxLW,
    flux_lw::FluxLW,
    band_flux,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    nlay,
    ncol,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        (; cloud_state, aerosol_state) = as
        FT = eltype(flux_lw.flux_up)
        _compute_aero_mask!(aerosol_state, gcol)
        n_cloudy_gpts = 0  # thread-local counter for LW cloud cover
        @inbounds for igpt in 1:n_gpt
            cloudy = lw_2stream_gpt_col!(
                igpt,
                gcol,
                flux,
                flux_lw,
                band_flux,
                src_lw,
                bcs_lw,
                op,
                as,
                state_cache,
                lookup_lw,
                lookup_lw_cld,
                lookup_lw_aero,
                major_gpt2bnd[igpt],
                nlev,
                ncol,
            )
            n_cloudy_gpts += cloudy ? 1 : 0
        end
        @inbounds begin
            compute_net_flux!(flux_lw, gcol, nlev)
            # write out LW cloud cover
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_lw)
                cloud_state.cld_cover_lw[gcol] = FT(n_cloudy_gpts) / n_gpt
            end
        end
    end
    return nothing
end

"""
    rte_lw_2stream_solve_both_skies!(device::CUDADevice, ...)

Host side of the fused longwave solve: one thread per column, launching
`rte_lw_2stream_solve_both_skies_CUDA!`. Spectral only -- the identically
named gray methods above take a `GrayAtmosphericState` and have no clouds.
"""
function rte_lw_2stream_solve_both_skies!(
    device::ClimaComms.CUDADevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    flux_lw_clear::FluxLW,
    band_flux,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
    gpt_groups = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    set_band_flux_to_zero!(band_flux)
    # The g-point loop splits only when no per-band fluxes are being kept: those
    # accumulate per g-point too, so they would need per-group partials of their
    # own and a wider reduction.
    if !isnothing(gpt_groups) && isnothing(band_flux)
        return _solve_lw_both_skies_grouped!(
            gpt_groups,
            flux_lw,
            flux_lw_clear,
            bcs_lw,
            as,
            state_cache,
            lookup_lw,
            lookup_lw_cld,
            lookup_lw_aero,
        )
    end
    args = (
        flux,
        flux_lw,
        flux_lw_clear,
        band_flux,
        src_lw,
        bcs_lw,
        op,
        nlay,
        ncol,
        as,
        state_cache,
        lookup_lw,
        lookup_lw_cld,
        lookup_lw_aero,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_solve_both_skies_CUDA!(
        args...,
    )
    return nothing
end

"""
    rte_lw_2stream_solve_both_skies_CUDA!(...)

The fused longwave kernel: one thread per column, looping g-points. Per g-point
it computes the optics once, sweeps and accumulates the clear sky, then adds the
cloud increment to the same optics and sweeps again. Per-column body is
`lw_2stream_gpt_col_both_skies!` in `src/rte/longwave_2stream.jl`, shared with the
CPU driver.
"""
# Duplicated for the CPU in src/rte/longwave_2stream.jl: the cloud-cover tally
# and net-flux finalize are in both. Change both.
function rte_lw_2stream_solve_both_skies_CUDA!(
    flux::FluxLW,
    flux_lw::FluxLW,
    flux_lw_clear::FluxLW,
    band_flux,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    nlay,
    ncol,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        (; cloud_state, aerosol_state) = as
        FT = eltype(flux_lw.flux_up)
        _compute_aero_mask!(aerosol_state, gcol)
        n_cloudy_gpts = 0  # thread-local counter for LW cloud cover
        @inbounds for igpt in 1:n_gpt
            cloudy = lw_2stream_gpt_col_both_skies!(
                igpt,
                gcol,
                flux,
                flux_lw,
                flux_lw_clear,
                band_flux,
                src_lw,
                bcs_lw,
                op,
                as,
                state_cache,
                lookup_lw,
                lookup_lw_cld,
                lookup_lw_aero,
                major_gpt2bnd[igpt],
                nlev,
                ncol,
            )
            n_cloudy_gpts += cloudy ? 1 : 0
        end
        @inbounds begin
            compute_net_flux!(flux_lw, gcol, nlev)
            compute_net_flux!(flux_lw_clear, gcol, nlev)
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_lw)
                cloud_state.cld_cover_lw[gcol] = FT(n_cloudy_gpts) / n_gpt
            end
        end
    end
    return nothing
end

#####
##### The g-point-grouped fused longwave solve.
#####

"""
    _solve_lw_both_skies_grouped!(groups, flux_lw, flux_lw_clear, bcs_lw, as, ...)

Run the fused longwave solve as one kernel over `ncol * n_gpt_groups(groups)`
threads, each taking every `ngroups`-th g-point of one column, then reduce
their partial fluxes.

The gain is waves, not occupancy. The kernel's 255 registers pin it to one
block per SM either way, so the grouped grid has exactly as much work resident
at any instant as the ungrouped one; what changes is that blocks are now queued
behind it, so an SM that finishes a cheap column is handed another instead of
idling until the launch's slowest column is done.

One grid rather than one launch per group is the whole point. At one thread per
column the grid is already 0.89 of the device's resident capacity, so separate
concurrent launches cannot co-reside -- measured overlap was 1.97x for sixteen
of them -- and the parallelism has to come from more blocks in the one grid.

A thread's wide column index is its own id: `wcol = gcol + (grp - 1) * ncol`
selects its slice of every workspace, and consecutive threads are consecutive
columns of one group, so the state reads and the scratch writes both stay
coalesced.

The fluxes are **not** bit-for-bit with the ungrouped solve. Summing a column's
g-points as `ngroups` interleaved partials and adding those is a different
order of the same floating-point sum. Everything else is unchanged, including
which mask each g-point draws -- the McICA key indexes a stateless draw in
`(key, column, g-point, layer)`, so splitting the loop cannot shift it.
"""
function _solve_lw_both_skies_grouped!(
    groups::LWGptGroups,
    flux_lw::FluxLW,
    flux_lw_clear::FluxLW,
    bcs_lw::LwBCs,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_lw.band_data.major_gpt2bnd)
    ngroups = n_gpt_groups(groups)
    (; aerosol_state, cloud_state) = as

    # The aerosol mask depends on the column alone, so it is computed once here
    # instead of identically by every group.
    if aerosol_state isa AerosolState
        atx, abx = _configure_threadblock(ncol)
        @cuda always_inline = true threads = (atx) blocks = (abx) _lw_aero_mask_CUDA!(
            aerosol_state,
            ncol,
        )
    end

    tx, bx = _configure_threadblock(ncol * ngroups)
    args = (
        groups.fluxb,
        groups.acc,
        groups.acc_clear,
        groups.src,
        bcs_lw,
        groups.op,
        groups.mask,
        groups.cloudy,
        nlay,
        ncol,
        as,
        state_cache,
        lookup_lw,
        lookup_lw_cld,
        lookup_lw_aero,
        ngroups,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_solve_grouped_CUDA!(
        args...,
    )

    rtx, rbx = _configure_threadblock(ncol * nlev)
    @cuda always_inline = true threads = (rtx) blocks = (rbx) _reduce_lw_groups_CUDA!(
        flux_lw,
        flux_lw_clear,
        groups.acc_data,
        groups.acc_clear_data,
        groups.cloudy,
        cloud_state isa CloudState ? cloud_state.cld_cover_lw : nothing,
        ngroups,
        n_gpt,
        ncol,
        nlev,
    )
    return nothing
end

# Aerosol mask for every column, hoisted out of the grouped kernel.
function _lw_aero_mask_CUDA!(aerosol_state::AerosolState, ncol)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    gcol ≤ ncol && _compute_aero_mask!(aerosol_state, gcol)
    return nothing
end

# The fused longwave g-point loop over `(column, group)` pairs. The per-column
# body is the same `lw_2stream_gpt_col_both_skies!` the ungrouped kernel calls,
# handed this thread's wide column index and told which g-point starts its
# stride, so the accumulators are assigned on the first one and added to after.
function rte_lw_2stream_solve_grouped_CUDA!(
    flux::FluxLW,
    acc::FluxLW,
    acc_clear::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    mask,
    cloudy_count,
    nlay,
    ncol,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
    ngroups,
)
    wcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    if wcol ≤ ncol * ngroups
        # Columns vary fastest, so a warp spans 32 consecutive columns of one
        # group and every read it makes is still contiguous.
        gcol = (wcol - 1) % ncol + 1
        grp = (wcol - 1) ÷ ncol + 1
        FT = eltype(acc.flux_up)
        n_cloudy_gpts = 0
        if grp > n_gpt
            # More groups than g-points: this one has no work, and the
            # reduction reads its slice regardless.
            @inbounds for ilev in 1:nlev
                acc.flux_up[wcol, ilev] = FT(0)
                acc.flux_dn[wcol, ilev] = FT(0)
                acc_clear.flux_up[wcol, ilev] = FT(0)
                acc_clear.flux_dn[wcol, ilev] = FT(0)
            end
        else
            @inbounds for igpt in grp:ngroups:n_gpt
                cloudy = lw_2stream_gpt_col_both_skies!(
                    igpt,
                    gcol,
                    flux,
                    acc,
                    acc_clear,
                    nothing,
                    src_lw,
                    bcs_lw,
                    op,
                    as,
                    state_cache,
                    lookup_lw,
                    lookup_lw_cld,
                    lookup_lw_aero,
                    major_gpt2bnd[igpt],
                    nlev,
                    ncol,
                    mask,
                    grp,
                    wcol,
                )
                n_cloudy_gpts += cloudy ? 1 : 0
            end
        end
        @inbounds cloudy_count[wcol] = Int32(n_cloudy_gpts)
    end
    return nothing
end

# Sum the groups' partial fluxes into the broadband ones, finish the net flux,
# and total the cloud-cover tally. One thread per (column, level).
function _reduce_lw_groups_CUDA!(
    flux_lw::FluxLW,
    flux_lw_clear::FluxLW,
    acc_data,
    acc_clear_data,
    cloudy_count,
    cld_cover_lw,
    ngroups,
    n_gpt,
    ncol,
    nlev,
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i ≤ ncol * nlev
        gcol = (i - 1) % ncol + 1
        ilev = (i - 1) ÷ ncol + 1
        FT = eltype(flux_lw.flux_up)
        up, dn, cup, cdn = FT(0), FT(0), FT(0), FT(0)
        @inbounds for g in 1:ngroups
            wcol = gcol + (g - 1) * ncol
            up += acc_data[wcol, ilev, 1]
            dn += acc_data[wcol, ilev, 2]
            cup += acc_clear_data[wcol, ilev, 1]
            cdn += acc_clear_data[wcol, ilev, 2]
        end
        @inbounds begin
            flux_lw.flux_up[gcol, ilev] = up
            flux_lw.flux_dn[gcol, ilev] = dn
            flux_lw.flux_net[gcol, ilev] = up - dn
            flux_lw_clear.flux_up[gcol, ilev] = cup
            flux_lw_clear.flux_dn[gcol, ilev] = cdn
            flux_lw_clear.flux_net[gcol, ilev] = cup - cdn
            if ilev == 1 && !isnothing(cld_cover_lw)
                n = 0
                for g in 1:ngroups
                    n += cloudy_count[gcol + (g - 1) * ncol]
                end
                cld_cover_lw[gcol] = FT(n) / n_gpt
            end
        end
    end
    return nothing
end
