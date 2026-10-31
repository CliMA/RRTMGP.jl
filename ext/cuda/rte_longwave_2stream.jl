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

# One stream per group plus the events that fork and join them, created on first
# use and kept, since the solve runs at every radiation step. The cache is tied
# to whatever context was current when it was filled, which is why it assumes
# what the rest of this driver assumes: one device, one solve in flight.
const _GPT_STREAMS = CUDA.CuStream[]
const _GPT_JOINS = CUDA.CuEvent[]
const _GPT_FORK = CUDA.CuEvent[]

function _gpt_concurrency(ngroups)
    while length(_GPT_STREAMS) < ngroups
        push!(_GPT_STREAMS, CUDA.CuStream())
        push!(_GPT_JOINS, CUDA.CuEvent())
    end
    isempty(_GPT_FORK) &&
        push!(_GPT_FORK, CUDA.CuEvent())
    return _GPT_STREAMS, _GPT_JOINS, _GPT_FORK[1]
end

"""
    _solve_lw_both_skies_grouped!(groups, flux_lw, flux_lw_clear, bcs_lw, as, ...)

Run the fused longwave solve as `n_gpt_groups(groups)` concurrent kernels, each
taking every `ngroups`-th g-point, and reduce their partial fluxes.

The point is resident parallelism: one thread per column is 96 blocks on 108
SMs, a single wave with no second wave to hide its latency behind. The groups
are separate launches on separate streams rather than one wider grid because
that way each group's workspaces arrive as kernel arguments already pointing at
its own memory, so the solve kernel compiles to exactly what the ungrouped one
does -- a kernel this close to the register ceiling cannot afford to compute
workspace offsets per thread.

Fork and join through events rather than `synchronize`: the host must not block
here, since it is otherwise free to queue the next steps' work.

The fluxes are **not** bit-for-bit with the ungrouped solve. Summing a column's
g-points in 16 interleaved partials and then adding the partials is a different
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
    tx, bx = _configure_threadblock(ncol)
    (; aerosol_state, cloud_state) = as

    # The aerosol mask depends on the column alone, so it is computed once here
    # instead of identically by every group.
    if aerosol_state isa AerosolState
        @cuda always_inline = true threads = (tx) blocks = (bx) _lw_aero_mask_CUDA!(
            aerosol_state,
            ncol,
        )
    end

    streams, joins, fork = _gpt_concurrency(ngroups)
    CUDA.record(fork)
    for g in 1:ngroups
        s = streams[g]
        CUDA.wait(fork, s)
        args = (
            groups.fluxb[g],
            groups.acc[g],
            groups.acc_clear[g],
            groups.src[g],
            bcs_lw,
            groups.op[g],
            groups.mask[g],
            view(groups.cloudy, :, g),
            nlay,
            ncol,
            as,
            state_cache,
            lookup_lw,
            lookup_lw_cld,
            lookup_lw_aero,
            g,
            ngroups,
        )
        @cuda always_inline = true threads = (tx) blocks = (bx) stream = s rte_lw_2stream_solve_group_CUDA!(
            args...,
        )
        CUDA.record(joins[g], s)
    end
    for g in 1:ngroups
        CUDA.wait(joins[g])
    end

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

# Aerosol mask for every column, hoisted out of the group kernels.
function _lw_aero_mask_CUDA!(aerosol_state::AerosolState, ncol)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    gcol ≤ ncol && _compute_aero_mask!(aerosol_state, gcol)
    return nothing
end

# One group's share of the fused longwave g-point loop. The per-column body is
# the same `lw_2stream_gpt_col_both_skies!` the ungrouped kernel calls, handed
# this group's workspaces and told which g-point starts its stride, so the
# accumulators are assigned on the first one and added to afterwards.
function rte_lw_2stream_solve_group_CUDA!(
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
    grp,
    ngroups,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        FT = eltype(acc.flux_up)
        n_cloudy_gpts = 0
        if grp > n_gpt
            # More groups than g-points: this one has no work, and the reduction
            # reads its slice regardless.
            @inbounds for ilev in 1:nlev
                acc.flux_up[gcol, ilev] = FT(0)
                acc.flux_dn[gcol, ilev] = FT(0)
                acc_clear.flux_up[gcol, ilev] = FT(0)
                acc_clear.flux_dn[gcol, ilev] = FT(0)
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
                )
                n_cloudy_gpts += cloudy ? 1 : 0
            end
        end
        @inbounds cloudy_count[gcol] = Int32(n_cloudy_gpts)
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
            up += acc_data[gcol, ilev, 1, g]
            dn += acc_data[gcol, ilev, 2, g]
            cup += acc_clear_data[gcol, ilev, 1, g]
            cdn += acc_clear_data[gcol, ilev, 2, g]
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
                    n += cloudy_count[gcol, g]
                end
                cld_cover_lw[gcol] = FT(n) / n_gpt
            end
        end
    end
    return nothing
end
