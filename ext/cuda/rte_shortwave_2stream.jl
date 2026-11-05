function rte_sw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    args = (flux_sw, op, bcs_sw, src_sw, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay,
    ncol,
    as::GrayAtmosphericState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt, ibnd = 1, 1, 1
    FT = eltype(bcs_sw.cos_zenith)
    solar_frac = FT(1)
    if gcol ≤ ncol
        μ₀ = bcs_sw.cos_zenith[gcol]
        @inbounds begin
            compute_optical_props!(op, as, gcol)
            # call shortwave rte solver
            rte_sw_2stream!(
                op,
                src_sw,
                bcs_sw,
                flux_sw,
                solar_frac,
                igpt,
                n_gpt,
                ibnd,
                nlev,
                gcol,
            )
            compute_net_flux!(flux_sw, gcol, nlev)
        end
        if μ₀ ≤ 0 # zero out columns with zenith angle ≥ π/2
            set_flux_to_zero!(flux_sw, gcol, nlev)
        end
    end
    return nothing
end

function rte_sw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    set_band_flux_to_zero!(band_flux)
    tx, bx = _configure_threadblock(ncol)
    args = (
        flux,
        flux_sw,
        band_flux,
        op,
        bcs_sw,
        src_sw,
        nlay,
        ncol,
        as,
        state_cache,
        lookup_sw,
        lookup_sw_cld,
        lookup_sw_aero,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux::FluxSW,
    flux_sw::FluxSW,
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay,
    ncol,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.band_data.major_gpt2bnd)
    if gcol ≤ ncol
        FT = eltype(flux_sw.flux_up)
        (; cloud_state, aerosol_state) = as
        μ₀ = @inbounds bcs_sw.cos_zenith[gcol]
        n_cloudy_gpts = 0  # thread-local counter for cloud cover
        @inbounds begin
            _compute_aero_mask!(aerosol_state, gcol)
            for igpt in 1:n_gpt
                cloudy = sw_2stream_gpt_col!(
                    igpt,
                    gcol,
                    flux,
                    flux_sw,
                    band_flux,
                    op,
                    bcs_sw,
                    src_sw,
                    as,
                    state_cache,
                    lookup_sw,
                    lookup_sw_cld,
                    lookup_sw_aero,
                    μ₀,
                    lookup_sw.band_data.major_gpt2bnd[igpt],
                    n_gpt,
                    nlev,
                )
                n_cloudy_gpts += cloudy ? 1 : 0
            end
            if μ₀ ≤ 0 # zero out columns with zenith angle ≥ π/2
                set_flux_to_zero!(flux_sw, gcol, nlev)
            else
                compute_net_flux!(flux_sw, gcol, nlev)

            end
            # write out SW cloud cover
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_sw)
                cloud_state.cld_cover_sw[gcol] = FT(n_cloudy_gpts) / n_gpt
            end
        end
    end
    return nothing
end

"""
    rte_sw_2stream_solve_both_skies!(device::CUDADevice, ...)

Host side of the fused shortwave solve: one thread per column, launching
`rte_sw_2stream_solve_both_skies_CUDA!`. Spectral only -- the identically
named gray methods above take a `GrayAtmosphericState` and have no clouds.
"""
function rte_sw_2stream_solve_both_skies!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    flux_sw_clear::FluxSW,
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
    gpt_groups = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    set_band_flux_to_zero!(band_flux)
    # The g-point loop splits only when no per-band fluxes are being kept: those
    # accumulate per g-point too, so they would need per-group partials of their
    # own and a wider reduction.
    if !isnothing(gpt_groups) && isnothing(band_flux)
        return _solve_sw_both_skies_grouped!(
            gpt_groups,
            flux_sw,
            flux_sw_clear,
            bcs_sw,
            as,
            state_cache,
            lookup_sw,
            lookup_sw_cld,
            lookup_sw_aero,
        )
    end
    args = (
        flux,
        flux_sw,
        flux_sw_clear,
        band_flux,
        op,
        bcs_sw,
        src_sw,
        nlay,
        ncol,
        as,
        state_cache,
        lookup_sw,
        lookup_sw_cld,
        lookup_sw_aero,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_both_skies_CUDA!(
        args...,
    )
    return nothing
end

"""
    rte_sw_2stream_solve_both_skies_CUDA!(...)

The fused shortwave kernel: one thread per column, looping g-points. Per g-point
it computes the optics once, sweeps and accumulates the clear sky, then adds the
cloud increment to the same optics and sweeps again. Per-column body is
`sw_2stream_gpt_col_both_skies!` in `src/rte/shortwave_2stream.jl`, shared with the
CPU driver.
"""
# Duplicated for the CPU in src/rte/shortwave_2stream.jl: the cloud-cover tally,
# net-flux finalize and night-column zeroing are in both. Change both.
function rte_sw_2stream_solve_both_skies_CUDA!(
    flux::FluxSW,
    flux_sw::FluxSW,
    flux_sw_clear::FluxSW,
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay,
    ncol,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.band_data.major_gpt2bnd)
    if gcol ≤ ncol
        FT = eltype(flux_sw.flux_up)
        (; cloud_state, aerosol_state) = as
        μ₀ = @inbounds bcs_sw.cos_zenith[gcol]
        n_cloudy_gpts = 0  # thread-local counter for cloud cover
        @inbounds begin
            _compute_aero_mask!(aerosol_state, gcol)
            for igpt in 1:n_gpt
                cloudy = sw_2stream_gpt_col_both_skies!(
                    igpt,
                    gcol,
                    flux,
                    flux_sw,
                    flux_sw_clear,
                    band_flux,
                    op,
                    bcs_sw,
                    src_sw,
                    as,
                    state_cache,
                    lookup_sw,
                    lookup_sw_cld,
                    lookup_sw_aero,
                    μ₀,
                    lookup_sw.band_data.major_gpt2bnd[igpt],
                    n_gpt,
                    nlev,
                )
                n_cloudy_gpts += cloudy ? 1 : 0
            end
            # Night columns never swept, so zero both skies rather than leaving
            # whatever the accumulators last held
            if μ₀ ≤ 0
                set_flux_to_zero!(flux_sw, gcol, nlev)
                set_flux_to_zero!(flux_sw_clear, gcol, nlev)
            else
                compute_net_flux!(flux_sw, gcol, nlev)
                compute_net_flux!(flux_sw_clear, gcol, nlev)
            end
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_sw)
                cloud_state.cld_cover_sw[gcol] = FT(n_cloudy_gpts) / n_gpt
            end
        end
    end
    return nothing
end

"""
    _solve_sw_both_skies_grouped!(groups, flux_sw, flux_sw_clear, bcs_sw, as, ...)

Run the fused shortwave solve as one kernel over `ncol * n_gpt_groups(groups)`
threads, each taking every `ngroups`-th g-point of one column, then reduce their
partial fluxes. The longwave counterpart `_solve_lw_both_skies_grouped!` carries
the reasoning for why this is one grid rather than one launch per group.

Two things are specific to the shortwave. The partials carry a direct-beam
profile as well as up and down, so the reduction sums three. And a night column
never sweeps: `sw_2stream_gpt_col_both_skies!` leaves its accumulators untouched
where `μ₀ ≤ 0`, which the ungrouped kernel then overwrites with zeros after the
loop, so here the kernel zeroes that thread's slice instead and the reduction
sums zeros to the same result.

The work done is the ungrouped work. Optics are still computed for night columns
and the cloud masks are still drawn for them, because the ungrouped kernel does
both and the cover tally depends on the draws. Skipping the optics where the sun
is down is a separate change, and section 4l of the profiling study is why it is
worth making once this one is in: it measured that skip at -2.63% and diagnosed
the ceiling as the solve being one resident wave, which this splits into some
fourteen.

The fluxes are **not** bit-for-bit with the ungrouped solve, for the same reason
as the longwave: a column's g-points are summed as `ngroups` interleaved
partials. The 550 nm aerosol optical depth **is** bit-for-bit, and needs care to
be: it is zeroed and re-accumulated per g-point, so it is a read-modify-write on
one per-column address that concurrent groups would interleave. Every g-point of
the band yields the same total and the ungrouped loop lets the last one stand, so
the kernel nominates that same g-point as the only writer.
"""
function _solve_sw_both_skies_grouped!(
    groups::SWGptGroups,
    flux_sw::FluxSW,
    flux_sw_clear::FluxSW,
    bcs_sw::SwBCs,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_sw.band_data.major_gpt2bnd)
    ngroups = n_gpt_groups(groups)
    (; aerosol_state, cloud_state) = as

    # The aerosol mask depends on the column alone, so it is computed once here
    # instead of identically by every group.
    if aerosol_state isa AerosolState
        atx, abx = _configure_threadblock(ncol)
        @cuda always_inline = true threads = (atx) blocks = (abx) _sw_aero_mask_CUDA!(
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
        bcs_sw,
        groups.op,
        groups.mask,
        groups.cloudy,
        nlay,
        ncol,
        as,
        state_cache,
        lookup_sw,
        lookup_sw_cld,
        lookup_sw_aero,
        ngroups,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_grouped_CUDA!(
        args...,
    )

    rtx, rbx = _configure_threadblock(ncol * nlev)
    @cuda always_inline = true threads = (rtx) blocks = (rbx) _reduce_sw_groups_CUDA!(
        flux_sw,
        flux_sw_clear,
        groups.acc_data,
        groups.acc_clear_data,
        groups.cloudy,
        cloud_state isa CloudState ? cloud_state.cld_cover_sw : nothing,
        ngroups,
        n_gpt,
        ncol,
        nlev,
    )
    return nothing
end

# Aerosol mask for every column, hoisted out of the grouped kernel.
function _sw_aero_mask_CUDA!(aerosol_state::AerosolState, ncol)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    gcol <= ncol && _compute_aero_mask!(aerosol_state, gcol)
    return nothing
end

# The fused shortwave g-point loop over `(column, group)` pairs. The per-column
# body is the same `sw_2stream_gpt_col_both_skies!` the ungrouped kernel calls,
# handed this thread's wide column index and told which g-point starts its
# stride, so the accumulators are assigned on the first one and added to after.
function rte_sw_2stream_solve_grouped_CUDA!(
    flux::FluxSW,
    acc::FluxSW,
    acc_clear::FluxSW,
    src_sw::SourceSW2Str,
    bcs_sw::SwBCs,
    op::TwoStream,
    mask,
    cloudy_count,
    nlay,
    ncol,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
    ngroups,
)
    wcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_sw.band_data
    n_gpt = length(major_gpt2bnd)
    if wcol <= ncol * ngroups
        # Columns vary fastest, so a warp spans 32 consecutive columns of one
        # group and every read it makes is still contiguous.
        gcol = (wcol - 1) % ncol + 1
        grp = (wcol - 1) ÷ ncol + 1
        FT = eltype(acc.flux_up)
        μ₀ = @inbounds bcs_sw.cos_zenith[gcol]
        n_cloudy_gpts = 0
        if grp > n_gpt
            # More groups than g-points: this one has no work, and the
            # reduction reads its slice regardless.
            _zero_sw_partials!(acc, acc_clear, wcol, nlev, FT)
        else
            @inbounds for igpt in grp:ngroups:n_gpt
                cloudy = sw_2stream_gpt_col_both_skies!(
                    igpt,
                    gcol,
                    flux,
                    acc,
                    acc_clear,
                    nothing,
                    op,
                    bcs_sw,
                    src_sw,
                    as,
                    state_cache,
                    lookup_sw,
                    lookup_sw_cld,
                    lookup_sw_aero,
                    μ₀,
                    major_gpt2bnd[igpt],
                    n_gpt,
                    nlev,
                    mask,
                    grp,
                    wcol,
                    # G-points are contiguous within a band, so this is the
                    # last one of its own band: exactly one thread writes the
                    # 550 nm optical depth, and it is the g-point the ungrouped
                    # loop leaves standing.
                    igpt == n_gpt ||
                    major_gpt2bnd[igpt + 1] != major_gpt2bnd[igpt],
                )
                n_cloudy_gpts += cloudy ? 1 : 0
            end
            # The sun is down, so the loop above drew this column's masks but
            # never swept and never assigned a partial.
            μ₀ <= 0 && _zero_sw_partials!(acc, acc_clear, wcol, nlev, FT)
        end
        @inbounds cloudy_count[wcol] = Int32(n_cloudy_gpts)
    end
    return nothing
end

@inline function _zero_sw_partials!(acc, acc_clear, wcol, nlev, ::Type{FT}) where {FT}
    @inbounds for ilev in 1:nlev
        acc.flux_up[wcol, ilev] = FT(0)
        acc.flux_dn[wcol, ilev] = FT(0)
        acc.flux_dn_dir[wcol, ilev] = FT(0)
        acc_clear.flux_up[wcol, ilev] = FT(0)
        acc_clear.flux_dn[wcol, ilev] = FT(0)
        acc_clear.flux_dn_dir[wcol, ilev] = FT(0)
    end
    return nothing
end

# Sum the groups' partial fluxes into the broadband ones, finish the net flux,
# and total the cloud-cover tally. One thread per (column, level).
function _reduce_sw_groups_CUDA!(
    flux_sw::FluxSW,
    flux_sw_clear::FluxSW,
    acc_data,
    acc_clear_data,
    cloudy_count,
    cld_cover_sw,
    ngroups,
    n_gpt,
    ncol,
    nlev,
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= ncol * nlev
        gcol = (i - 1) % ncol + 1
        ilev = (i - 1) ÷ ncol + 1
        FT = eltype(flux_sw.flux_up)
        up, dn, dir = FT(0), FT(0), FT(0)
        cup, cdn, cdir = FT(0), FT(0), FT(0)
        @inbounds for g in 1:ngroups
            wcol = gcol + (g - 1) * ncol
            up += acc_data[wcol, ilev, 1]
            dn += acc_data[wcol, ilev, 2]
            dir += acc_data[wcol, ilev, 4]
            cup += acc_clear_data[wcol, ilev, 1]
            cdn += acc_clear_data[wcol, ilev, 2]
            cdir += acc_clear_data[wcol, ilev, 4]
        end
        @inbounds begin
            flux_sw.flux_up[gcol, ilev] = up
            flux_sw.flux_dn[gcol, ilev] = dn
            flux_sw.flux_dn_dir[gcol, ilev] = dir
            flux_sw.flux_net[gcol, ilev] = up - dn
            flux_sw_clear.flux_up[gcol, ilev] = cup
            flux_sw_clear.flux_dn[gcol, ilev] = cdn
            flux_sw_clear.flux_dn_dir[gcol, ilev] = cdir
            flux_sw_clear.flux_net[gcol, ilev] = cup - cdn
            if ilev == 1 && !isnothing(cld_cover_sw)
                n = 0
                for g in 1:ngroups
                    n += cloudy_count[gcol + (g - 1) * ncol]
                end
                cld_cover_sw[gcol] = FT(n) / n_gpt
            end
        end
    end
    return nothing
end
