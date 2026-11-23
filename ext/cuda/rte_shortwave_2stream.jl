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
    lookup_sw_precip::Union{LookUpPrecip, Nothing} = nothing,
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
        lookup_sw_precip,
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
    lookup_sw_precip::Union{LookUpPrecip, Nothing} = nothing,
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
                    lookup_sw_precip,
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
    lookup_sw_precip::Union{LookUpPrecip, Nothing} = nothing,
    gpt_groups = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    set_band_flux_to_zero!(band_flux)
    # The g-point loop splits only when no per-band fluxes are being kept: those
    # accumulate per g-point too, so they would need per-group partials of their
    # own and a wider reduction.
    if !isnothing(gpt_groups) && isnothing(band_flux)
        return _solve_both_skies_grouped!(
            gpt_groups,
            flux_sw,
            flux_sw_clear,
            bcs_sw,
            as,
            state_cache,
            lookup_sw,
            lookup_sw_cld,
            lookup_sw_aero,
            lookup_sw_precip,
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
        lookup_sw_precip,
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
# net-flux finalize and night-column zeroing are in both, and the grouped
# reduction below holds a third copy of the finalize and the zeroing. Change all
# three.
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
    lookup_sw_precip::Union{LookUpPrecip, Nothing} = nothing,
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
                    lookup_sw_precip,
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

#####
##### The g-point-grouped shortwave solve, as optics and then sweeps.
#####

# This round's gas and aerosol optics for every `(column, group)` pair whose sun
# is up, and for the one thread per dark column that owes it a 550 nm optical
# depth.
@inline _is_550nm_band(::Nothing, ibnd) = false
@inline _is_550nm_band(lkp::LookUpAerosolMerra, ibnd) = ibnd == lkp.iband_550nm

function _sw_group_optics_CUDA!(
    op::TwoStream,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
    ncol,
    ngroups,
    n_gpt,
    round,
    cos_zenith,
)
    wcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    (; major_gpt2bnd) = lookup_sw.band_data
    if wcol <= ncol * ngroups
        gcol = (wcol - 1) % ncol + 1
        grp = (wcol - 1) ÷ ncol + 1
        igpt = grp + (round - 1) * ngroups
        if igpt <= n_gpt
            @inbounds begin
                ibnd = major_gpt2bnd[igpt]
                # G-points are contiguous within a band, so this is the last
                # one of its own band: exactly one thread writes the 550 nm
                # optical depth, and it is the g-point the ungrouped loop
                # leaves standing.
                collect_aod = igpt == n_gpt || major_gpt2bnd[igpt + 1] != ibnd
                # Nothing else reads a dark column's optics: its thread in the
                # sweeps zeroes the partials instead of sweeping, and the cloud
                # increment is behind the same test. So of that column's 224
                # g-points only the one owing it an optical depth runs. `igpt`
                # is uniform across a block, so a block skips or runs whole.
                if cos_zenith[gcol] > 0 ||
                   (collect_aod && _is_550nm_band(lookup_sw_aero, ibnd))
                    compute_optical_props!(
                        op,
                        as,
                        state_cache,
                        gcol,
                        igpt,
                        lookup_sw,
                        nothing,
                        lookup_sw_aero,
                        nothing,
                        wcol,
                        collect_aod,
                    )
                end
            end
        end
    end
    return nothing
end

# The cloud mask and both sweeps for this round, over the optics the kernel
# above left in `op`. The partials of a night column and of a group with no
# g-point are zeroed on the first round rather than after the last: neither
# ever assigns one, so zeroing first is the same result and costs no second
# pass over the rounds.
function _sw_group_sweeps_CUDA!(
    flux::FluxSW,
    acc::NamedTuple,
    acc_clear::NamedTuple,
    src_sw::SourceSW2Str,
    bcs_sw::SwBCs,
    op::TwoStream,
    mask,
    cloudy_count,
    nlay,
    ncol,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_precip::Union{LookUpPrecip, Nothing},
    ngroups,
    n_gpt,
    round,
)
    wcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_sw.band_data
    if wcol <= ncol * ngroups
        gcol = (wcol - 1) % ncol + 1
        grp = (wcol - 1) ÷ ncol + 1
        FT = eltype(acc.flux_up)
        μ₀ = @inbounds bcs_sw.cos_zenith[gcol]
        igpt = grp + (round - 1) * ngroups
        if round == 1
            @inbounds cloudy_count[wcol] = Int32(0)
            (grp > n_gpt || μ₀ <= 0) &&
                _zero_partials!(acc, acc_clear, wcol, nlev)
        end
        if igpt <= n_gpt
            @inbounds begin
                cloudy = _build_cloud_mask!(
                    as.cloud_state,
                    Val(:mask_sw),
                    gcol,
                    igpt,
                    mask,
                    wcol,
                )
                sw_2stream_gpt_col_sweeps!(
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
                    lookup_sw,
                    lookup_sw_cld,
                    lookup_sw_precip,
                    μ₀,
                    major_gpt2bnd[igpt],
                    n_gpt,
                    nlev,
                    mask,
                    grp,
                    wcol,
                )
                cloudy && (cloudy_count[wcol] += Int32(1))
            end
        end
    end
    return nothing
end

# Sum the groups' partial fluxes into the broadband ones, finish the net flux,
# and total the cloud-cover tally. One thread per (column, level).
function _reduce_sw_groups_CUDA!(
    flux_sw::FluxSW,
    flux_sw_clear::FluxSW,
    acc,
    acc_clear,
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
            up += acc.flux_up[wcol, ilev]
            dn += acc.flux_dn[wcol, ilev]
            dir += acc.flux_dn_dir[wcol, ilev]
            cup += acc_clear.flux_up[wcol, ilev]
            cdn += acc_clear.flux_dn[wcol, ilev]
            cdir += acc_clear.flux_dn_dir[wcol, ilev]
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
