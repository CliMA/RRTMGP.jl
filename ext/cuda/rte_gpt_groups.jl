#####
##### The g-point-grouped both-skies solve: one fused kernel per band.
#####

# The fused kernels are register-bound, and capping their registers buys more
# warps than the spills cost. At this cap, a multiprocessor holds one block of
# 1024 threads, which was faster than more blocks of fewer threads.
const FUSED_BLOCK_THREADS = 1024
const FUSED_MAXREGS = 64

"""
    _solve_both_skies_grouped!(groups, flux, flux_clear, bcs, as, state_cache,
                               lookup, lookup_cld, lookup_aero, lookup_precip)

Run the both-skies solve of one band (`groups` is an `LWGptGroups` or an
`SWGptGroups`) with `N = n_gpt_groups(groups)` threads per column, which share a
warp. Lane `l` of a column takes the column's g-points `l`, `l + N`, ..., into
its own slice `wcol` of the wide workspaces, and the lanes sum their fluxes over
each level after every sweep. Summing in that order instead of the serial one
changes the fluxes at roundoff; the McICA masks are unchanged, being stateless
draws in `(key, column, g-point, layer)`.
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
    nthreads = ncol * n_gpt_groups(groups)
    tx, bx = _configure_threadblock(FUSED_BLOCK_THREADS, nthreads)
    args = (
        groups.fluxb,
        groups.src,
        bcs,
        groups.op,
        groups.mask,
        flux,
        flux_clear,
        Int32(nlay),
        Int32(ncol),
        Val(n_gpt_groups(groups)),
        as,
        state_cache,
        lookup,
        lookup_cld,
        lookup_aero,
        lookup_precip,
    )
    @cuda always_inline = true maxregs = FUSED_MAXREGS threads = tx blocks = bx fused_kernel(
        groups,
    )(
        int32_arrays(args)...,
    )
    return nothing
end

@inline fused_kernel(::LWGptGroups) = lw_fused_CUDA!
@inline fused_kernel(::SWGptGroups) = sw_fused_CUDA!

# The lanes of this thread's column: `N` consecutive lanes of its warp, `N` a
# power of two (ClimaCore's `subwarp_lanes`).
@inline column_lanes(::Val{N}) where {N} =
    N == 32 ? CUDA.FULL_MASK :
    ((one(UInt32) << N) - one(UInt32)) <<
    ((CUDA.laneid() - 1i32) % UInt32 & ~UInt32(N - 1))

# The sum of `x` over the lanes of a column, in each of them (ClimaCore's
# `shuffle_reduce`).
@inline function lane_sum(x, lanes::Val{N}) where {N}
    for offset in ntuple(n -> one(UInt32) << (n - 1), Val(trailing_zeros(N)))
        x += CUDA.shfl_xor_sync(column_lanes(lanes), x, offset, N % UInt32)
    end
    return x
end

@inline flux_planes(flux::FluxLW) = (flux.flux_up, flux.flux_dn)
@inline flux_planes(flux::FluxSW) =
    (flux.flux_up, flux.flux_dn, flux.flux_dn_dir)

# Add the fluxes that a column's lanes left in `fluxb` to the column's fluxes in
# `acc`, level by level, or overwrite them on the first round.
@inline function sum_lanes_into!(
    acc,
    fluxb,
    active,
    gcol,
    wcol,
    nlev,
    first,
    lanes,
    lane,
)
    for ilev in 1i32:nlev
        foreach(flux_planes(acc), flux_planes(fluxb)) do a, b
            v = active ? (@inbounds b[wcol, ilev]) : zero(eltype(b))
            v = lane_sum(v, lanes)
            lane == 1i32 &&
                (@inbounds a[gcol, ilev] = first ? v : a[gcol, ilev] + v)
        end
    end
    return nothing
end

# The fused longwave solve: `N` consecutive threads per column, each taking every
# `N`-th g-point of it into its own slice (`wcol`) of the wide workspaces, and
# the column's fluxes summed over its lanes after every sweep.
function lw_fused_CUDA!(
    fluxb::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    mask,
    flux_lw::FluxLW,
    flux_lw_clear::FluxLW,
    nlay,
    ncol,
    lanes::Val{N},
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
    lookup_lw_precip::Union{LookUpPrecip, Nothing},
) where {N}
    nlanes = Int32(N)
    wcol = threadIdx().x + (blockIdx().x - 1i32) * blockDim().x
    gcol = ((wcol - 1i32) >> trailing_zeros(nlanes)) + 1i32
    lane = ((wcol - 1i32) & (nlanes - 1i32)) + 1i32
    if gcol ≤ ncol
        nlev = nlay + 1i32
        (; major_gpt2bnd) = lookup_lw.band_data
        n_gpt = length(major_gpt2bnd) % Int32
        (; cloud_state, aerosol_state) = as
        FT = eltype(flux_lw.flux_up)
        lane == 1i32 && _compute_aero_mask!(aerosol_state, gcol)
        CUDA.sync_warp(column_lanes(lanes))
        n_cloudy = 0i32
        @inbounds for round in 1i32:cld(n_gpt, nlanes)
            igpt = (round - 1i32) * nlanes + lane
            active = igpt ≤ n_gpt
            first = round == 1i32
            ibnd = active ? major_gpt2bnd[igpt] : 1
            if active
                cloudy = _build_cloud_mask!(
                    cloud_state,
                    Val(:mask_lw),
                    gcol,
                    igpt,
                    mask,
                    wcol,
                )
                n_cloudy += Int32(cloudy)
                compute_optical_props!(
                    op,
                    as,
                    state_cache,
                    src_lw,
                    gcol,
                    igpt,
                    lookup_lw,
                    nothing,
                    lookup_lw_aero,
                    nothing,
                    wcol,
                )
                rte_lw_2stream!(
                    op,
                    fluxb,
                    src_lw,
                    bcs_lw,
                    gcol,
                    igpt,
                    ibnd,
                    nlev,
                    ncol,
                    wcol,
                )
            end
            sum_lanes_into!(
                flux_lw_clear,
                fluxb,
                active,
                gcol,
                wcol,
                nlev,
                first,
                lanes,
                lane,
            )
            if active
                isnothing(lookup_lw_cld) || add_cloud_optics_lw!(
                    op,
                    as,
                    gcol,
                    lookup_lw_cld,
                    ibnd,
                    lookup_lw_precip,
                    wcol,
                    mask,
                )
                rte_lw_2stream!(
                    op,
                    fluxb,
                    src_lw,
                    bcs_lw,
                    gcol,
                    igpt,
                    ibnd,
                    nlev,
                    ncol,
                    wcol,
                )
            end
            sum_lanes_into!(
                flux_lw,
                fluxb,
                active,
                gcol,
                wcol,
                nlev,
                first,
                lanes,
                lane,
            )
        end
        n_cloudy = lane_sum(n_cloudy, lanes)
        if lane == 1i32
            @inbounds compute_net_flux!(flux_lw, gcol, nlev)
            @inbounds compute_net_flux!(flux_lw_clear, gcol, nlev)
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_lw)
                @inbounds cloud_state.cld_cover_lw[gcol] = FT(n_cloudy) / n_gpt
            end
        end
    end
    return nothing
end

# The fused shortwave solve (see `lw_fused_CUDA!`). A night column skips its
# sweeps and is zeroed at the end, but still draws its cloud masks for the cloud
# cover, and the last g-point of the 550 nm band still collects the aerosol
# optical depth.
function sw_fused_CUDA!(
    fluxb::FluxSW,
    src_sw::SourceSW2Str,
    bcs_sw::SwBCs,
    op::TwoStream,
    mask,
    flux_sw::FluxSW,
    flux_sw_clear::FluxSW,
    nlay,
    ncol,
    lanes::Val{N},
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
    lookup_sw_precip::Union{LookUpPrecip, Nothing},
) where {N}
    nlanes = Int32(N)
    wcol = threadIdx().x + (blockIdx().x - 1i32) * blockDim().x
    gcol = ((wcol - 1i32) >> trailing_zeros(nlanes)) + 1i32
    lane = ((wcol - 1i32) & (nlanes - 1i32)) + 1i32
    if gcol ≤ ncol
        nlev = nlay + 1i32
        (; major_gpt2bnd) = lookup_sw.band_data
        n_gpt = length(major_gpt2bnd) % Int32
        (; cloud_state, aerosol_state) = as
        FT = eltype(flux_sw.flux_up)
        μ₀ = @inbounds bcs_sw.cos_zenith[gcol]
        lane == 1i32 && _compute_aero_mask!(aerosol_state, gcol)
        CUDA.sync_warp(column_lanes(lanes))
        n_cloudy = 0i32
        @inbounds for round in 1i32:cld(n_gpt, nlanes)
            igpt = (round - 1i32) * nlanes + lane
            active = igpt ≤ n_gpt
            first = round == 1i32
            ibnd = active ? major_gpt2bnd[igpt] : 1
            if active
                cloudy = _build_cloud_mask!(
                    cloud_state,
                    Val(:mask_sw),
                    gcol,
                    igpt,
                    mask,
                    wcol,
                )
                n_cloudy += Int32(cloudy)
                collect_aod =
                    igpt == n_gpt || major_gpt2bnd[igpt + 1i32] != ibnd
                if μ₀ > 0 ||
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
                if μ₀ > 0
                    solar_frac = lookup_sw.solar_src_scaled[igpt]
                    rte_sw_2stream!(
                        op,
                        src_sw,
                        bcs_sw,
                        fluxb,
                        solar_frac,
                        igpt,
                        n_gpt,
                        ibnd,
                        nlev,
                        gcol,
                        wcol,
                    )
                end
            end
            μ₀ > 0 && sum_lanes_into!(
                flux_sw_clear,
                fluxb,
                active,
                gcol,
                wcol,
                nlev,
                first,
                lanes,
                lane,
            )
            if active && μ₀ > 0
                isnothing(lookup_sw_cld) || add_cloud_optics_sw!(
                    op,
                    as,
                    gcol,
                    lookup_sw_cld,
                    ibnd,
                    lookup_sw_precip,
                    wcol,
                    mask,
                )
                solar_frac = lookup_sw.solar_src_scaled[igpt]
                rte_sw_2stream!(
                    op,
                    src_sw,
                    bcs_sw,
                    fluxb,
                    solar_frac,
                    igpt,
                    n_gpt,
                    ibnd,
                    nlev,
                    gcol,
                    wcol,
                )
            end
            μ₀ > 0 && sum_lanes_into!(
                flux_sw,
                fluxb,
                active,
                gcol,
                wcol,
                nlev,
                first,
                lanes,
                lane,
            )
        end
        n_cloudy = lane_sum(n_cloudy, lanes)
        if lane == 1i32
            if μ₀ ≤ 0
                @inbounds set_flux_to_zero!(flux_sw, gcol, nlev)
                @inbounds set_flux_to_zero!(flux_sw_clear, gcol, nlev)
            else
                @inbounds compute_net_flux!(flux_sw, gcol, nlev)
                @inbounds compute_net_flux!(flux_sw_clear, gcol, nlev)
            end
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_sw)
                @inbounds cloud_state.cld_cover_sw[gcol] = FT(n_cloudy) / n_gpt
            end
        end
    end
    return nothing
end


@inline _is_550nm_band(::Nothing, ibnd) = false
@inline _is_550nm_band(lkp::LookUpAerosolMerra, ibnd) = ibnd == lkp.iband_550nm
