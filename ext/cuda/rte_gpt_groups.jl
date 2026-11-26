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
warp. The lanes of a column take its spectral bands in turn: they compute the
band's aerosol and cloud optics together, a layer each, and then take its
g-points `l`, `l + N`, ..., into their own slices `wcol` of the wide workspaces,
summing their fluxes over each level after every sweep. Summing in that order
instead of the serial one changes the fluxes at roundoff; the McICA masks are
unchanged, being stateless draws in `(key, column, g-point, layer)`.
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
        groups.band,
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

# The aerosol optics of band `ibnd` in this thread's column, filled by the
# column's lanes a layer each into planes 1 to 3 of `band`, and the optical depth
# and scattering optical depth of this lane's layers.
@inline function fill_band_aerosol!(
    band,
    as,
    state_cache,
    lkp_aero,
    ibnd,
    gcol,
    lane,
    lanes::Val{N},
    nlay,
    delta_scaling,
) where {N}
    (; aero_mask, aero_size, aero_mass) = as.aerosol_state
    rel_hum = AtmosphericStates.rel_hum_view(state_cache, as, gcol)
    FT = eltype(band)
    τ_ext = τ_sca = zero(FT)
    @inbounds for glay in lane:Int32(N):nlay
        if aero_mask[glay, gcol]
            (τ, ssa, g), (ext, sca) = Optics.aerosol_optics_2stream(
                lkp_aero,
                ibnd,
                view(aero_mass, :, :, gcol),
                view(aero_size, :, :, gcol),
                rel_hum[glay],
                glay,
                delta_scaling,
            )
            band[glay, 1, gcol], band[glay, 2, gcol], band[glay, 3, gcol] =
                τ, ssa, g
            τ_ext += ext
            τ_sca += sca
        end
    end
    return τ_ext, τ_sca
end

# The cloud optics of band `ibnd` in this thread's column, filled by the column's
# lanes a layer each into planes 4 to 6 of `band`, in every layer that any
# g-point's mask can select.
@inline function fill_band_clouds!(
    band,
    as,
    lkp_cld,
    lkp_precip,
    ibnd,
    gcol,
    lane,
    lanes::Val{N},
    nlay,
    delta_scaling,
) where {N}
    cs = as.cloud_state
    @inbounds for glay in lane:Int32(N):nlay
        if cs.cld_frac[glay, gcol] > 0
            band[glay, 4, gcol], band[glay, 5, gcol], band[glay, 6, gcol] =
                Optics.cloud_optics_2stream(
                    view(cs.cld_r_eff_liq, :, gcol),
                    view(cs.cld_r_eff_ice, :, gcol),
                    view(cs.cld_path_liq, :, gcol),
                    view(cs.cld_path_ice, :, gcol),
                    cs.ice_rgh,
                    lkp_cld,
                    ibnd,
                    glay;
                    lkp_precip,
                    cld_r_eff_sno = Optics._view_or_nothing(
                        cs.cld_r_eff_sno,
                        gcol,
                    ),
                    cld_path_sno = Optics._view_or_nothing(
                        cs.cld_path_sno,
                        gcol,
                    ),
                    delta_scaling,
                )
        end
    end
    return nothing
end

# Add the optics in planes `k + 1` to `k + 3` of `band` to this thread's optics,
# in the layers where `mask` is set.
@inline function add_band_optics!(op, band, k, mask, gcol, wcol, nlay)
    @inbounds for glay in 1i32:nlay
        if mask[glay]
            op.τ[wcol, glay], op.ssa[wcol, glay], op.g[wcol, glay] =
                Optics.increment_2stream(
                    op.τ[wcol, glay],
                    op.ssa[wcol, glay],
                    op.g[wcol, glay],
                    band[glay, k + 1i32, gcol],
                    band[glay, k + 2i32, gcol],
                    band[glay, k + 3i32, gcol],
                )
        end
    end
    return nothing
end

# The fused longwave solve (see `_solve_both_skies_grouped!`).
function lw_fused_CUDA!(
    fluxb::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    mask,
    band,
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
        (; major_gpt2bnd, bnd_lims_gpt) = lookup_lw.band_data
        n_gpt = length(major_gpt2bnd) % Int32
        (; cloud_state, aerosol_state) = as
        FT = eltype(flux_lw.flux_up)
        lane == 1i32 && _compute_aero_mask!(aerosol_state, gcol)
        n_cloudy = 0i32
        first = true
        @inbounds for ibnd in 1i32:(size(bnd_lims_gpt, 2) % Int32)
            CUDA.sync_warp(column_lanes(lanes))
            isnothing(lookup_lw_aero) || fill_band_aerosol!(
                band,
                as,
                state_cache,
                lookup_lw_aero,
                ibnd,
                gcol,
                lane,
                lanes,
                nlay,
                false,
            )
            isnothing(lookup_lw_cld) || fill_band_clouds!(
                band,
                as,
                lookup_lw_cld,
                lookup_lw_precip,
                ibnd,
                gcol,
                lane,
                lanes,
                nlay,
                false,
            )
            CUDA.sync_warp(column_lanes(lanes))
            gpt_last = bnd_lims_gpt[2, ibnd] % Int32
            for igpt_first in (bnd_lims_gpt[1, ibnd] % Int32):nlanes:gpt_last
                igpt = igpt_first + lane - 1i32
                active = igpt ≤ gpt_last
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
                        nothing,
                        nothing,
                        wcol,
                    )
                    isnothing(lookup_lw_aero) || add_band_optics!(
                        op,
                        band,
                        0i32,
                        view(aerosol_state.aero_mask, :, gcol),
                        gcol,
                        wcol,
                        nlay,
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
                    isnothing(lookup_lw_cld) || add_band_optics!(
                        op,
                        band,
                        3i32,
                        view(mask, :, wcol),
                        gcol,
                        wcol,
                        nlay,
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
                first = false
            end
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

# The fused shortwave solve (see `_solve_both_skies_grouped!`). A night column
# skips its optics and sweeps and is zeroed at the end, but still draws its cloud
# masks for the cloud cover, and its 550 nm band still collects the aerosol
# optical depth.
function sw_fused_CUDA!(
    fluxb::FluxSW,
    src_sw::SourceSW2Str,
    bcs_sw::SwBCs,
    op::TwoStream,
    mask,
    band,
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
        (; major_gpt2bnd, bnd_lims_gpt) = lookup_sw.band_data
        n_gpt = length(major_gpt2bnd) % Int32
        (; cloud_state, aerosol_state) = as
        FT = eltype(flux_sw.flux_up)
        day = @inbounds bcs_sw.cos_zenith[gcol] > 0
        lane == 1i32 && _compute_aero_mask!(aerosol_state, gcol)
        n_cloudy = 0i32
        first = true
        @inbounds for ibnd in 1i32:(size(bnd_lims_gpt, 2) % Int32)
            CUDA.sync_warp(column_lanes(lanes))
            if !isnothing(lookup_sw_aero) &&
               (day || ibnd == lookup_sw_aero.iband_550nm)
                τ_ext, τ_sca = fill_band_aerosol!(
                    band,
                    as,
                    state_cache,
                    lookup_sw_aero,
                    ibnd,
                    gcol,
                    lane,
                    lanes,
                    nlay,
                    true,
                )
                if ibnd == lookup_sw_aero.iband_550nm
                    τ_ext = lane_sum(τ_ext, lanes)
                    τ_sca = lane_sum(τ_sca, lanes)
                    if lane == 1i32
                        aerosol_state.aod_sw_ext[gcol] = τ_ext
                        aerosol_state.aod_sw_sca[gcol] = τ_sca
                    end
                end
            end
            day &&
                !isnothing(lookup_sw_cld) &&
                fill_band_clouds!(
                    band,
                    as,
                    lookup_sw_cld,
                    lookup_sw_precip,
                    ibnd,
                    gcol,
                    lane,
                    lanes,
                    nlay,
                    true,
                )
            CUDA.sync_warp(column_lanes(lanes))
            gpt_last = bnd_lims_gpt[2, ibnd] % Int32
            for igpt_first in (bnd_lims_gpt[1, ibnd] % Int32):nlanes:gpt_last
                igpt = igpt_first + lane - 1i32
                active = igpt ≤ gpt_last
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
                    if day
                        compute_optical_props!(
                            op,
                            as,
                            state_cache,
                            gcol,
                            igpt,
                            lookup_sw,
                            nothing,
                            nothing,
                            nothing,
                            wcol,
                        )
                        isnothing(lookup_sw_aero) || add_band_optics!(
                            op,
                            band,
                            0i32,
                            view(aerosol_state.aero_mask, :, gcol),
                            gcol,
                            wcol,
                            nlay,
                        )
                        rte_sw_2stream!(
                            op,
                            src_sw,
                            bcs_sw,
                            fluxb,
                            lookup_sw.solar_src_scaled[igpt],
                            igpt,
                            n_gpt,
                            ibnd,
                            nlev,
                            gcol,
                            wcol,
                        )
                    end
                end
                day && sum_lanes_into!(
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
                if active && day
                    isnothing(lookup_sw_cld) || add_band_optics!(
                        op,
                        band,
                        3i32,
                        view(mask, :, wcol),
                        gcol,
                        wcol,
                        nlay,
                    )
                    rte_sw_2stream!(
                        op,
                        src_sw,
                        bcs_sw,
                        fluxb,
                        lookup_sw.solar_src_scaled[igpt],
                        igpt,
                        n_gpt,
                        ibnd,
                        nlev,
                        gcol,
                        wcol,
                    )
                end
                day && sum_lanes_into!(
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
                first = false
            end
        end
        n_cloudy = lane_sum(n_cloudy, lanes)
        if lane == 1i32
            if day
                @inbounds compute_net_flux!(flux_sw, gcol, nlev)
                @inbounds compute_net_flux!(flux_sw_clear, gcol, nlev)
            else
                @inbounds set_flux_to_zero!(flux_sw, gcol, nlev)
                @inbounds set_flux_to_zero!(flux_sw_clear, gcol, nlev)
            end
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_sw)
                @inbounds cloud_state.cld_cover_sw[gcol] = FT(n_cloudy) / n_gpt
            end
        end
    end
    return nothing
end
