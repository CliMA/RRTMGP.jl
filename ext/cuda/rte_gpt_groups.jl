#####
##### The g-point-grouped both-skies solve: one fused kernel per band.
#####

# The register cap and threads per block of each fused kernel. The longwave
# kernel keeps more values live through its sweep up than the shortwave one, and
# runs fastest with enough registers to spill little, in small blocks that fill a
# multiprocessor at that cap; the shortwave kernel runs fastest at a lower cap,
# in blocks of 256 threads, four to a multiprocessor (larger and smaller blocks
# at that cap are slower). (Measured on an A100 with the flagship AMIP
# configuration.)
@inline fused_maxregs(::LWGptGroups) = 96
@inline fused_threads(::LWGptGroups) = 128
@inline fused_maxregs(::SWGptGroups) = 64
@inline fused_threads(::SWGptGroups) = 256

"""
    _solve_both_skies_grouped!(groups, flux, flux_clear, bcs, as, state_cache,
                               lookup, lookup_cld, lookup_aero, lookup_precip)

Run the both-skies solve (`groups` is an `LWGptGroups` or an `SWGptGroups`) with
`N = n_gpt_groups(groups)` threads per column, which share a warp. The lanes of
a column take its spectral bands in turn. For each band, they first compute
together, a layer each, what the band's g-points share: its aerosol and cloud
optics and the parts of its gas optics that do not depend on the g-point
(`fill_band_aerosol!`, `fill_band_clouds!` and `fill_band_gas!`). Each lane then
solves its g-points `l`, `l + N`, ..., computing each layer's optics during its
sweeps of both skies (`lw_stream!` and `sw_stream!`), and the lanes sum their
fluxes in registers during the sweep down (`sweep_down_sums!`). Summing in that
order instead of the serial one changes the fluxes at roundoff; the McICA masks
are unchanged, being stateless draws in `(key, column, g-point, layer)`.
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
    tx, bx = _configure_threadblock(fused_threads(groups), nthreads)
    args = (
        groups.fluxb,
        groups.src,
        bcs,
        groups.op,
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
    kargs = int32_arrays(args)
    kernel = @cuda launch = false always_inline = true maxregs =
        fused_maxregs(groups) fused_CUDA!(kargs...)
    # the lanes' cloud masks, a bit for each layer (see `lane_masks`)
    shmem = tx * cld(nlay, 32) * sizeof(UInt32)
    if shmem > 48 * 1024
        CUDA.attributes(kernel.fun)[CUDA.FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES] =
            shmem
    end
    kernel(kargs...; threads = tx, blocks = bx, shmem)
    return nothing
end


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

# The parts of the aerosol optics of this thread's column that do not depend on
# the band (see `Optics.aerosol_searches`), filled by the column's lanes a layer
# each into the planes from `BAND_AEROSOL_SEARCHES` of `band`, which the bands'
# fills leave alone.
@inline function fill_aerosol_searches!(
    band,
    as,
    state_cache,
    lkp_aero,
    gcol,
    lane,
    ::Val{N},
    nlay,
) where {N}
    (; aero_mask, aero_size, aero_mass) = as.aerosol_state
    rel_hum = AtmosphericStates.rel_hum_view(state_cache, as, gcol)
    FT = eltype(band)
    k = BAND_AEROSOL_SEARCHES
    @inbounds for glay in lane:Int32(N):nlay
        if aero_mask[glay, gcol]
            loc, factor, bins = Optics.aerosol_searches(
                lkp_aero,
                view(aero_mass, :, :, gcol),
                view(aero_size, :, :, gcol),
                rel_hum[glay],
                glay,
            )
            band[glay, k, gcol] = pack_bits(FT, UInt32(loc))
            band[glay, k + 1i32, gcol] = factor
            band[glay, k + 2i32, gcol] = pack_bits(FT, bins)
        end
    end
    return nothing
end

# The aerosol optics of band `ibnd` in this thread's column, filled by the
# column's lanes a layer each into planes 1 to 3 of `band` from the searches
# that `fill_aerosol_searches!` left, and the layers' optical depths and
# scattering optical depths into planes 7 and 8.
@inline function fill_band_aerosol!(
    band,
    as,
    lkp_aero,
    ibnd,
    gcol,
    lane,
    ::Val{N},
    nlay,
    delta_scaling,
) where {N}
    (; aero_mask, aero_mass) = as.aerosol_state
    k = BAND_AEROSOL_SEARCHES
    @inbounds for glay in lane:Int32(N):nlay
        if aero_mask[glay, gcol]
            searches = (
                Int(unpack_bits(band[glay, k, gcol])),
                band[glay, k + 1i32, gcol],
                unpack_bits(band[glay, k + 2i32, gcol]),
            )
            (τ, ssa, g), (ext, sca) = Optics.aerosol_optics_2stream(
                lkp_aero,
                ibnd,
                view(aero_mass, :, :, gcol),
                glay,
                searches,
                delta_scaling,
            )
            band[glay, 1, gcol], band[glay, 2, gcol], band[glay, 3, gcol] =
                τ, ssa, g
            band[glay, 7, gcol], band[glay, 8, gcol] = ext, sca
        end
    end
    return nothing
end

# The aerosol optical depth and scattering optical depth of this thread's column,
# summed in layer order from planes 7 and 8 of `band`, as the serial solve sums
# them.
@inline function aerosol_optical_depths(band, aero_mask, gcol, nlay)
    FT = eltype(band)
    τ_ext = τ_sca = zero(FT)
    @inbounds for glay in 1i32:nlay
        if aero_mask[glay, gcol]
            τ_ext += band[glay, 7, gcol]
            τ_sca += band[glay, 8, gcol]
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

# The interpolation indices of a layer, packed into the bits of one float of the
# band workspace: `jtemp` (4 bits), `jpress` (6), `jη1` and `jη2` (4 each), the
# troposphere flag (1), and the first minor contributor of the band in that part
# of the atmosphere (7) and their number (4).
@inline pack_bits(::Type{Float32}, u::UInt32) = reinterpret(Float32, u)
@inline pack_bits(::Type{Float64}, u::UInt32) = reinterpret(Float64, UInt64(u))
@inline unpack_bits(x::Float32) = reinterpret(UInt32, x)
@inline unpack_bits(x::Float64) = reinterpret(UInt64, x) % UInt32

# The parts of the gas optics of band `ibnd` in this thread's column that its
# g-points share, filled by the column's lanes a layer (and level) each into the
# band workspace (see `BAND_GAS_INDICES` and the planes after it).
@inline function fill_band_gas!(
    band,
    as,
    state_cache,
    lkp,
    ibnd,
    gcol,
    lane,
    ::Val{N},
    nlay,
) where {N}
    (; vmr) = as
    FT = eltype(band)
    col_dry_col = AtmosphericStates.col_dry_view(state_cache, as, gcol)
    p_lay_col = AtmosphericStates.p_lay_view(state_cache, as, gcol)
    t_lay_col = AtmosphericStates.t_lay_view(state_cache, as, gcol)
    @inbounds for glay in lane:Int32(N):nlay
        col_dry, p_lay, t_lay =
            col_dry_col[glay], p_lay_col[glay], t_lay_col[glay]
        tropo,
        vmr_h2o,
        (jtemp, ftemp),
        (jpress, fpress),
        (jη1, jη2, fη1, fη2),
        (col_mix1, col_mix2) =
            Optics.gas_layer_interp(lkp, vmr, ibnd, p_lay, t_lay, glay, gcol)
        lkp_minor = tropo == 1 ? lkp.minor_lower : lkp.minor_upper
        st, n = LookUpTables.minor_bounds(lkp_minor, ibnd)
        bits =
            UInt32(jtemp) | UInt32(jpress) << 4 | UInt32(jη1) << 10 |
            UInt32(jη2) << 14 | UInt32(tropo - 1) << 18 | UInt32(st) << 19 |
            UInt32(n) << 26
        band[glay, BAND_GAS_INDICES, gcol] = pack_bits(FT, bits)
        k = BAND_GAS_FRACTIONS
        band[glay, k, gcol] = ftemp
        band[glay, k + 1i32, gcol] = fpress
        band[glay, k + 2i32, gcol] = fη1
        band[glay, k + 3i32, gcol] = fη2
        band[glay, k + 4i32, gcol] = col_mix1
        band[glay, k + 5i32, gcol] = col_mix2
        lkp isa LookUpSW && (band[glay, BAND_PLANCK, gcol] = vmr_h2o)
        for i in 0i32:(Int32(n) - 1i32)
            band[glay, BAND_MINOR + i, gcol] = Optics.minor_scaling(
                lkp_minor,
                vmr,
                vmr_h2o,
                col_dry,
                p_lay,
                t_lay,
                st + i,
                glay,
                gcol,
            )
        end
    end
    if lkp isa LookUpLW
        (; t_planck, tot_planck) = lkp.planck
        totplnk = view(tot_planck, :, ibnd)
        t_lev_col = AtmosphericStates.t_lev_view(state_cache, as, gcol)
        @inbounds for glev in lane:Int32(N):(nlay + 1i32)
            band[glev, BAND_PLANCK, gcol] =
                Optics.interp1d_equispaced(t_lev_col[glev], t_planck, totplnk)
        end
        if lane == 1i32
            @inbounds band[1, BAND_PLANCK_SFC, gcol] =
                Optics.interp1d_equispaced(as.t_sfc[gcol], t_planck, totplnk)
        end
    end
    return nothing
end

# The band data of layer `glay` that the gas optics of each of its g-points reads
# (besides the minor absorber scalings), which a sweep can load ahead of use:
# packed, so that it holds few registers until `gas_layer` unpacks it.
@inline function band_gas_data(band, glay, gcol, col_dry)
    k = BAND_GAS_FRACTIONS
    @inbounds return (
        unpack_bits(band[glay, BAND_GAS_INDICES, gcol]),
        band[glay, k, gcol],
        band[glay, k + 1i32, gcol],
        band[glay, k + 2i32, gcol],
        band[glay, k + 3i32, gcol],
        band[glay, k + 4i32, gcol],
        band[glay, k + 5i32, gcol],
        col_dry,
    )
end

# The interpolation of a layer (see `Optics.gas_layer_interp`, with the water
# vapor mixing ratio `vmr_h2o`) and its band's minor contributors, unpacked from
# the layer's band data.
@inline function gas_layer(
    (bits, ftemp, fpress, fη1, fη2, col_mix1, col_mix2, _),
    vmr_h2o,
)
    layer = (
        Int32((bits >> 18) & 0x1) + 1i32,
        vmr_h2o,
        (Int(bits & 0xf), ftemp),
        (Int((bits >> 4) & 0x3f), fpress),
        (Int((bits >> 10) & 0xf), Int((bits >> 14) & 0xf), fη1, fη2),
        (col_mix1, col_mix2),
    )
    return layer, (Int32((bits >> 19) & 0x7f), Int32((bits >> 26) & 0xf))
end

# The major and minor optical depths of the `loc`-th g-point `igpt` of its band
# in layer `glay` (see `Optics.gas_τ`), with the minor scalings that
# `fill_band_gas!` left in the band workspace, and the layer's interpolation.
@inline function band_gas_τ(lkp, band, data, igpt, loc, glay, gcol, vmr_h2o)
    layer, (st, n) = gas_layer(data, vmr_h2o)
    τ = Optics.gas_τ(
        lkp,
        layer,
        last(data),
        igpt,
        loc,
        (st, n),
        imnr -> @inbounds(band[glay, BAND_MINOR + (imnr - st), gcol]),
        Val(false),
    )
    return τ, layer
end

# The longwave gas optical depth and Planck fraction of g-point `igpt` in layer
# `glay` (see `Optics.compute_gas_optics(::LookUpLW, ...)`).
@inline lw_band_gas(lkp, band, igpt, loc, glay, gcol, col_dry) = lw_band_gas(
    lkp,
    band,
    band_gas_data(band, glay, gcol, col_dry),
    igpt,
    loc,
    glay,
    gcol,
)
@inline function lw_band_gas(lkp, band, data::Tuple, igpt, loc, glay, gcol)
    (τ_major, τ_minor), layer =
        band_gas_τ(lkp, band, data, igpt, loc, glay, gcol, zero(last(data)))
    pfrac = Optics.planck_fraction(lkp, layer, igpt)
    return max(τ_major + τ_minor, zero(τ_major)), pfrac
end

# The shortwave gas optical depth and single-scattering albedo of g-point `igpt`
# in layer `glay` (see `Optics.compute_gas_optics(::LookUpSW, ...)`).
@inline function sw_band_gas(lkp, band, igpt, loc, glay, gcol, col_dry)
    @inbounds vmr_h2o = band[glay, BAND_PLANCK, gcol]
    data = band_gas_data(band, glay, gcol, col_dry)
    (τ_major, τ_minor), layer =
        band_gas_τ(lkp, band, data, igpt, loc, glay, gcol, vmr_h2o)
    τ_ray = Optics.rayleigh_τ(lkp, layer, col_dry, igpt)
    return Optics.gas_τ_ssa(τ_major, τ_minor, τ_ray)
end

# Whether layer `ilay` has aerosols and whether the g-point's cloud mask selects
# it, with their optics from the band workspace (zeros where absent), which the
# longwave sweep loads a layer ahead of use.
@inline function layer_extras(
    band,
    as,
    mask,
    lkp_aero,
    lkp_cld,
    ilay,
    gcol,
    wcol,
    bounds,
)
    FT = eltype(band)
    aero =
        !isnothing(lkp_aero) && @inbounds as.aerosol_state.aero_mask[ilay, gcol]
    cld = !isnothing(lkp_cld) && masked(mask, ilay, wcol, bounds)
    ka, kc = BAND_AEROSOL, BAND_CLOUD
    @inbounds aero_optics =
        aero ?
        (
            band[ilay, ka + 1i32, gcol],
            band[ilay, ka + 2i32, gcol],
            band[ilay, ka + 3i32, gcol],
        ) : (FT(0), FT(0), FT(0))
    @inbounds cld_optics =
        cld ?
        (
            band[ilay, kc + 1i32, gcol],
            band[ilay, kc + 2i32, gcol],
            band[ilay, kc + 3i32, gcol],
        ) : (FT(0), FT(0), FT(0))
    return aero, cld, aero_optics, cld_optics
end

# The longwave two-stream solve of g-point `igpt` for both skies, in one sweep up
# that computes each layer's optics as it reaches it, and one sweep down. The
# clear sky's coefficients serve the cloudy sky in the layers that the g-point's
# cloud mask leaves clear. Each sky keeps the two coefficients of each layer
# that the sweep down needs (`a` and `b` below) and the albedo and source of
# each level, which the sweep down overwrites with the upward and downward
# fluxes: the clear sky in `fluxb.flux_up`, `fluxb.flux_dn`, `src_lw.albedo` and
# `src_lw.src`, and the cloudy sky in `op.τ`, `op.ssa`, `fluxb.flux_net` and
# `src_lw.lev_source`. Below the lowest layer that the cloud mask selects, the
# cloudy sky is the clear sky, so the cloudy sky's values are stored only from
# that layer up, and the sweep returns the layer (`nlay + 1` if the mask selects
# none) for the sweep down, which reads the clear sky's values for both below it.
# The sweep loads each layer's band data a layer ahead of use (which, unlike in
# the shortwave, measurably speeds it up). The arithmetic is that of the serial
# solve's `compute_optical_props!`, aerosol and cloud increments, and
# `rte_lw_2stream!`.
@inline function lw_stream!(
    op,
    fluxb,
    src_lw,
    bcs_lw,
    mask,
    band,
    as,
    state_cache,
    lkp,
    lkp_cld,
    lkp_aero,
    gcol,
    wcol,
    igpt,
    ibnd,
    loc,
    bounds,
    nlay,
)
    FT = eltype(band)
    nlev = nlay + 1i32
    (; inc_flux, sfc_emis) = bcs_lw
    col_dry_col = AtmosphericStates.col_dry_view(state_cache, as, gcol)
    a_clr, b_clr, up_clr, dn_clr =
        map(streamed, (fluxb.flux_up, fluxb.flux_dn, src_lw.albedo, src_lw.src))
    a_all, b_all, up_all, dn_all =
        map(streamed, (op.τ, op.ssa, fluxb.flux_net, src_lw.lev_source))
    @inbounds begin
        flux_dn_top = isnothing(inc_flux) ? FT(0) : inc_flux[gcol, igpt]
        emis = sfc_emis[ibnd, gcol]
        τ_gas, pfrac =
            lw_band_gas(lkp, band, igpt, loc, 1i32, gcol, col_dry_col[1])
        # loaded a layer ahead of use: the band data of the gas optics of the
        # next layer, and the Planck function, aerosol optics and cloud mask of
        # this layer and the next
        n2 = min(2i32, nlay)
        data_next = band_gas_data(band, n2, gcol, col_dry_col[n2])
        extras = layer_extras(
            band,
            as,
            mask,
            lkp_aero,
            lkp_cld,
            1i32,
            gcol,
            wcol,
            bounds,
        )
        planck_top_next = band[2, BAND_PLANCK, gcol]
        albedo_clr = FT(1) - emis
        src_clr = FT(π) * emis * (band[1, BAND_PLANCK_SFC, gcol] * pfrac)
        albedo_all, src_all = albedo_clr, src_clr
        up_clr[wcol, 1], dn_clr[wcol, 1] = albedo_clr, src_clr
        lev_src_bot = band[1, BAND_PLANCK, gcol] * pfrac
    end
    lo = nlev
    @inbounds for ilay in 1i32:nlay
        data_cur, planck_top = data_next, planck_top_next
        if ilay + 2i32 ≤ nlay
            data_next =
                band_gas_data(band, ilay + 2i32, gcol, col_dry_col[ilay + 2i32])
        end
        if ilay < nlay
            planck_top_next = band[ilay + 2i32, BAND_PLANCK, gcol]
            extras_next = layer_extras(
                band,
                as,
                mask,
                lkp_aero,
                lkp_cld,
                ilay + 1i32,
                gcol,
                wcol,
                bounds,
            )
        else
            extras_next = extras
        end
        lev_src_inc = planck_top * pfrac
        if ilay < nlay
            τ_next, pfrac_next =
                lw_band_gas(lkp, band, data_cur, igpt, loc, ilay + 1i32, gcol)
            lev_src_top = sqrt(lev_src_inc * (planck_top * pfrac_next))
        else
            τ_next, pfrac_next = τ_gas, pfrac
            lev_src_top = lev_src_inc
        end
        τ, ssa, g = τ_gas, FT(0), FT(0)
        aero, cld, aero_optics, cld_optics = extras
        if aero
            τ, ssa, g = Optics.increment_2stream(τ, ssa, g, aero_optics...)
        end
        Rdif, Tdif, src_up, src_dn =
            lw_2stream_coeffs(τ, ssa, g, lev_src_bot, lev_src_top)
        a_clr[wcol, ilay], b_clr[wcol, ilay], albedo_clr, src_clr =
            adding_layer(Rdif, Tdif, src_up, src_dn, albedo_clr, src_clr)
        up_clr[wcol, ilay + 1i32], dn_clr[wcol, ilay + 1i32] =
            albedo_clr, src_clr
        # Below the lowest layer that the cloud mask selects, the cloudy sky is
        # the clear sky, whose planes the sweep down reads for both.
        lo = cld ? min(lo, ilay) : lo
        if lo ≤ ilay
            if cld
                τ, ssa, g = Optics.increment_2stream(τ, ssa, g, cld_optics...)
                Rdif, Tdif, src_up, src_dn =
                    lw_2stream_coeffs(τ, ssa, g, lev_src_bot, lev_src_top)
            end
            a_all[wcol, ilay], b_all[wcol, ilay], albedo_all, src_all =
                adding_layer(Rdif, Tdif, src_up, src_dn, albedo_all, src_all)
            up_all[wcol, ilay + 1i32], dn_all[wcol, ilay + 1i32] =
                albedo_all, src_all
        else
            albedo_all, src_all = albedo_clr, src_clr
        end
        lev_src_bot = lev_src_top
        τ_gas, pfrac = τ_next, pfrac_next
        extras = extras_next
    end
    # Eq 12 at the top of the domain; `sweep_down_sums!` sweeps down
    @inbounds begin
        up_clr[wcol, nlev] = flux_dn_top * albedo_clr + src_clr
        dn_clr[wcol, nlev] = flux_dn_top
        lo < nlev && (up_all[wcol, nlev] = flux_dn_top * albedo_all + src_all)
    end
    return lo
end

# The shortwave two-stream solve of g-point `igpt` for both skies, in a sweep
# down that computes each layer's optics as it reaches it, with the direct beam
# of each sky, a sweep up that computes the layers' coefficients, reusing the
# clear sky's in layers that the g-point's cloud mask leaves clear, and a second
# sweep down for the fluxes. The sweep down keeps the clear sky's optics in
# `op`, which the sweep up overwrites with the two coefficients of each layer
# that the second sweep down needs (`a` and `b` below), those of the clear sky in
# `op.τ` and `op.ssa` and those of the cloudy sky in `op.g` and the third plane
# of `src_sw.leveldata`. Each sky keeps the direct beam (in `fluxb.flux_dn_dir`
# and `fluxb.flux_net`), and the albedo and source of each level, which the
# second sweep down overwrites with the upward and downward fluxes (the clear
# sky in `src_sw.albedo` and `src_sw.src`, the cloudy sky in `fluxb.flux_up` and
# `fluxb.flux_dn`). The cloudy sky's values are stored only where they can differ
# from the clear sky's: nowhere for a g-point that the cloud mask leaves clear,
# and otherwise everywhere but the albedos and coefficients `a` below the lowest
# layer that the mask selects and the direct beam above the highest. The sweep
# returns those two layers (`nlay + 1` and 0 if the mask selects none) for the
# second sweep down, which reads the clear sky's values for both there. The
# arithmetic is that of the serial solve's `compute_optical_props!`, aerosol and
# cloud increments, and `rte_sw_2stream!`.
@inline function sw_stream!(
    op,
    fluxb,
    src_sw,
    bcs_sw,
    mask,
    band,
    as,
    state_cache,
    lkp,
    lkp_cld,
    lkp_aero,
    gcol,
    wcol,
    igpt,
    ibnd,
    loc,
    bounds,
    nlay,
)
    FT = eltype(band)
    nlev = nlay + 1i32
    col_dry_col = AtmosphericStates.col_dry_view(state_cache, as, gcol)
    a_clr, b_clr, a_all, b_all =
        map(streamed, (op.τ, op.ssa, op.g, src_sw.leveldata))
    up_clr, dn_clr, dir_clr =
        map(streamed, (src_sw.albedo, src_sw.src, fluxb.flux_dn_dir))
    up_all, dn_all, dir_all =
        map(streamed, (fluxb.flux_up, fluxb.flux_dn, fluxb.flux_net))
    @inbounds begin
        toa_flux = bcs_sw.toa_flux[gcol]
        sfc_alb_direct = bcs_sw.sfc_alb_direct[ibnd, gcol]
        sfc_alb_diffuse = bcs_sw.sfc_alb_diffuse[ibnd, gcol]
        μ₀ = bcs_sw.cos_zenith[gcol]
        solar_frac = lkp.solar_src_scaled[igpt]
    end
    flux_dn_dir_top = toa_flux * solar_frac * μ₀
    inv_μ₀ = inv(max(μ₀, Numerics.μ₀_min(FT)))
    @inbounds dir_clr[wcol, nlev] = flux_dn_dir_top
    τ_cum_clr = τ_cum_all = FT(0)
    lo, hi = nlev, 0i32
    ilay = nlay
    @inbounds while ilay ≥ 1i32
        τ, ssa =
            sw_band_gas(lkp, band, igpt, loc, ilay, gcol, col_dry_col[ilay])
        g = FT(0)
        if !isnothing(lkp_aero) && as.aerosol_state.aero_mask[ilay, gcol]
            k = BAND_AEROSOL
            τ, ssa, g = Optics.increment_2stream(
                τ,
                ssa,
                g,
                band[ilay, k + 1i32, gcol],
                band[ilay, k + 2i32, gcol],
                band[ilay, k + 3i32, gcol],
            )
        end
        a_clr[wcol, ilay], b_clr[wcol, ilay], a_all[wcol, ilay] = τ, ssa, g
        τ_cum_clr += τ
        dir_clr[wcol, ilay] = flux_dn_dir_top * exp(-τ_cum_clr * inv_μ₀)
        if !isnothing(lkp_cld) && masked(mask, ilay, wcol, bounds)
            τ += band[ilay, BAND_CLOUD + 1i32, gcol]
            lo, hi = ilay, max(hi, ilay)
        end
        τ_cum_all += τ
        # Above the highest layer that the cloud mask selects, the cloudy sky's
        # direct beam is the clear sky's, which the sweeps read for both.
        if ilay ≤ hi
            dir_all[wcol, ilay] = flux_dn_dir_top * exp(-τ_cum_all * inv_μ₀)
        end
        ilay -= 1i32
    end
    # A g-point that the cloud mask leaves clear has the clear sky's fluxes in
    # the cloudy sky. Otherwise, the cloudy sky's albedos and coefficients `a`
    # are the clear sky's below the lowest layer that the mask selects, which
    # the sweep down reads for both.
    cloudy = lo < nlev
    @inbounds begin
        albedo_clr = albedo_all = sfc_alb_diffuse
        src_clr = dir_clr[wcol, 1] * sfc_alb_direct
        up_clr[wcol, 1], dn_clr[wcol, 1] = albedo_clr, src_clr
        if cloudy
            src_all = dir_all[wcol, 1] * sfc_alb_direct
            dn_all[wcol, 1] = src_all
        else
            src_all = src_clr
        end
    end
    @inbounds for ilay in 1i32:nlay
        τ, ssa, g = a_clr[wcol, ilay], b_clr[wcol, ilay], a_all[wcol, ilay]
        (Rdir, Tdir, _, Rdif, Tdif) = sw_2stream_coeffs(τ, ssa, g, μ₀)
        dir = dir_c = dir_clr[wcol, ilay + 1i32]
        src_up, src_dn = Rdir * dir, Tdir * dir
        a_clr[wcol, ilay], b_clr[wcol, ilay], albedo_clr, src_clr =
            adding_layer(Rdif, Tdif, src_up, src_dn, albedo_clr, src_clr)
        up_clr[wcol, ilay + 1i32], dn_clr[wcol, ilay + 1i32] =
            albedo_clr, src_clr
        cloudy || continue
        if !isnothing(lkp_cld) && masked(mask, ilay, wcol, bounds)
            k = BAND_CLOUD
            τ, ssa, g = Optics.increment_2stream(
                τ,
                ssa,
                g,
                band[ilay, k + 1i32, gcol],
                band[ilay, k + 2i32, gcol],
                band[ilay, k + 3i32, gcol],
            )
            (Rdir, Tdir, _, Rdif, Tdif) = sw_2stream_coeffs(τ, ssa, g, μ₀)
        end
        dir = ilay + 1i32 ≤ hi ? dir_all[wcol, ilay + 1i32] : dir_c
        src_up, src_dn = Rdir * dir, Tdir * dir
        a, b_all[wcol, ilay, 3i32], albedo_all, src_all =
            adding_layer(Rdif, Tdif, src_up, src_dn, albedo_all, src_all)
        ilay ≥ lo && (a_all[wcol, ilay] = a)
        ilay ≥ lo && (up_all[wcol, ilay + 1i32] = albedo_all)
        dn_all[wcol, ilay + 1i32] = src_all
    end
    # Eq 12 at the top of the domain, where the diffuse flux down is zero;
    # `sweep_down_sums!` sweeps down
    @inbounds begin
        up_clr[wcol, nlev] = FT(0) * albedo_clr + src_clr
        dn_clr[wcol, nlev] = FT(0) + flux_dn_dir_top
        cloudy && (up_all[wcol, nlev] = FT(0) * albedo_all + src_all)
    end
    return lo, hi
end

# Keep `h` of the values `v` that each lane holds and add those that the lane
# `bit` apart sends: a step of a reduce-scatter by recursive halving (none if `h`
# is zero).
@inline function halve_group(
    v::NTuple{M, T},
    ::Val{h},
    ::Val{bit},
    lanes::Val{N},
    lane,
) where {M, T, h, bit, N}
    h == 0 && return v
    upper = (lane - 1i32) & Int32(bit) != 0i32
    mask = column_lanes(lanes)
    return ntuple(Val(h)) do i
        keep, send = upper ? (v[h + i], v[i]) : (v[i], v[h + i])
        keep + CUDA.shfl_xor_sync(mask, send, bit % UInt32, N % UInt32)
    end
end

# The sum over the `N` lanes of a column of one of the `G` values that each lane
# holds (`G ≤ N`, both powers of two): a reduce-scatter by recursive halving over
# the top bits of the lane, then a butterfly over the rest, which leaves the sum
# of value `((lane - 1) >> log2(N ÷ G)) + 1` in the lane. Each value is summed by
# lanes `N/2`, `N/4`, ..., `1` apart in turn, whatever `G` is.
@inline function group_sum(v::NTuple{G, T}, lanes::Val{N}, lane) where {G, N, T}
    Base.Cartesian.@nexprs 5 s ->
        (v = halve_group(v, Val(G >> s), Val(N >> s), lanes, lane))
    x = only(v)
    mask = column_lanes(lanes)
    Base.Cartesian.@nexprs 5 s -> begin
        b_s = (N ÷ G) >> s
        if b_s > 0
            x += CUDA.shfl_xor_sync(mask, x, b_s % UInt32, N % UInt32)
        end
    end
    return x
end

# Add the sums over a column's lanes of quantities `q` of `vals`, where
# `vals[k][q]` is that quantity at level `ilev - k + 1` in this lane, to `outs[q]`
# (or overwrite them on the first round).
@inline function add_group_sums!(
    outs::NTuple{Q},
    vals::NTuple{G},
    gcol,
    ilev,
    first,
    lanes::Val{N},
    lane,
) where {Q, G, N}
    shift = trailing_zeros(N) - trailing_zeros(G)
    level = ilev - ((lane - 1i32) >> shift)
    leader = (lane - 1i32) & Int32((1 << shift) - 1) == 0i32
    ntuple(Val(Q)) do q
        @inline
        s = group_sum(ntuple(k -> vals[k][q], Val(G)), lanes, lane)
        if leader
            out = outs[q]
            @inbounds out[gcol, level] = first ? s : out[gcol, level] + s
        end
        return nothing
    end
    return nothing
end

# The longwave fluxes at level `ilev` of both skies: the upward and downward flux
# of the clear sky and of the cloudy sky, from the coefficients, albedos and
# sources that `lw_stream!` left (Eqs 13 and 12), and the downward fluxes, which
# the next level down continues from. Below `lo`, the lowest layer that the
# g-point's cloud mask selects, the cloudy sky's coefficients, albedos and
# sources are the clear sky's, which `lw_stream!` left for both. An inactive
# lane has zero fluxes.
@inline function lw_level_fluxes(
    planes,
    wcol,
    ilev,
    nlev,
    flux_dn_clr,
    flux_dn_all,
    (active, lo, _),
)
    a_clr, b_clr, up_clr, dn_clr, a_all, b_all, up_all, dn_all = planes
    FT = typeof(flux_dn_clr)
    active || return (FT(0), FT(0), FT(0), FT(0)), flux_dn_clr, flux_dn_all
    @inbounds if ilev == nlev
        flux_dn_clr = flux_dn_all = dn_clr[wcol, ilev]
        up_c = up_clr[wcol, ilev]
        up_a = ilev > lo ? up_all[wcol, ilev] : up_c
        fluxes = (up_c, flux_dn_clr, up_a, flux_dn_all)
    else
        a_c, b_c = a_clr[wcol, ilev], b_clr[wcol, ilev]
        up_c, dn_c = up_clr[wcol, ilev], dn_clr[wcol, ilev]
        a_a, b_a =
            ilev ≥ lo ? (a_all[wcol, ilev], b_all[wcol, ilev]) : (a_c, b_c)
        up_a, dn_a =
            ilev > lo ? (up_all[wcol, ilev], dn_all[wcol, ilev]) : (up_c, dn_c)
        flux_dn_clr = a_c * flux_dn_clr + b_c
        flux_dn_all = a_a * flux_dn_all + b_a
        fluxes = (
            flux_dn_clr * up_c + dn_c,
            flux_dn_clr,
            flux_dn_all * up_a + dn_a,
            flux_dn_all,
        )
    end
    return fluxes, flux_dn_clr, flux_dn_all
end

# The shortwave fluxes at level `ilev` of both skies: the upward, downward and
# direct downward flux of the clear sky and of the cloudy sky (see
# `lw_level_fluxes`); the downward fluxes that the next level continues from are
# the diffuse ones, which are zero at the top. Where `sw_stream!` left the
# cloudy sky's values to the clear sky -- everywhere for a g-point that the
# cloud mask leaves clear, the albedos and coefficients `a` below `lo`, the
# lowest layer that the mask selects, and the direct beam above `hi`, the highest
# -- they are read from the clear sky's planes.
@inline function sw_level_fluxes(
    planes,
    wcol,
    ilev,
    nlev,
    flux_dn_clr,
    flux_dn_all,
    (active, lo, hi),
)
    a_clr,
    b_clr,
    up_clr,
    dn_clr,
    dir_clr,
    a_all,
    b_all,
    up_all,
    dn_all,
    dir_all = planes
    FT = typeof(flux_dn_clr)
    active || return ntuple(_ -> FT(0), Val(6)), flux_dn_clr, flux_dn_all
    cloudy = lo < nlev
    @inbounds if ilev == nlev
        up_c, dn_c = up_clr[wcol, ilev], dn_clr[wcol, ilev]
        dir_c = dir_clr[wcol, ilev]
        up_a = cloudy ? up_all[wcol, ilev] : up_c
        fluxes = (up_c, dn_c, dir_c, up_a, dn_c, dir_c)
    else
        a_c, b_c = a_clr[wcol, ilev], b_clr[wcol, ilev]
        up_c, dn_c = up_clr[wcol, ilev], dn_clr[wcol, ilev]
        dir_c = dir_clr[wcol, ilev]
        a_a = ilev ≥ lo ? a_all[wcol, ilev] : a_c
        b_a, dn_a =
            cloudy ? (b_all[wcol, ilev, 3i32], dn_all[wcol, ilev]) : (b_c, dn_c)
        up_a = ilev > lo ? up_all[wcol, ilev] : up_c
        dir_a = ilev ≤ hi ? dir_all[wcol, ilev] : dir_c
        flux_dn_clr = a_c * flux_dn_clr + b_c
        flux_dn_all = a_a * flux_dn_all + b_a
        fluxes = (
            flux_dn_clr * up_c + dn_c,
            flux_dn_clr + dir_c,
            dir_c,
            flux_dn_all * up_a + dn_a,
            flux_dn_all + dir_a,
            dir_a,
        )
    end
    return fluxes, flux_dn_clr, flux_dn_all
end

# Levels whose fluxes a column's lanes sum together in registers (see
# `group_sum`): four in the longwave, and two in the shortwave, whose six fluxes
# a level would otherwise spill at its register cap (10% faster on an A100 with
# the flagship AMIP configuration, where four are fastest in the longwave).
@inline group_levels(::typeof(lw_level_fluxes), ::Val{N}) where {N} = min(N, 4)
@inline group_levels(::typeof(sw_level_fluxes), ::Val{N}) where {N} = min(N, 2)

# The fluxes of `G` levels down from `ilev` (see `lw_level_fluxes`), of a lane
# in state `lane_state`.
@inline function level_group(
    level_fluxes::F,
    ::Val{G},
    planes,
    wcol,
    ilev,
    nlev,
    flux_dn_clr,
    flux_dn_all,
    lane_state,
) where {F, G}
    Base.Cartesian.@nexprs 4 k -> begin
        if k ≤ G
            v_k, flux_dn_clr, flux_dn_all = level_fluxes(
                planes,
                wcol,
                ilev - Int32(k - 1),
                nlev,
                flux_dn_clr,
                flux_dn_all,
                lane_state,
            )
        else
            v_k = v_1
        end
    end
    vs = (v_1, v_2, v_3, v_4)
    return ntuple(k -> vs[k], Val(G)), flux_dn_clr, flux_dn_all
end

# The sweep down of both skies from the top level (see `lw_level_fluxes` and
# `sw_level_fluxes`, which take the lane's state `lane_state`), with the fluxes
# of each level summed over the lanes of the column into `outs` as the sweep
# goes, a group of levels at a time.
@inline function sweep_down_sums!(
    level_fluxes::F,
    outs,
    planes,
    lane_state,
    gcol,
    wcol,
    nlay,
    first,
    flux_dn_top,
    lanes::Val{N},
    lane,
) where {F, N}
    planes = map(streamed, planes)
    nlev = nlay + 1i32
    n = group_levels(level_fluxes, lanes)
    G, g = Val(n), Int32(n)
    flux_dn_clr = flux_dn_all = flux_dn_top
    ilev = nlev
    while ilev ≥ g
        vals, flux_dn_clr, flux_dn_all = level_group(
            level_fluxes,
            G,
            planes,
            wcol,
            ilev,
            nlev,
            flux_dn_clr,
            flux_dn_all,
            lane_state,
        )
        add_group_sums!(outs, vals, gcol, ilev, first, lanes, lane)
        ilev -= g
    end
    while ilev ≥ 1i32
        vals, flux_dn_clr, flux_dn_all = level_group(
            level_fluxes,
            Val(1),
            planes,
            wcol,
            ilev,
            nlev,
            flux_dn_clr,
            flux_dn_all,
            lane_state,
        )
        add_group_sums!(outs, vals, gcol, ilev, first, lanes, lane)
        ilev -= 1i32
    end
    return nothing
end

# The per-level planes that a sweep keeps for its sweep down are written once and
# read once, far apart, so `Streamed` accesses them with the cache-streaming
# (evict-first) qualifier, which keeps them from evicting the lookup tables from
# L2 (6% faster on an A100 with the flagship AMIP configuration).
for (T, ty, reg) in ((Float32, "float", "f"), (Float64, "double", "d"))
    sfx = T === Float32 ? "f32" : "f64"
    @eval @inline function st_cs(p::Core.LLVMPtr{$T, 1}, v::$T)
        Base.llvmcall(
            (
                $"""
               define void @entry(ptr addrspace(1) %p, $ty %v) #0 {
                   %i = ptrtoint ptr addrspace(1) %p to i64
                   call void asm sideeffect "st.global.cs.$sfx [\$0], \$1;", "l,$reg,~{memory}"(i64 %i, $ty %v)
                   ret void
               }
               attributes #0 = { alwaysinline }
               """,
                "entry",
            ),
            Cvoid,
            Tuple{Core.LLVMPtr{$T, 1}, $T},
            p,
            v,
        )
    end
    @eval @inline function ld_cs(p::Core.LLVMPtr{$T, 1})
        Base.llvmcall(
            (
                $"""
               define $ty @entry(ptr addrspace(1) %p) #0 {
                   %i = ptrtoint ptr addrspace(1) %p to i64
                   %v = call $ty asm sideeffect "ld.global.cs.$sfx \$0, [\$1];", "=$reg,l,~{memory}"(i64 %i)
                   ret $ty %v
               }
               attributes #0 = { alwaysinline }
               """,
                "entry",
            ),
            $T,
            Tuple{Core.LLVMPtr{$T, 1}},
            p,
        )
    end
end
@inline element_pointer(
    A::Int32Array{T, N, <:CUDA.CuDeviceArray},
    I...,
) where {T, N} = pointer(A.array, linear_index(A, I))
@inline element_pointer(A::SubArray, I...) =
    element_pointer(parent(A), Base.reindex(A.indices, I)...)
struct Streamed{T, N, A <: AbstractArray{T, N}} <: AbstractArray{T, N}
    array::A
end
Base.size(s::Streamed) = size(s.array)
@inline Base.getindex(s::Streamed, I::Integer...) =
    ld_cs(element_pointer(s.array, I...))
@inline Base.setindex!(s::Streamed{T}, v, I::Integer...) where {T} =
    (st_cs(element_pointer(s.array, I...), convert(T, v)); s)
const DevicePlane{T} = Union{
    Int32Array{T, <:Any, <:CUDA.CuDeviceArray},
    SubArray{T, <:Any, <:Int32Array{T, <:Any, <:CUDA.CuDeviceArray}},
}
@inline streamed(A::DevicePlane{T}) where {T} =
    Streamed{T, ndims(A), typeof(A)}(A)
@inline streamed(A) = A # elsewhere, such as an emulation on the CPU

# The first and last cloudy layers of this thread's column, which every g-point's
# cloud mask shares (see `Optics.cloud_layer_bounds`), or zeros.
@inline cloud_bounds(::Nothing, gcol) = (0i32, 0i32)
@inline function cloud_bounds(cloud_state, gcol)
    start, finish =
        Optics.cloud_layer_bounds(view(cloud_state.cld_frac, :, gcol))
    return (start % Int32, finish % Int32)
end
# Whether the cloud mask selects layer `ilay`: it is set only within the bounds.
@inline masked(mask, ilay, wcol, (start, finish)) =
    start ≤ ilay ≤ finish && @inbounds mask[ilay, wcol]

# The lanes' cloud masks, kept as bits in the block's dynamic shared memory (in
# place of a byte for each layer in global memory, which a sweep would read back
# from L2): a row of `words` for each thread, bit `(ilay - 1) % 32` of word
# `(ilay - 1) ÷ 32 + 1` for layer `ilay`. A thread reads and writes only its own
# row, so `view(mask, :, wcol)` and `mask[ilay, wcol]` ignore `wcol`.
struct LaneMasks{A}
    words::A
end
struct LaneMask{A} <: AbstractVector{Bool}
    words::A
    tid::Int32
end
@inline lane_masks(nlay) = LaneMasks(
    CUDA.CuDynamicSharedArray(UInt32, (Int(blockDim().x), Int(cld(nlay, 32)))),
)
Base.size(m::LaneMask) = (32 * size(m.words, 2),)
@inline Base.view(m::LaneMasks, ::Colon, wcol) =
    LaneMask(m.words, threadIdx().x)
@inline Base.getindex(m::LaneMasks, ilay::Integer, wcol::Integer) =
    view(m, :, wcol)[ilay]
@inline mask_bit(ilay) =
    (Int32((ilay - 1) >> 5) + 1i32, one(UInt32) << (UInt32(ilay - 1) & 0x1f))
@inline function Base.getindex(m::LaneMask, ilay::Integer)
    w, b = mask_bit(ilay)
    return @inbounds (m.words[m.tid, w] & b) != zero(UInt32)
end
@inline function Base.setindex!(m::LaneMask, x::Bool, ilay::Integer)
    w, b = mask_bit(ilay)
    @inbounds word = m.words[m.tid, w]
    @inbounds m.words[m.tid, w] = x ? word | b : word & ~b
    return m
end

# The fused solve of either spectrum (see `_solve_both_skies_grouped!`). A night
# column (in the shortwave) skips its optics and sweeps and is zeroed at the end,
# but still draws its cloud masks for the cloud cover, and its 550 nm band still
# collects the aerosol optical depth.
function fused_CUDA!(
    fluxb,
    src,
    bcs,
    op::TwoStream,
    band,
    flux,
    flux_clear,
    nlay,
    ncol,
    lanes::Val{N},
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup::Union{LookUpLW, LookUpSW},
    lookup_cld::Union{LookUpCld, Nothing},
    lookup_aero::Union{LookUpAerosolMerra, Nothing},
    lookup_precip::Union{LookUpPrecip, Nothing},
) where {N}
    nlanes = Int32(N)
    wcol = threadIdx().x + (blockIdx().x - 1i32) * blockDim().x
    gcol = ((wcol - 1i32) >> trailing_zeros(nlanes)) + 1i32
    lane = ((wcol - 1i32) & (nlanes - 1i32)) + 1i32
    if gcol ≤ ncol
        nlev = nlay + 1i32
        (; major_gpt2bnd, bnd_lims_gpt) = lookup.band_data
        n_gpt = length(major_gpt2bnd) % Int32
        (; cloud_state, aerosol_state) = as
        FT = eltype(flux.flux_up)
        mask = lane_masks(nlay)
        sw = lookup isa LookUpSW # delta-scaled clouds and aerosols, and the AOD
        day = sunlit(bcs, gcol)
        aod_band =
            sw && !isnothing(lookup_aero) ? lookup_aero.iband_550nm % Int32 :
            0i32
        lane == 1i32 && _compute_aero_mask!(aerosol_state, gcol)
        if !isnothing(lookup_aero)
            CUDA.sync_warp(column_lanes(lanes))
            fill_aerosol_searches!(
                band,
                as,
                state_cache,
                lookup_aero,
                gcol,
                lane,
                lanes,
                nlay,
            )
        end
        n_cloudy = 0i32
        bounds = cloud_bounds(cloud_state, gcol)
        first = true
        @inbounds for ibnd in 1i32:(size(bnd_lims_gpt, 2) % Int32)
            CUDA.sync_warp(column_lanes(lanes))
            if !isnothing(lookup_aero) && (day || ibnd == aod_band)
                fill_band_aerosol!(
                    band,
                    as,
                    lookup_aero,
                    ibnd,
                    gcol,
                    lane,
                    lanes,
                    nlay,
                    sw,
                )
            end
            if day
                fill_band_gas!(
                    band,
                    as,
                    state_cache,
                    lookup,
                    ibnd,
                    gcol,
                    lane,
                    lanes,
                    nlay,
                )
                isnothing(lookup_cld) || fill_band_clouds!(
                    band,
                    as,
                    lookup_cld,
                    lookup_precip,
                    ibnd,
                    gcol,
                    lane,
                    lanes,
                    nlay,
                    sw,
                )
            end
            CUDA.sync_warp(column_lanes(lanes))
            if lane == 1i32 && ibnd == aod_band
                aerosol_state.aod_sw_ext[gcol], aerosol_state.aod_sw_sca[gcol] =
                    aerosol_optical_depths(
                        band,
                        aerosol_state.aero_mask,
                        gcol,
                        nlay,
                    )
            end
            gpt_first = bnd_lims_gpt[1, ibnd] % Int32
            gpt_last = bnd_lims_gpt[2, ibnd] % Int32
            for igpt_first in gpt_first:nlanes:gpt_last
                igpt = igpt_first + lane - 1i32
                active = igpt ≤ gpt_last
                lo, hi = nlev, 0i32
                if active
                    # A night column only counts its cloudy g-points.
                    cloudy =
                        day ?
                        _build_cloud_mask!(
                            cloud_state,
                            mask_name(lookup),
                            gcol,
                            igpt,
                            mask,
                            wcol,
                            bounds,
                        ) :
                        _build_cloud_mask!(
                            cloud_state,
                            mask_name(lookup),
                            gcol,
                            igpt,
                            mask,
                            wcol,
                            bounds,
                            Val(false),
                        )
                    n_cloudy += Int32(cloudy)
                    if day
                        lo, hi = stream!(
                            op,
                            fluxb,
                            src,
                            bcs,
                            mask,
                            band,
                            as,
                            state_cache,
                            lookup,
                            lookup_cld,
                            lookup_aero,
                            gcol,
                            wcol,
                            igpt,
                            ibnd,
                            igpt - gpt_first + 1i32,
                            bounds,
                            nlay,
                        )
                    end
                end
                day && sweep_down_sums!(
                    level_fluxes(src),
                    down_outs(flux, flux_clear),
                    down_planes(op, fluxb, src),
                    (active, lo, hi),
                    gcol,
                    wcol,
                    nlay,
                    first,
                    zero(FT),
                    lanes,
                    lane,
                )
                first = false
            end
        end
        n_cloudy = lane_sum(n_cloudy, lanes)
        if lane == 1i32
            if day
                @inbounds compute_net_flux!(flux, gcol, nlev)
                @inbounds compute_net_flux!(flux_clear, gcol, nlev)
            else
                @inbounds set_flux_to_zero!(flux, gcol, nlev)
                @inbounds set_flux_to_zero!(flux_clear, gcol, nlev)
            end
            cover = cloud_cover(cloud_state, lookup)
            isnothing(cover) || @inbounds cover[gcol] = FT(n_cloudy) / n_gpt
        end
    end
    return nothing
end

# What sets the two spectra's fused solves apart.
@inline sunlit(::LwBCs, gcol) = true
@inline sunlit(bcs::SwBCs, gcol) = @inbounds bcs.cos_zenith[gcol] > 0
@inline mask_name(::LookUpLW) = Val(:mask_lw)
@inline mask_name(::LookUpSW) = Val(:mask_sw)
@inline cloud_cover(cs, lookup) = nothing
@inline cloud_cover(cs::CloudState, ::LookUpLW) = cs.cld_cover_lw
@inline cloud_cover(cs::CloudState, ::LookUpSW) = cs.cld_cover_sw
@inline stream!(op, fluxb, src::SourceLW2Str, args...) =
    (lw_stream!(op, fluxb, src, args...), 0i32)
@inline stream!(op, fluxb, src::SourceSW2Str, args...) =
    sw_stream!(op, fluxb, src, args...)
@inline level_fluxes(::SourceLW2Str) = lw_level_fluxes
@inline level_fluxes(::SourceSW2Str) = sw_level_fluxes
# The fluxes that the sweep down sums into, clear sky first, and the planes that
# `lw_stream!` and `sw_stream!` leave for it.
@inline down_outs(flux::FluxLW, flux_clear::FluxLW) =
    (flux_clear.flux_up, flux_clear.flux_dn, flux.flux_up, flux.flux_dn)
@inline down_outs(flux::FluxSW, flux_clear::FluxSW) = (
    flux_clear.flux_up,
    flux_clear.flux_dn,
    flux_clear.flux_dn_dir,
    flux.flux_up,
    flux.flux_dn,
    flux.flux_dn_dir,
)
@inline down_planes(op, fluxb::FluxLW, src::SourceLW2Str) = (
    fluxb.flux_up,
    fluxb.flux_dn,
    src.albedo,
    src.src,
    op.τ,
    op.ssa,
    fluxb.flux_net,
    src.lev_source,
)
@inline down_planes(op, fluxb::FluxSW, src::SourceSW2Str) = (
    op.τ,
    op.ssa,
    src.albedo,
    src.src,
    fluxb.flux_dn_dir,
    op.g,
    src.leveldata,
    fluxb.flux_up,
    fluxb.flux_dn,
    fluxb.flux_net,
)
