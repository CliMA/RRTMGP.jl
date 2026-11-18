"""
    add_cloud_optics_1scalar!(
        τ,
        cld_mask,
        cld_r_eff_liq,
        cld_r_eff_ice,
        cld_path_liq,
        cld_path_ice,
        ice_rgh,
        lkp_cld,
        ibnd;
        lkp_precip = nothing,
        cld_r_eff_sno = nothing,
        cld_path_sno = nothing,
    )

Compute the `OneScalar` (longwave) cloud and snow optical depths and add them to
the gas optical depth in layers where `cld_mask` is `true` (sampled from
`cloud_fraction`). When a snow water path is supplied (`cld_path_sno`
non-`nothing`), `lkp_precip::LookUpPrecip` must also be provided; snow optical
depth is evaluated from `lkp_precip` (`snodata`, log-uniform in effective radius
from `radsno_lwr` to `radsno_upr`) via [`compute_lookup_precip_props`](@ref).
Rain optics are not yet wired here.
"""
@inline function add_cloud_optics_1scalar!(
    τ,
    cld_mask,
    cld_r_eff_liq,
    cld_r_eff_ice,
    cld_path_liq,
    cld_path_ice,
    ice_rgh,
    lkp_cld::LookUpCld,
    ibnd;
    lkp_precip::Union{LookUpPrecip, Nothing} = nothing,
    cld_r_eff_sno = nothing,
    cld_path_sno = nothing,
)
    use_snow = !isnothing(cld_path_sno)
    if use_snow && isnothing(lkp_precip)
        error(
            "snow optics requested (cld_path_sno !== nothing), but no LookUpPrecip lookup table was provided",
        )
    end
    @inbounds begin
        nlay = length(τ)
        extliq, ssaliq, asyliq = LookUpTables.getview_liqdata(lkp_cld, ibnd)
        extice, ssaice, asyice =
            LookUpTables.getview_icedata(lkp_cld, ibnd, ice_rgh)
        _, _, nsize_liq, nsize_ice, _ = lkp_cld.dims
        radliq_lwr, radliq_upr, radice_lwr, radice_upr = lkp_cld.bounds

        for glay in 1:nlay
            if cld_mask[glay]
                τl, τl_ssa, _ = compute_lookup_cld_liq_props(
                    nsize_liq,
                    radliq_lwr,
                    radliq_upr,
                    extliq,
                    ssaliq,
                    asyliq,
                    cld_r_eff_liq[glay],
                    cld_path_liq[glay],
                )
                τi, τi_ssa, _ = compute_lookup_cld_ice_props(
                    nsize_ice,
                    radice_lwr,
                    radice_upr,
                    extice,
                    ssaice,
                    asyice,
                    cld_r_eff_ice[glay],
                    cld_path_ice[glay],
                )
                τ[glay] += (τl - τl_ssa) + (τi - τi_ssa)
                if use_snow
                    extsno, ssasno, asysno =
                        LookUpTables.getview_snodata(lkp_precip, ibnd)
                    nsize_sno = LookUpTables.get_nsize_sno(lkp_precip)
                    radsno_lwr, radsno_upr =
                        lkp_precip.bounds[1], lkp_precip.bounds[2]
                    τs, τs_ssa, _ = compute_lookup_precip_props(
                        nsize_sno,
                        radsno_lwr,
                        radsno_upr,
                        extsno,
                        ssasno,
                        asysno,
                        cld_r_eff_sno[glay],
                        cld_path_sno[glay],
                    )
                    τ[glay] += (τs - τs_ssa)
                end
            end
        end
    end
    return nothing
end

"""
    add_cloud_optics_2stream!(
        τ,
        ssa,
        g,
        cld_mask,
        cld_r_eff_liq,
        cld_r_eff_ice,
        cld_path_liq,
        cld_path_ice,
        ice_rgh,
        lkp_cld,
        ibnd;
        lkp_precip = nothing,
        cld_r_eff_sno = nothing,
        cld_path_sno = nothing,
        delta_scaling = false,
    )

Compute the `TwoStream` cloud and snow optical properties and add them to the
`TwoStream` gas/aerosol optical properties in layers where `cld_mask` is `true`
(sampled from `cloud_fraction`). When a snow water path is supplied
(`cld_path_sno` non-`nothing`), `lkp_precip::LookUpPrecip` must also be
provided; snow optical properties (`τ`, `τ*ssa`, `τ*ssa*g`) are computed from
`lkp_precip` (`snodata`, log-uniform in effective radius from `radsno_lwr` to
`radsno_upr`) via [`compute_lookup_precip_props`](@ref) and combined with cloud
liquid and ice before delta-scaling. Rain optics are not yet wired here.
"""
@inline function add_cloud_optics_2stream!(
    τ,
    ssa,
    g,
    cld_mask,
    cld_r_eff_liq,
    cld_r_eff_ice,
    cld_path_liq,
    cld_path_ice,
    ice_rgh,
    lkp_cld::LookUpCld,
    ibnd;
    lkp_precip::Union{LookUpPrecip, Nothing} = nothing,
    cld_r_eff_sno = nothing,
    cld_path_sno = nothing,
    delta_scaling = false,
)
    use_snow = !isnothing(cld_path_sno)
    if use_snow && isnothing(lkp_precip)
        error(
            "snow optics requested (cld_path_sno !== nothing), but no LookUpPrecip lookup table was provided",
        )
    end
    nlay = length(τ)
    FT = eltype(τ)
    extliq, ssaliq, asyliq = LookUpTables.getview_liqdata(lkp_cld, ibnd)
    extice, ssaice, asyice =
        LookUpTables.getview_icedata(lkp_cld, ibnd, ice_rgh)
    @inbounds _, _, nsize_liq, nsize_ice, _ = lkp_cld.dims
    @inbounds radliq_lwr, radliq_upr, radice_lwr, radice_upr = lkp_cld.bounds

    @inbounds begin
        for glay in 1:nlay
            if cld_mask[glay]
                # cloud liquid + ice
                τl, τl_ssa, τl_ssag = compute_lookup_cld_liq_props(
                    nsize_liq,
                    radliq_lwr,
                    radliq_upr,
                    extliq,
                    ssaliq,
                    asyliq,
                    cld_r_eff_liq[glay],
                    cld_path_liq[glay],
                )
                τi, τi_ssa, τi_ssag = compute_lookup_cld_ice_props(
                    nsize_ice,
                    radice_lwr,
                    radice_upr,
                    extice,
                    ssaice,
                    asyice,
                    cld_r_eff_ice[glay],
                    cld_path_ice[glay],
                )
                τ_cl = τl + τi
                ssa_cl = τl_ssa + τi_ssa      # running Σ τ*ssa
                ssag_cl = τl_ssag + τi_ssag   # running Σ τ*ssa*g
                if use_snow
                    extsno, ssasno, asysno =
                        LookUpTables.getview_snodata(lkp_precip, ibnd)
                    nsize_sno = LookUpTables.get_nsize_sno(lkp_precip)
                    radsno_lwr, radsno_upr =
                        lkp_precip.bounds[1], lkp_precip.bounds[2]
                    τs, τs_ssa, τs_ssag = compute_lookup_precip_props(
                        nsize_sno,
                        radsno_lwr,
                        radsno_upr,
                        extsno,
                        ssasno,
                        asysno,
                        cld_r_eff_sno[glay],
                        cld_path_sno[glay],
                    )
                    τ_cl += τs
                    ssa_cl += τs_ssa
                    ssag_cl += τs_ssag
                end
                # g = Σ(τ*ssa*g) / Σ(τ*ssa); ssa = Σ(τ*ssa) / Σ(τ)  (same order as before)
                g_cl = ssag_cl / max(eps(FT), ssa_cl)
                ssa_cl /= max(eps(FT), τ_cl)
                if delta_scaling # shortwave only
                    τ_cl, ssa_cl, g_cl = delta_scale(τ_cl, ssa_cl, g_cl)
                end
                τ[glay], ssa[glay], g[glay] = increment_2stream(
                    τ[glay],
                    ssa[glay],
                    g[glay],
                    τ_cl,
                    ssa_cl,
                    g_cl,
                )
            end
        end
    end
    return nothing
end

"""
    compute_lookup_cld_liq_props(
        nsize_liq,
        radliq_lwr,
        radliq_upr,
        extliq,
        ssaliq,
        asyliq,
        re_liq,
        cld_path_liq,
    )

Compute `TwoStream` optical properties (τ, τ*ssa, τ*ssa*g) for cloud liquid via
the `LookUpTable` method.
"""
@inline function compute_lookup_cld_liq_props(
    nsize_liq,
    radliq_lwr,
    radliq_upr,
    extliq,
    ssaliq,
    asyliq,
    re_liq,
    cld_path_liq,
)
    FT = eltype(re_liq)
    τl, τl_ssa, τl_ssag = FT(0), FT(0), FT(0)
    # cloud liquid particles
    if cld_path_liq > eps(FT)
        Δr_liq = (radliq_upr - radliq_lwr) / FT(nsize_liq - 1)
        re_liq = max(min(re_liq, radliq_upr), radliq_lwr)
        loc = Int(
            max(
                min(
                    unsafe_trunc(Int, (re_liq - radliq_lwr) / Δr_liq) + 1,
                    nsize_liq - 1,
                ),
                1,
            ),
        )
        fac = (re_liq - radliq_lwr - (loc - 1) * Δr_liq) / Δr_liq
        fc1 = FT(1) - fac
        @inbounds begin
            τl = max(
                (fc1 * extliq[loc] + fac * extliq[loc + 1]) * cld_path_liq,
                FT(0),
            )
            τl_ssa = (fc1 * ssaliq[loc] + fac * ssaliq[loc + 1]) * τl
            τl_ssag = (fc1 * asyliq[loc] + fac * asyliq[loc + 1]) * τl_ssa
        end
    end
    return (τl, τl_ssa, τl_ssag)
end

"""
    compute_lookup_cld_ice_props(
        nsize_ice,
        radice_lwr,
        radice_upr,
        extice,
        ssaice,
        asyice,
        re_ice,
        cld_path_ice,
    )

Compute the `TwoStream` cloud ice properties using the `LookUpTable` method.
"""
@inline function compute_lookup_cld_ice_props(
    nsize_ice,
    radice_lwr,
    radice_upr,
    extice,
    ssaice,
    asyice,
    re_ice,
    cld_path_ice,
)
    FT = eltype(re_ice)
    τi, τi_ssa, τi_ssag = FT(0), FT(0), FT(0)
    # cloud ice particles
    if cld_path_ice > eps(FT)
        Δr_ice = (radice_upr - radice_lwr) / FT(nsize_ice - 1)
        re_ice = max(min(re_ice, radice_upr), radice_lwr)
        loc = Int(
            max(
                min(
                    unsafe_trunc(Int, (re_ice - radice_lwr) / Δr_ice) + 1,
                    nsize_ice - 1,
                ),
                1,
            ),
        )
        fac = (re_ice - radice_lwr - (loc - 1) * Δr_ice) / Δr_ice
        fc1 = FT(1) - fac
        @inbounds begin
            τi = max(
                (fc1 * extice[loc] + fac * extice[loc + 1]) * cld_path_ice,
                FT(0),
            )
            τi_ssa = (fc1 * ssaice[loc] + fac * ssaice[loc + 1]) * τi
            τi_ssag = (fc1 * asyice[loc] + fac * asyice[loc + 1]) * τi_ssa
        end
    end
    return (τi, τi_ssa, τi_ssag)
end

"""
    compute_lookup_precip_props(
        nsize,
        rad_lwr,
        rad_upr,
        ext_tbl,
        ssa_tbl,
        asy_tbl,
        re,
        precip_path,
    )

Compute `TwoStream` optical properties (`τ`, `τ*ssa`, `τ*ssa*g`) for a
precipitation species (such as snow, or rain when wired) with per-layer mass
path `precip_path` [g/m²], tabulated in [`LookUpPrecip`](@ref) on a log-uniform
effective-radius grid `r_k = rad_lwr * (rad_upr / rad_lwr)^((k - 1) / (nsize - 1))`.
Because `log(r_k)` is uniformly spaced with step
`Δlog_r = log(rad_upr / rad_lwr) / (nsize - 1)`, the table index and
interpolation weight are computed in `O(1)` without binary search.
"""
@inline function compute_lookup_precip_props(
    nsize,
    rad_lwr,
    rad_upr,
    ext_tbl,
    ssa_tbl,
    asy_tbl,
    re,
    precip_path,
)
    FT = eltype(re)
    τp, τp_ssa, τp_ssag = FT(0), FT(0), FT(0)
    if precip_path > eps(FT)
        log_lwr = log(rad_lwr)
        Δlog_r = (log(rad_upr) - log_lwr) / FT(nsize - 1)
        re_clamped = max(min(re, rad_upr), rad_lwr)
        log_re = log(re_clamped)
        loc = Int(
            max(
                min(
                    unsafe_trunc(Int, (log_re - log_lwr) / Δlog_r) + 1,
                    nsize - 1,
                ),
                1,
            ),
        )
        fac = (log_re - log_lwr - (loc - 1) * Δlog_r) / Δlog_r
        fc1 = FT(1) - fac
        @inbounds begin
            τp = max(
                (fc1 * ext_tbl[loc] + fac * ext_tbl[loc + 1]) * precip_path,
                FT(0),
            )
            τp_ssa = (fc1 * ssa_tbl[loc] + fac * ssa_tbl[loc + 1]) * τp
            τp_ssag = (fc1 * asy_tbl[loc] + fac * asy_tbl[loc + 1]) * τp_ssa
        end
    end
    return (τp, τp_ssa, τp_ssag)
end

# Uniform draw in [0, 1) from the key and the coordinates alone -- no RNG state,
# no allocation -- so the same (key, column, g-point, layer) always gives the
# same number, on any device and in any launch order.
#
# The four multipliers mix each coordinate into its own part of the word before
# they are combined; the two that follow are splitmix64's finalizer, which is
# what makes the combined word avalanche. They are from Steele, Lea and Flood,
# "Fast splittable pseudorandom number generators" (OOPSLA 2014), via
# https://prng.di.unimi.it/splitmix64.c -- chosen for their mixing behavior,
# with no arithmetic meaning beyond being odd. 0x9E3779B97F4A7C15 is 2^64/phi,
# the usual golden-ratio choice.
@inline function _keyed_uniform(
    ::Type{FT},
    key::Unsigned,
    gcol::Integer,
    igpt::Integer,
    ilay::Integer,
) where {FT}
    z =
        (UInt64(key) * 0xD1342543DE82EF95) ⊻
        (UInt64(gcol) * 0x9E3779B97F4A7C15) ⊻
        (UInt64(igpt) * 0xBF58476D1CE4E5B9) ⊻
        (UInt64(ilay) * 0x94D049BB133111EB)
    z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
    z = z ⊻ (z >> 31)
    # Use the mantissa the format actually has: 24 bits is exact in Float32,
    # and taking only 24 for Float64 would quantize cloud probability in steps
    # of 2^-24 and make fractions below ~6e-8 unreachable.
    return _uniform_from_bits(FT, z)
end

@inline _uniform_from_bits(::Type{Float32}, z::UInt64) =
    Float32(z >> 40) * Float32(1 / (1 << 24))
@inline function _uniform_from_bits(::Type{FT}, z::UInt64) where {FT}
    return FT(z >> 11) * FT(1 / (UInt64(1) << 53))
end

"""
    build_cloud_mask!(cld_mask, cld_frac, ::MaxRandomOverlap, key, gcol, igpt)

Build a McICA-sampled cloud mask from cloud fraction data for maximum-random overlap

Reference: https://github.com/AER-RC/RRTMG_SW/

Determinism: the mask comes from `_keyed_uniform` on `(key, gcol, igpt, ilay)`, so it is
reproducible on any device and independent of launch order and thread count. The key is
`cloud_state.mcica_key`, which `update_fluxes!` sets from its `seedval`; hosts pass the timestep
index, so each radiation step samples afresh. No `seedval` means a random key per call.
"""
function build_cloud_mask!(
    cld_mask::AbstractArray{Bool, 1},
    cld_frac::AbstractArray{FT, 1},
    ::MaxRandomOverlap,
    key::Unsigned,
    gcol::Integer,
    igpt::Integer,
) where {FT}
    nlay = size(cld_frac, 1)
    start = _get_start(cld_frac) # first cloudy layer

    if start > 0
        finish = _get_finish(cld_frac) # last cloudy layer
        # set cloud mask for non-cloudy layers
        _mask_outer_non_cloudy_layers!(cld_mask, start, finish)
        # RRTMG uses random_arr[finish] > (FT(1) - cld_frac[finish]),
        # we change > to >= to address edge cases
        @inbounds cld_frac_ilayplus1 = cld_frac[finish]
        random_ilayplus1 = _keyed_uniform(FT, key, gcol, igpt, finish)
        @inbounds cld_mask[finish] =
            cld_mask_ilayplus1 =
                random_ilayplus1 >= (FT(1) - cld_frac_ilayplus1)
        ilay = finish - 1
        while ilay ≥ start
            @inbounds cld_frac_ilay = cld_frac[ilay]
            if cld_frac_ilay > FT(0)
                # use same random number from the layer above if layer above is cloudy
                # update random numbers if layer above is not cloudy
                random_ilay =
                    cld_mask_ilayplus1 ? random_ilayplus1 :
                    _keyed_uniform(FT, key, gcol, igpt, ilay) *
                    (FT(1) - cld_frac_ilayplus1)
                # RRTMG uses random_arr[ilay] > (FT(1) - cld_frac[ilay]), we change > to >= to address edge cases
                cld_mask_ilay = random_ilay >= (FT(1) - cld_frac_ilay)
                random_ilayplus1 = random_ilay
            else
                cld_mask_ilay = false
            end
            @inbounds cld_mask[ilay] = cld_mask_ilay
            cld_frac_ilayplus1 = cld_frac_ilay
            cld_mask_ilayplus1 = cld_mask_ilay
            ilay -= 1
        end
    else
        map!(x -> false, cld_mask, cld_mask)
    end
    return nothing
end

function _get_finish(cld_frac)
    @inbounds for ilay in reverse(eachindex(cld_frac))
        cld_frac[ilay] > 0 && return ilay
    end
    return 0
end

function _get_start(cld_frac)
    @inbounds for ilay in eachindex(cld_frac)
        cld_frac[ilay] > 0 && return ilay
    end
    return 0
end

function _mask_outer_non_cloudy_layers!(cld_mask, start, finish)
    if start > 0
        for ilay in 1:(start - 1)
            @inbounds cld_mask[ilay] = false
        end
        nlay = length(cld_mask)
        for ilay in (finish + 1):nlay
            @inbounds cld_mask[ilay] = false
        end
    end
    return nothing
end
