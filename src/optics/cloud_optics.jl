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
)
    @inbounds begin
        nlay = length(τ)
        extliq, ssaliq, asyliq = LookUpTables.getview_liqdata(lkp_cld, ibnd)
        extice, ssaice, asyice =
            LookUpTables.getview_icedata(lkp_cld, ibnd, ice_rgh)
        _, _, nsize_liq, nsize_ice, _ = lkp_cld.dims
        radliq_lwr, radliq_upr, radice_lwr, radice_upr = lkp_cld.bounds

        for glay in 1:nlay
            if cld_mask[glay]
                # cloud liquid particles
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
                # cloud ice particles
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
                # add cloud optical optics
                τ[glay] += (τl - τl_ssa) + (τi - τi_ssa)
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
        delta_scaling = false,
    )
Compute the `TwoStream` cloud optics properties and add them
to the `TwoStream` gas optics properties.
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
    delta_scaling = false,
)
    nlay = length(τ)
    FT = eltype(τ)
    extliq, ssaliq, asyliq = LookUpTables.getview_liqdata(lkp_cld, ibnd)
    extice, ssaice, asyice =
        LookUpTables.getview_icedata(lkp_cld, ibnd, ice_rgh)
    _, _, nsize_liq, nsize_ice, _ = lkp_cld.dims
    radliq_lwr, radliq_upr, radice_lwr, radice_upr = lkp_cld.bounds
    @inbounds begin
        for glay in 1:nlay
            if cld_mask[glay]
                # cloud liquid particles
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
                # cloud ice particles
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
                ssa_cl = τl_ssa + τi_ssa
                g_cl = (τl_ssag + τi_ssag) / max(eps(FT), ssa_cl)
                ssa_cl /= max(eps(FT), τ_cl)

                if delta_scaling # delta scaling is applied for shortwave problem
                    τ_cl, ssa_cl, g_cl = delta_scale(τ_cl, ssa_cl, g_cl)
                end
                # compute cloud optical optics
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

Compute the `TwoStream` cloud liquid properties using the `LookUpTable` method.
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

# A uniform draw in [0, 1) determined entirely by the key and where it is used,
# with no RNG state and no allocation: the same (key, column, g-point, layer)
# always yields the same number, on any device, in any launch order, at any
# thread count. The key comes from the `seedval` passed to `update_fluxes!`, so
# it varies between radiation steps and repeats exactly on a rerun or restart.
#
# This is what makes McICA sampling reproducible on the GPU. The device RNG is
# keyed per kernel launch, so drawing from it makes the sample depend on how many
# kernels ran rather than on the physics (RRTMGP.jl#316). Mixing is splitmix64's
# finalizer, a few integer operations -- so the reason build_cloud_mask! gave for
# not doing this, keeping the kernel allocation-free, does not apply.
#
# Internal, and deliberately not a docstring: Documenter runs with checkdocs, so
# a docstring here would have to be carried in the public manual.
@inline function _keyed_uniform(
    ::Type{FT},
    key::Unsigned,
    gcol::Integer,
    igpt::Integer,
    ilay::Integer,
) where {FT}
    z = (UInt64(key) * 0xD1342543DE82EF95) ⊻
        (UInt64(gcol) * 0x9E3779B97F4A7C15) ⊻
        (UInt64(igpt) * 0xBF58476D1CE4E5B9) ⊻
        (UInt64(ilay) * 0x94D049BB133111EB)
    z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
    z = z ⊻ (z >> 31)
    # Top 24 bits, so the result is exact in Float32 as well as Float64
    return FT(z >> 40) * FT(1 / (1 << 24))
end

"""
    build_cloud_mask!(cld_mask, cld_frac, ::MaxRandomOverlap, key, gcol, igpt)

Build a McICA-sampled cloud mask from cloud fraction data for maximum-random overlap

Reference: https://github.com/AER-RC/RRTMG_SW/

Determinism: the mask is drawn by `_keyed_uniform` from `(key, gcol, igpt, ilay)`, so it
is reproducible everywhere -- one CPU thread, many threads, or the GPU -- and independent of
launch order, launch count and work assignment. Nothing is seeded and nothing is allocated.

This replaces a bare `Random.rand()`, which took the device RNG on the GPU. That RNG is keyed
per kernel launch, so the sample depended on how many kernels had run: two identical solves
disagreed, `reset_rng_seed` had no effect on the GPU, and restarts were not reproducible
(RRTMGP.jl#316, #544).

The key is `cloud_state.seed`, set from the `seedval` given to `update_fluxes!` -- hosts
typically pass the timestep index, which gives a fresh sample every radiation step and the same
sample on a rerun. A host that passes no `seedval` gets a random key per call, i.e. the previous
unreproducible behavior.
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
