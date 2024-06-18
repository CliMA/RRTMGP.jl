
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

This function computes the TwoStream clouds optics properties and adds them
to the TwoStream gas optics properties.
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
    lkp_cld::PadeCld,
    ibnd;
    delta_scaling = false,
)
    nlay = length(τ)
    @inbounds begin
        for glay in 1:nlay
            if cld_mask[glay]
                # compute cloud optical optics
                τ_cl, ssa_cl, g_cl = compute_pade_cld_props(
                    lkp_cld,
                    cld_r_eff_liq[glay],
                    cld_r_eff_ice[glay],
                    ice_rgh,
                    cld_path_liq[glay],
                    cld_path_ice[glay],
                    ibnd,
                )
                if delta_scaling # delta scaling is applied for shortwave problem
                    τ_cl, ssa_cl, g_cl = delta_scale(τ_cl, ssa_cl, g_cl)
                end
                # compute cloud optical optics
                τ[glay], ssa[glay], g[glay] = increment_2stream(τ[glay], ssa[glay], g[glay], τ_cl, ssa_cl, g_cl)
            end
        end
    end
    return nothing
end

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
    lut_extliq, lut_ssaliq, lut_asyliq = LookUpTables.getview_liqdata(lkp_cld, ibnd)
    lut_extice, lut_ssaice, lut_asyice = LookUpTables.getview_icedata(lkp_cld, ibnd, ice_rgh)
    _, _, nsize_liq, nsize_ice, _ = lkp_cld.dims
    radliq_lwr, radliq_upr, _, radice_lwr, radice_upr, _ = lkp_cld.bounds
    @inbounds begin
        for glay in 1:nlay
            if cld_mask[glay]
                # cloud liquid particles
                τl, τl_ssa, τl_ssag = compute_lookup_cld_liq_props(
                    nsize_liq,
                    radliq_lwr,
                    radliq_upr,
                    lut_extliq,
                    lut_ssaliq,
                    lut_asyliq,
                    cld_r_eff_liq[glay],
                    cld_path_liq[glay],
                )
                # cloud ice particles
                τi, τi_ssa, τi_ssag = compute_lookup_cld_ice_props(
                    nsize_ice,
                    radice_lwr,
                    radice_upr,
                    lut_extice,
                    lut_ssaice,
                    lut_asyice,
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
                τ[glay], ssa[glay], g[glay] = increment_2stream(τ[glay], ssa[glay], g[glay], τ_cl, ssa_cl, g_cl)
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
        lut_extliq,
        lut_ssaliq,
        lut_asyliq,
        re_liq,
        cld_path_liq,
    )

This function computes the `TwoStream` cloud liquid properties using the `LookUpTable` method.
"""
@inline function compute_lookup_cld_liq_props(
    nsize_liq,
    radliq_lwr,
    radliq_upr,
    lut_extliq,
    lut_ssaliq,
    lut_asyliq,
    re_liq,
    cld_path_liq,
)
    FT = eltype(re_liq)
    τl, τl_ssa, τl_ssag = FT(0), FT(0), FT(0)
    # cloud liquid particles
    if cld_path_liq > eps(FT)
        Δr_liq = (radliq_upr - radliq_lwr) / FT(nsize_liq - 1)
        loc = Int(max(min(unsafe_trunc(Int, (re_liq - radliq_lwr) / Δr_liq) + 1, nsize_liq - 1), 1))
        fac = (re_liq - radliq_lwr - (loc - 1) * Δr_liq) / Δr_liq
        fc1 = FT(1) - fac
        @inbounds begin
            τl = max((fc1 * lut_extliq[loc] + fac * lut_extliq[loc + 1]) * cld_path_liq, FT(0))
            τl_ssa = (fc1 * lut_ssaliq[loc] + fac * lut_ssaliq[loc + 1]) * τl
            τl_ssag = (fc1 * lut_asyliq[loc] + fac * lut_asyliq[loc + 1]) * τl_ssa
        end
    end
    return (τl, τl_ssa, τl_ssag)
end

"""
    compute_lookup_cld_ice_props(
        nsize_ice,
        radice_lwr,
        radice_upr,
        lut_extice,
        lut_ssaice,
        lut_asyice,
        re_ice,
        cld_path_ice,
    )

This function computes the `TwoStream` cloud ice properties using the `LookUpTable` method.
"""
@inline function compute_lookup_cld_ice_props(
    nsize_ice,
    radice_lwr,
    radice_upr,
    lut_extice,
    lut_ssaice,
    lut_asyice,
    re_ice,
    cld_path_ice,
)
    FT = eltype(re_ice)
    τi, τi_ssa, τi_ssag = FT(0), FT(0), FT(0)
    # cloud ice particles
    if cld_path_ice > eps(FT)
        Δr_ice = (radice_upr - radice_lwr) / FT(nsize_ice - 1)
        loc = Int(max(min(unsafe_trunc(Int, (re_ice - radice_lwr) / Δr_ice) + 1, nsize_ice - 1), 1))
        fac = (re_ice - radice_lwr - (loc - 1) * Δr_ice) / Δr_ice
        fc1 = FT(1) - fac
        @inbounds begin
            τi = max((fc1 * lut_extice[loc] + fac * lut_extice[loc + 1]) * cld_path_ice, FT(0))
            τi_ssa = (fc1 * lut_ssaice[loc] + fac * lut_ssaice[loc + 1]) * τi
            τi_ssag = (fc1 * lut_asyice[loc] + fac * lut_asyice[loc + 1]) * τi_ssa
        end
    end
    return (τi, τi_ssa, τi_ssag)
end

"""
    compute_pade_cld_props(
        lkp_cld::PadeCld,
        re_liq::FT,
        re_ice::FT,
        ice_rgh::Int,
        cld_path_liq::FT,
        cld_path_ice::FT,
        ibnd::UInt8,
    ) where {FT}

This function computes the `TwoSteam` cloud optics properties using the pade method.
"""
function compute_pade_cld_props(
    lkp_cld::PadeCld,
    re_liq::FT,
    re_ice::FT,
    ice_rgh::Int,
    cld_path_liq::FT,
    cld_path_ice::FT,
    ibnd::Int,
) where {FT}
    τl, τl_ssa, τl_ssag = FT(0), FT(0), FT(0)
    τi, τi_ssa, τi_ssag = FT(0), FT(0), FT(0)

    (;
        pade_extliq,
        pade_ssaliq,
        pade_asyliq,
        pade_extice,
        pade_ssaice,
        pade_asyice,
        pade_sizreg_extliq,
        pade_sizreg_ssaliq,
        pade_sizreg_asyliq,
        pade_sizreg_extice,
        pade_sizreg_ssaice,
        pade_sizreg_asyice,
    ) = lkp_cld
    m_ext, m_ssa_g = 3, 3
    n_ext, n_ssa_g = 3, 2
    # Finds index into size regime table
    # This works only if there are precisely three size regimes (four bounds) and it's
    # previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
    if cld_path_liq > eps(FT)
        @inbounds begin
            irad = Int(min(floor((re_liq - pade_sizreg_extliq[2]) / pade_sizreg_extliq[3]) + 2, 3))
            τl = pade_eval(ibnd, re_liq, irad, m_ext, n_ext, pade_extliq) * cld_path_liq

            irad = Int(min(floor((re_liq - pade_sizreg_ssaliq[2]) / pade_sizreg_ssaliq[3]) + 2, 3))
            τl_ssa = (FT(1) - max(FT(0), pade_eval(ibnd, re_liq, irad, m_ssa_g, n_ssa_g, pade_ssaliq))) * τl

            irad = Int(min(floor((re_liq - pade_sizreg_asyliq[2]) / pade_sizreg_asyliq[3]) + 2, 3))
            τl_ssag = pade_eval(ibnd, re_liq, irad, m_ssa_g, n_ssa_g, pade_asyliq) * τl_ssa
        end
    end

    if cld_path_ice > eps(FT)
        @inbounds begin
            irad = Int(min(floor((re_ice - pade_sizreg_extice[2]) / pade_sizreg_extice[3]) + 2, 3))

            τi = pade_eval(ibnd, re_ice, irad, m_ext, n_ext, pade_extice, ice_rgh) * cld_path_ice

            irad = Int(min(floor((re_ice - pade_sizreg_ssaice[2]) / pade_sizreg_ssaice[3]) + 2, 3))
            τi_ssa = (FT(1) - max(FT(0), pade_eval(ibnd, re_ice, irad, m_ssa_g, n_ssa_g, pade_ssaice, ice_rgh))) * τi

            irad = Int(min(floor((re_ice - pade_sizreg_asyice[2]) / pade_sizreg_asyice[3]) + 2, 3))
            τi_ssag = pade_eval(ibnd, re_ice, irad, m_ssa_g, n_ssa_g, pade_asyice, ice_rgh) * τi_ssa
        end
    end

    τ = τl + τi
    τ_ssa = τl_ssa + τi_ssa
    τ_ssag = (τl_ssag + τi_ssag) / max(eps(FT), τ_ssa)
    τ_ssa /= max(eps(FT), τ)

    return (τ, τ_ssa, τ_ssag)
end


"""
    pade_eval(
        ibnd,
        re,
        irad,
        m,
        n,
        pade_coeffs,
        irgh::Union{Int,Nothing} = nothing,
    )

Evaluate Pade approximant of order [m/n]
"""
function pade_eval(ibnd, re, irad, m, n, pade_coeffs, irgh::Union{Int, Nothing} = nothing)
    FT = eltype(re)
    if irgh isa Int
        coeffs = view(pade_coeffs, :, :, :, irgh)

    else
        coeffs = pade_coeffs
    end

    @inbounds denom = coeffs[ibnd, irad, n + m]
    i = (n + m - 1)
    @inbounds while i ≥ (1 + m)
        denom = coeffs[ibnd, irad, i] + re * denom
        i -= 1
    end
    denom = FT(1) + re * denom

    @inbounds numer = coeffs[ibnd, irad, m]
    i = (m - 1)
    @inbounds while i ≥ 2
        numer = coeffs[ibnd, irad, i] + re * numer
        i -= 1
    end
    @inbounds numer = coeffs[ibnd, irad, 1] + re * numer

    return (numer / denom)
end

"""
    build_cloud_mask!(cld_mask, cld_frac, ::MaxRandomOverlap)

Builds McICA-sampled cloud mask from cloud fraction data for maximum-random overlap

Reference: https://github.com/AER-RC/RRTMG_SW/
"""
function build_cloud_mask!(
    cld_mask::AbstractArray{Bool, 1},
    cld_frac::AbstractArray{FT, 1},
    ::MaxRandomOverlap,
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
        random_ilayplus1 = Random.rand()
        @inbounds cld_mask[finish] = cld_mask_ilayplus1 = random_ilayplus1 >= (FT(1) - cld_frac_ilayplus1)
        ilay = finish - 1
        while ilay ≥ start
            @inbounds cld_frac_ilay = cld_frac[ilay]
            if cld_frac_ilay > FT(0)
                # use same random number from the layer above if layer above is cloudy
                # update random numbers if layer above is not cloudy
                random_ilay = cld_mask_ilayplus1 ? random_ilayplus1 : Random.rand() * (FT(1) - cld_frac_ilayplus1)
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
