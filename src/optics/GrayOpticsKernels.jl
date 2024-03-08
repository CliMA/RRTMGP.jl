
# Calculate optical properties for gray radiation solver.

"""
    compute_optical_props!(
        op::AbstractOpticalProps,
        as::GrayAtmosphericState,
        sf::AbstractSourceLW,
        gcol::Int,
    )

Computes optical properties for the longwave gray radiation problem.
"""
function compute_optical_props!(op::OneScalar, as::GrayAtmosphericState, sf::SourceLWNoScat, gcol::Int)
    nlay = AtmosphericStates.get_nlay(as)
    (; p_lay, p_lev, t_lay, t_lev, otp) = as
    (; lay_source, lev_source_inc, lev_source_dec, sfc_source) = sf
    τ = op.τ
    FT = eltype(τ)
    sbc = FT(RP.Stefan(sf.param_set))
    @inbounds begin
        lat = as.lat[gcol]
        p0 = p_lev[1, gcol]
        p_lev_glay = p_lev[1, gcol]
        t_lev_dec = t_lev[1, gcol]
        t_sfc = as.t_sfc[gcol]
        sfc_source[gcol] = sbc * (t_sfc * t_sfc * t_sfc * t_sfc) / FT(π)   # computing sfc_source

        for glay in 1:nlay
            # compute optical thickness
            p_lev_glayplus1 = p_lev[glay + 1, gcol]
            Δp = p_lev_glayplus1 - p_lev_glay
            p = p_lay[glay, gcol]
            τ[glay, gcol] = compute_gray_optical_thickness_lw(otp, p0, Δp, p, lat)
            p_lev_glay = p_lev_glayplus1
            # compute longwave source terms
            t_lev_inc = t_lev[glay + 1, gcol]
            t_lay_loc = t_lay[glay, gcol]
            lay_source[glay, gcol] = sbc * (t_lay_loc * t_lay_loc * t_lay_loc * t_lay_loc) / FT(π)   # computing lay_source
            lev_source_inc[glay, gcol] = sbc * (t_lev_inc * t_lev_inc * t_lev_inc * t_lev_inc) / FT(π)
            lev_source_dec[glay, gcol] = sbc * (t_lev_dec * t_lev_dec * t_lev_dec * t_lev_dec) / FT(π)
            t_lev_dec = t_lev_inc
        end
    end
    return nothing
end

function compute_optical_props!(op::TwoStream, as::GrayAtmosphericState, sf::SourceLW2Str, gcol::Int)
    nlay = AtmosphericStates.get_nlay(as)
    (; p_lay, p_lev, t_lay, t_lev, otp) = as
    (; τ, ssa, g) = op
    (; lev_source, sfc_source) = sf
    FT = eltype(τ)
    sbc = FT(RP.Stefan(sf.param_set))
    @inbounds begin
        lat = as.lat[gcol]
        p0 = p_lev[1, gcol]
        p_lev_glay = p_lev[1, gcol]
        t_lev_dec = t_lev[1, gcol]
        t_sfc = as.t_sfc[gcol]
        sfc_source[gcol] = sbc * (t_sfc * t_sfc * t_sfc * t_sfc) / FT(π)   # computing sfc_source
        lev_src_inc_prev = FT(0)
        lev_src_dec_prev = FT(0)
        for glay in 1:nlay
            p_lev_glayplus1 = p_lev[glay + 1, gcol]
            Δp = p_lev_glayplus1 - p_lev_glay
            p = p_lay[glay, gcol]
            τ[glay, gcol] = compute_gray_optical_thickness_lw(otp, p0, Δp, p, lat)
            p_lev_glay = p_lev_glayplus1
            # compute longwave source terms
            t_lev_inc = t_lev[glay + 1, gcol]
            lev_src_inc = sbc * (t_lev_inc * t_lev_inc * t_lev_inc * t_lev_inc) / FT(π)
            lev_src_dec = sbc * (t_lev_dec * t_lev_dec * t_lev_dec * t_lev_dec) / FT(π)
            if glay == 1
                lev_source[glay, gcol] = lev_src_dec
            else
                lev_source[glay, gcol] = sqrt(lev_src_inc_prev * lev_src_dec)
            end
            lev_src_dec_prev = lev_src_dec
            lev_src_inc_prev = lev_src_inc
            t_lev_dec = t_lev_inc
        end
        lev_source[nlay + 1, gcol] = lev_src_inc_prev
    end
    zeroval = zero(FT)
    map!(x -> zeroval, view(ssa, :, gcol), view(ssa, :, gcol))
    map!(x -> zeroval, view(g, :, gcol), view(g, :, gcol))
    return nothing
end


"""
    compute_optical_props!(
        op::AbstractOpticalProps,
        as::GrayAtmosphericState,
        gcol::Int,
    )

Computes optical properties for the shortwave gray radiation problem.
"""
function compute_optical_props!(op::AbstractOpticalProps, as::GrayAtmosphericState, gcol::Int)
    nlay = AtmosphericStates.get_nlay(as)
    (; p_lay, p_lev, otp) = as
    τ = op.τ
    @inbounds lat = as.lat[gcol]
    @inbounds p0 = p_lev[1, gcol]
    @inbounds p_lev_glay = p_lev[1, gcol]
    for glay in 1:nlay
        @inbounds p_lev_glayplus1 = p_lev[glay + 1, gcol]
        @inbounds Δp = p_lev_glayplus1 - p_lev_glay
        @inbounds p = p_lay[glay, gcol]
        @inbounds τ[glay, gcol] = compute_gray_optical_thickness_sw(otp, p0, Δp, p, lat)
        p_lev_glay = p_lev_glayplus1
    end
    if op isa TwoStream
        (; ssa, g) = op
        FT = eltype(τ)
        zeroval = zero(FT)
        map!(x -> zeroval, view(ssa, :, gcol), view(ssa, :, gcol))
        map!(x -> zeroval, view(g, :, gcol), view(g, :, gcol))
    end
    return nothing
end

"""
    compute_gray_optical_thickness_lw(
        params::GrayOpticalThicknessSchneider2004{FT},
        p0,
        Δp,
        p,
        lat,
    ) where {FT<:AbstractFloat}

This functions calculates the optical thickness based on pressure 
and lapse rate for a gray atmosphere. 
See Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
DOI: https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2
"""
function compute_gray_optical_thickness_lw(params::GrayOpticalThicknessSchneider2004{FT}, p0, Δp, p, lat) where {FT}
    (; α, te, tt, Δt) = params
    # surface temp at a given latitude (K) / temp at top of atmosphere
    ts_by_tt = (te + Δt * (FT(1) / FT(3) - sin(lat / FT(180) * FT(π))^2)) / tt

    ts_by_tt_pow4 = ts_by_tt * ts_by_tt * ts_by_tt * ts_by_tt

    d0 = ts_by_tt_pow4 - FT(1) # optical depth

    return abs((α * d0 * pow_fast(p / p0, α) / p) * Δp)
end

compute_gray_optical_thickness_sw(params::GrayOpticalThicknessSchneider2004{FT}, rest...) where {FT} = FT(0)

"""
    compute_gray_optical_thickness_lw(
        params::GrayOpticalThicknessOGorman2008{FT},
        p0,
        Δp,
        p,
        lat,
    ) where {FT}

This functions calculates the optical thickness based on pressure 
and lapse rate for a gray atmosphere. 
See O'Gorman 2008, Journal of Climate Vol 21, Page(s): 3815–3832.
DOI: https://doi.org/10.1175/2007JCLI2065.1
"""
function compute_gray_optical_thickness_lw(params::GrayOpticalThicknessOGorman2008{FT}, p0, Δp, p, lat) where {FT}
    (; α, fₗ, τₑ, τₚ) = params
    σ = p / p0

    @inbounds τ = (α * Δp / p) * (fₗ * σ + (1 - fₗ) * 4 * σ^4) * (τₑ + (τₚ - τₑ) * sin(lat / FT(180) * FT(π))^2)

    return abs(τ)
end
"""
    compute_gray_optical_thickness_sw(
        params::GrayOpticalThicknessOGorman2008{FT},
        p0,
        Δp,
        p,
        rest...,
    ) where {FT}
 
This functions calculates the optical thickness based on pressure 
for a gray atmosphere. 
See O'Gorman 2008, Journal of Climate Vol 21, Page(s): 3815–3832.
DOI: https://doi.org/10.1175/2007JCLI2065.1
"""
function compute_gray_optical_thickness_sw(params::GrayOpticalThicknessOGorman2008{FT}, p0, Δp, p, rest...) where {FT}
    (; τ₀) = params
    @inbounds τ = 2 * τ₀ * (p / p0) * (Δp / p0)
    return abs(τ)
end
