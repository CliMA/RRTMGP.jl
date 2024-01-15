
# Calculate optical properties for gray radiation solver.

"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::GrayAtmosphericState{FT},
        sf::AbstractSourceLW{FT},
        gcol::Int,
    )

Computes optical properties for the gray radiation problem.
"""
@inline function compute_optical_props!(
    grayoptfun::Function,
    op::AbstractOpticalProps,
    as::GrayAtmosphericState,
    gcol::Int,
)
    (; p_lay, p_lev, otp, nlay) = as
    τ = op.τ
    @inbounds lat = as.lat[gcol]
    @inbounds p0 = p_lev[1, gcol]
    @inbounds p_lev_glay = p_lev[1, gcol]
    for glay in 1:nlay
        @inbounds p_lev_glayplus1 = p_lev[glay + 1, gcol]
        @inbounds Δp = p_lev_glayplus1 - p_lev_glay
        @inbounds p = p_lay[glay, gcol]
        @inbounds τ[glay, gcol] = grayoptfun(otp, p0, Δp, p, lat)
        p_lev_glay = p_lev_glayplus1
    end
    if op isa TwoStream
        (; ssa, g) = op
        FT = eltype(ssa)
        zeroval = zero(FT)
        map!(x -> zeroval, view(ssa, :, gcol), view(ssa, :, gcol))
        map!(x -> zeroval, view(g, :, gcol), view(g, :, gcol))
    end
    return nothing
end

"""
    compute_sources_gray_kernel!(
        as::GrayAtmosphericState{FT},
        sf::AbstractSourceLW{FT},
        glay, gcol,
    ) where {FT<:AbstractFloat}

This function computes the Planck sources for the gray longwave solver.
"""
@inline function compute_longwave_sources!(as::GrayAtmosphericState, sf::SourceLWNoScat, gcol::Int)
    # computing Planck sources
    (; t_lay, t_lev, t_sfc, nlay) = as
    (; lay_source, lev_source_inc, lev_source_dec, sfc_source) = sf
    FT = eltype(t_lay)
    sbc = FT(RP.Stefan(sf.param_set))
    for glay in 1:nlay
        @inbounds lay_source[glay, gcol] = sbc * t_lay[glay, gcol]^FT(4) / FT(π)   # computing lay_source
        @inbounds lev_source_inc[glay, gcol] = sbc * t_lev[glay + 1, gcol]^FT(4) / FT(π)
        @inbounds lev_source_dec[glay, gcol] = sbc * t_lev[glay, gcol]^FT(4) / FT(π)
        if glay == 1
            @inbounds sfc_source[gcol] = sbc * t_sfc[gcol]^FT(4) / FT(π)   # computing sfc_source
        end
    end
    return nothing
end

@inline function compute_longwave_sources!(as::GrayAtmosphericState, sf::SourceLW2Str, gcol::Int)
    # computing Planck sources
    (; t_lay, t_lev, t_sfc, nlay) = as
    #(; lay_source, lev_source, sfc_source) = sf
    (; lev_source, sfc_source) = sf
    FT = eltype(t_lay)
    sbc = FT(RP.Stefan(sf.param_set))
    lev_source_inc_prev = FT(0)
    lev_source_dec_prev = FT(0)
    @inbounds begin
        for glay in 1:nlay
            #lay_source[glay, gcol] = sbc * t_lay[glay, gcol]^FT(4) / FT(π)   # computing lay_source
            lev_source_inc = sbc * t_lev[glay + 1, gcol]^FT(4) / FT(π)
            lev_source_dec = sbc * t_lev[glay, gcol]^FT(4) / FT(π)
            if glay == 1
                sfc_source[gcol] = sbc * t_sfc[gcol]^FT(4) / FT(π)   # computing sfc_source
                lev_source[glay, gcol] = lev_source_dec
            else
                lev_source[glay, gcol] = sqrt(lev_source_dec * lev_source_inc_prev)
            end
            lev_source_inc_prev = lev_source_inc
            lev_source_dec_prev = lev_source_dec
        end
        lev_source[nlay + 1, gcol] = lev_source_inc_prev
    end
    return nothing
end

"""
    gray_optical_thickness_lw(
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
@inline function gray_optical_thickness_lw(
    params::GrayOpticalThicknessSchneider2004{FT},
    p0,
    Δp,
    p,
    lat,
) where {FT <: AbstractFloat}
    (; α, te, tt, Δt) = params
    ts = te + Δt * (FT(1) / FT(3) - sin(lat / FT(180) * FT(π))^2) # surface temp at a given latitude (K)
    d0 = FT((ts / tt)^FT(4) - FT(1)) # optical depth
    @inbounds τ = (α * d0 * pow_fast(p / p0, α) / p) * Δp
    return abs(τ)
end

@inline gray_optical_thickness_sw(params::GrayOpticalThicknessSchneider2004{FT}, rest...) where {FT} = FT(0)

"""
    gray_optical_thickness_lw(
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
@inline function gray_optical_thickness_lw(params::GrayOpticalThicknessOGorman2008{FT}, p0, Δp, p, lat) where {FT}
    (; α, fₗ, τₑ, τₚ) = params
    σ = p / p0

    @inbounds τ = (α * Δp / p) * (fₗ * σ + (1 - fₗ) * 4 * σ^4) * (τₑ + (τₚ - τₑ) * sin(lat / FT(180) * FT(π))^2)

    return abs(τ)
end
"""
    gray_optical_thickness_sw(
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
@inline function gray_optical_thickness_sw(params::GrayOpticalThicknessOGorman2008{FT}, p0, Δp, p, rest...) where {FT}
    (; τ₀) = params
    @inbounds τ = 2 * τ₀ * (p / p0) * (Δp / p0)
    return abs(τ)
end
