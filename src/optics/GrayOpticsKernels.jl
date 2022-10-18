
# This functions calculates the optical thickness based on pressure 
# and lapse rate for a gray atmosphere. 
# See Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
# https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2
"""
    compute_optical_props_kernel!(
        op::AbstractOpticalProps{FT},
        as::GrayAtmosphericState{FT},
        glaycol,
        source::AbstractSourceLW{FT},
    ) where {FT<:AbstractFloat}

This function computes the optical properties using the gray atmosphere assumption
for the longwave solver.
"""
function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::GrayAtmosphericState{FT},
    glaycol,
    source::AbstractSourceLW{FT},
) where {FT <: AbstractFloat}

    compute_optical_props_kernel_lw!(op, as, glaycol)     # computing optical thickness
    compute_sources_gray_kernel!(source, as, glaycol) # computing Planck sources
end

"""
    τ_lw_gray(p, pꜜ, pꜛ, p₀, τ₀, f)

Optical depth for longwave GrayRadiation.
Reference:
 - TODO: add reference
"""
function τ_lw_gray(p, pꜜ, pꜛ, p₀, τ₀, f)
    FT = eltype(p)
    return τ₀ * (f / p₀ + 4 * (1 - f) / p₀ * (p / p₀)^3) * (pꜜ - pꜛ)
end

function compute_optical_props_kernel_lw!(
    op::AbstractOpticalProps{FT},
    as::GrayAtmosphericState{FT},
    glaycol,
) where {FT <: AbstractFloat}
    glay, gcol = glaycol
    (; p_lay, p_lev, d0) = as
    (; f) = TODO_some_struct
    @inbounds op.τ[glay, gcol] = τ_lw_gray(
        p_lay[glay, gcol],
        p_lev[glay, gcol],
        p_lev[glay + 1, gcol],
        p_lev[1, gcol],
        d0[gcol],
        f, # = FT(0.2)
    )
    if op isa TwoStream
        op.ssa[glaycol...] = FT(0)
        op.g[glaycol...] = FT(0)
    end
end

"""
    τ_sw_gray(p, pꜜ, pꜛ, p₀, τ₀)

Optical depth for shortwave GrayRadiation.
Reference:
 - TODO: add reference
"""
τ_sw_gray(p, pꜜ, pꜛ, p₀, τ₀) = 2 * τ₀ * (p / p₀) / p₀ * (pꜜ - pꜛ)

"""
    compute_optical_props_kernel!(
        op::AbstractOpticalProps{FT},
        as::GrayAtmosphericState{FT},
        glaycol,
    ) where {FT<:AbstractFloat}


This function computes the optical properties using the gray atmosphere assumption
for the shortwave solver.
"""
function compute_optical_props_kernel!(
    op::AbstractOpticalProps{FT},
    as::GrayAtmosphericState{FT},
    glaycol,
) where {FT <: AbstractFloat}
    # setting references
    glay, gcol = glaycol
    (; p_lay, p_lev, d0, α) = as
    @inbounds p0 = p_lev[1, gcol]
    (; f, τ₀) = TODO_some_struct # TODO: where should this be unpacked from ? `Y` vs `Yₜ`
    @inbounds op.τ[glay, gcol] = τ_sw_gray(
        p_lay[glay, gcol],
        p_lev[glay, gcol],
        p_lev[glay + 1, gcol],
        p_lev[1, gcol],
        τ₀, # = FT(0.22) hardcode the value of τ₀ for shortwave gray radiation
    )

    if op isa TwoStream
        op.ssa[glaycol...] = FT(0)
        op.g[glaycol...] = FT(0)
    end
end
"""
    compute_sources_gray_kernel!(
        source::AbstractSourceLW{FT},
        as::GrayAtmosphericState{FT},
        glaycol,
    ) where {FT<:AbstractFloat}

This function computes the Planck sources for the gray longwave solver.
"""
function compute_sources_gray_kernel!(
    source::AbstractSourceLW{FT},
    as::GrayAtmosphericState{FT},
    glaycol,
) where {FT <: AbstractFloat}
    # computing Planck sources
    glay, gcol = glaycol
    (; t_lay, t_lev) = as
    (; lay_source, lev_source_inc, lev_source_dec, sfc_source) = source

    sbc = FT(RP.Stefan(source.param_set))
    @inbounds lay_source[glaycol...] = sbc * t_lay[glaycol...]^FT(4) / FT(π)   # computing lay_source
    @inbounds lev_source_inc[glaycol...] = sbc * t_lev[glay + 1, gcol]^FT(4) / FT(π)
    @inbounds lev_source_dec[glaycol...] = sbc * t_lev[glaycol...]^FT(4) / FT(π)
    if glay == 1
        @inbounds sfc_source[gcol] = sbc * as.t_sfc[gcol]^FT(4) / FT(π)   # computing sfc_source
    end
end
