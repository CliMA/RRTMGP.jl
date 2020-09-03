
# Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
# See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13

@kernel function rte_lw_noscat_gray_source_kernel!(
    lay_source::FTA3D,
    lev_source_up::FTA3D,
    lev_source_dn::FTA3D,
    τ::FTA2D,
    trans::FTA2D,
    source_up::FTA2D,
    source_dn::FTA2D,
    Ds::FTA1D,
    ::Val{nlay},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    nlay,
    ncol,
}
    glay, gcol = @index(Global, NTuple)  # global col & lay ids

    τ_thresh = sqrt(eps(FT)) # or abs(eps(FT))?
    τ_loc = τ[glay, gcol] * Ds[1]      # Optical path and transmission,
    trans[glay, gcol] = exp(-τ_loc)    # used in source function and transport calculations

    # Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
    # is of order epsilon (smallest difference from 1. in working precision)
    # Thanks to Peter Blossey

    fact = (τ_loc > τ_thresh) ?
        ((FT(1) - trans[glay, gcol]) / τ_loc - trans[glay, gcol]) :
        τ_loc * (FT(1 / 2) - FT(1 / 3) * τ_loc)

    # Equations below are developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13

    source_up[glay, gcol] =
        (FT(1) - trans[glay, gcol]) * lev_source_up[glay, gcol, 1] +
        FT(2) *
        fact *
        (lay_source[glay, gcol, 1] - lev_source_up[glay, gcol, 1])

    source_dn[glay, gcol] =
        (FT(1) - trans[glay, gcol]) * lev_source_dn[glay, gcol, 1] +
        FT(2) *
        fact *
        (lay_source[glay, gcol, 1] - lev_source_dn[glay, gcol, 1])

    @synchronize
end

@kernel function rte_lw_noscat_gray_transport_kernel!(
    flux_up::FTA2D,
    flux_dn::FTA2D,
    flux_net::FTA2D,
    trans::FTA2D,
    source_up::FTA2D,
    source_dn::FTA2D,
    sfc_source::FTA2D,
    sfc_emis::FTA1D,
    inc_flux::Union{Nothing,FTA1D},
    τ::FTA2D,
    n_μ::Int,
    Ds::FTA1D,
    w_μ::FTA1D,
    top_at_1::Bool,
    ::Val{nlay},
    ::Val{nlev},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    nlay,
    nlev,
    ncol,
}
    gcol = @index(Global, Linear) # global col ids
    i_lev_top = top_at_1 ? 1 : nlev # index for top level of column

    @inbounds for ilev = 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    if inc_flux ≠ nothing
        flux_dn[i_lev_top, gcol] = inc_flux[gcol]
    end
    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    flux_dn[i_lev_top, gcol] /= (FT(2) * FT(π) * w_μ[1])

    if top_at_1
        #-------Top of domain is index 1
        # Downward propagation
        @inbounds for ilev = 2:nlev
            flux_dn[ilev, gcol] =
                trans[ilev-1, gcol] * flux_dn[ilev-1, gcol] +
                source_dn[ilev-1, gcol]
        end
        # Surface reflection and emission
        flux_up[nlev, gcol] =
            flux_dn[nlev, gcol] * (FT(1) - sfc_emis[gcol]) + sfc_source[gcol, 1]

        # Upward propagation
        @inbounds for ilev = nlay:-1:1
            flux_up[ilev, gcol] =
                trans[ilev, gcol] * flux_up[ilev+1, gcol] +
                source_up[ilev, gcol]
        end
        #-----------------------------------
    else
        #-------Top of domain is index nlev
        # Downward propagation
        @inbounds for ilev = nlay:-1:1
            flux_dn[ilev, gcol] =
                trans[ilev, gcol] * flux_dn[ilev+1, gcol] +
                source_dn[ilev, gcol]
        end

        # Surface reflection and emission
        flux_up[1, gcol] =
            flux_dn[1, gcol] * (FT(1) - sfc_emis[gcol]) + sfc_source[gcol, 1]

        # Upward propagation
        @inbounds for ilev = 2:nlay+1
            flux_up[ilev, gcol] =
                trans[ilev-1, gcol] * flux_up[ilev-1, gcol] +
                source_up[ilev-1, gcol]
        end
        #-----------------------------------
    end

    @inbounds for ilev = 1:nlev
        flux_up[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_dn[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
    @synchronize
end

@kernel function rte_lw_2stream_gray_combine_sources_kernel!(
    lev_source_inc::FTA3D,
    lev_source_dec::FTA3D,
    lev_source::FTA2D,
    ::Val{nlev},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    nlev,
    ncol,
}
    # RRTMGP provides source functions at each level using the spectral mapping
    # of each adjacent layer. Combine these for two-stream calculations
    # lw combine sources
    glev, gcol = @index(Global, NTuple)  # global col & lev ids

    nlay = nlev - 1

    if glev == 1
        lev_source[1, gcol] = lev_source_dec[1, gcol]
    elseif glev == nlev
        lev_source[nlev, gcol] = lev_source_inc[nlay, gcol]
    else
        lev_source[glev, gcol] =
            sqrt(lev_source_dec[glev, gcol] * lev_source_inc[glev-1, gcol])
    end
    @synchronize
end

@kernel function rte_lw_2stream_gray_source_kernel!(
    τ::FTA2D,
    ssa::FTA2D,
    g::FTA2D,
    lev_source::FTA2D,
    src_up::FTA2D,
    src_dn::FTA2D,
    Tdif::FTA2D,
    Rdif::FTA2D,
    lw_diff_sec::FT,
    top_at_1::Bool,
    ::Val{nlay},
    ::Val{ncol},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2},nlay,ncol}
    # Cell properties: reflection, transmission for diffuse radiation
    # Coupling coefficients needed for source function
    # -------------------------------------------------------------------------------------------------
    #
    # Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
    #    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
    #
    # Equations are developed in Meador and Weaver, 1980,
    #    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
    #
    # -------------------------------------------------------------------------------------------------    
    glay, gcol = @index(Global, NTuple)  # global col & layer ids

    γ1 = lw_diff_sec * (1 - FT(0.5) * ssa[glay, gcol] * (1 + g[glay, gcol]))
    γ2 = lw_diff_sec * FT(0.5) * ssa[glay, gcol] * (1 - g[glay, gcol])
    k = sqrt(max((γ1 + γ2) * (γ1 - γ2), FT(1e-12)))

    # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    RT_term =
        1 / (
            k * (1 + exp(-2 * τ[glay, gcol] * k)) +
            γ1 * (1 - exp(-2 * τ[glay, gcol] * k))
        )

    Rdif[glay, gcol] = RT_term * γ2 * (1 - exp(-2 * τ[glay, gcol] * k)) # Equation 25
    Tdif[glay, gcol] = RT_term * 2 * k * exp(-τ[glay, gcol] * k) # Equation 26

    # Source function for diffuse radiation
    # Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
    # This version straight from ECRAD
    # Source is provided as W/m2-str; factor of pi converts to flux units
    # lw_source_2str
    if top_at_1
        top = glay
        bot = glay + 1
    else
        top = glay + 1
        bot = glay
    end

    if τ[glay, gcol] > 1.0e-8
        # Toon et al. (JGR 1989) Eqs 26-27
        Z =
            (lev_source[bot, gcol] - lev_source[top, gcol]) /
            (τ[glay, gcol] * (γ1 + γ2))
        Zup_top = Z + lev_source[top, gcol]
        Zup_bottom = Z + lev_source[bot, gcol]
        Zdn_top = -Z + lev_source[top, gcol]
        Zdn_bottom = -Z + lev_source[bot, gcol]
        src_up[glay, gcol] =
            pi * (
                Zup_top - Rdif[glay, gcol] * Zdn_top -
                Tdif[glay, gcol] * Zup_bottom
            )
        src_dn[glay, gcol] =
            pi * (
                Zdn_bottom - Rdif[glay, gcol] * Zup_bottom -
                Tdif[glay, gcol] * Zdn_top
            )
    else
        src_up[glay, gcol] = FT(0)
        src_dn[glay, gcol] = FT(0)
    end

    @synchronize
end

# ---------------------------------------------------------------
#
# Transport of diffuse radiation through a vertically layered atmosphere.
#   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
#   This routine is shared by longwave and shortwave
#
# ---------------------------------------------------------------
@kernel function adding_gray_kernel!(
    top_at_1::B,
    sfc_emis::FTA1D,
    Rdif::FTA2D,
    Tdif::FTA2D,
    src_up::FTA2D,
    src_dn::FTA2D,
    sfc_source::FTA2D,
    flux_up::FTA2D,
    flux_dn::FTA2D,
    flux_net::FTA2D,
    albedo::FTA2D,
    src::FTA2D,
    ::Val{nlev},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    B<:Bool,
    nlev,
    ncol,
}
    gcol = @index(Global, Linear)  # global col ids

    nlay = nlev - 1

    @inbounds for ilev = 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    if top_at_1
        #--------------------------------------------------------------------------------------
        # Albedo of lowest level is the surface albedo...
        albedo[nlev, gcol] = FT(1) - sfc_emis[gcol]
        # ... and source of diffuse radiation is surface emission
        src[nlev, gcol] = FT(π) * sfc_emis[gcol] * sfc_source[gcol]
        # From bottom to top of atmosphere --
        # compute albedo and source of upward radiation
        @inbounds for ilev = nlay:-1:1
            denom = FT(1) / (FT(1) - Rdif[ilev, gcol] * albedo[ilev+1, gcol])
            albedo[ilev, gcol] =
                Rdif[ilev, gcol] + (
                    Tdif[ilev, gcol] *
                    Tdif[ilev, gcol] *
                    albedo[ilev+1, gcol] *
                    denom
                )
            # Equation 11 -- source is emitted upward radiation at top of layer plus
            # radiation emitted at bottom of layer,
            # transmitted through the layer and reflected from layers below (tdiff*src*albedo)
            src[ilev, gcol] =
                src_up[ilev, gcol] +
                Tdif[ilev, gcol] *
                denom *
                (src[ilev+1, gcol] + albedo[ilev+1, gcol] * src_dn[ilev, gcol])
        end

        # Eq 12, at the top of the domain upwelling diffuse is due to ...
        flux_up[1, gcol] = flux_dn[1, gcol] * albedo[1, gcol] + src[1, gcol] # reflection of incident diffuse and
        # emission from below

        # From the top of the atmosphere downward -- compute fluxes
        @inbounds for ilev = 2:nlay+1
            denom = FT(1) / (FT(1) - Rdif[ilev-1, gcol] * albedo[ilev, gcol])
            flux_dn[ilev, gcol] =
                (
                    Tdif[ilev-1, gcol] * flux_dn[ilev-1, gcol] + # Equation 13
                    Rdif[ilev-1, gcol] * src[ilev, gcol] +
                    src_dn[ilev-1, gcol]
                ) * denom
            flux_up[ilev, gcol] =
                flux_dn[ilev, gcol] * albedo[ilev, gcol] + # Equation 12
                src[ilev, gcol]
        end
        #--------------------------------------------------------------------------------------
    else
        #--------------------------------------------------------------------------------------
        # Albedo of lowest level is the surface albedo...
        albedo[1, gcol] = FT(1) - sfc_emis[gcol]
        # ... and source of diffuse radiation is surface emission
        src[1, gcol] = FT(π) * sfc_emis[gcol] * sfc_source[gcol]
        #--------------------------------------------------------
        # From bottom to top of atmosphere --
        #   compute albedo and source of upward radiation
        @inbounds for ilev = 1:nlay
            denom = FT(1) / (FT(1) - Rdif[ilev, gcol] * albedo[ilev, gcol])  # Eq 10
            albedo[ilev+1, gcol] =
                Rdif[ilev, gcol] +
                Tdif[ilev, gcol] * Tdif[ilev, gcol] * albedo[ilev, gcol] * denom # Equation 9
            # 
            # Equation 11 -- source is emitted upward radiation at top of layer plus
            # radiation emitted at bottom of layer,
            # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
            src[ilev+1, gcol] =
                src_up[ilev, gcol] +
                Tdif[ilev, gcol] *
                denom *
                (src[ilev, gcol] + albedo[ilev, gcol] * src_dn[ilev, gcol])
        end

        # Eq 12, at the top of the domain upwelling diffuse is due to ...
        flux_up[nlev, gcol] =
            flux_dn[nlev, gcol] * albedo[nlev, gcol] + # ... reflection of incident diffuse and
            src[nlev, gcol]                          # scattering by the direct beam below

        # From the top of the atmosphere downward -- compute fluxes
        @inbounds for ilev = nlay:-1:1
            denom = FT(1) / (FT(1) - Rdif[ilev, gcol] * albedo[ilev, gcol])  # Eq 10
            flux_dn[ilev, gcol] =
                (
                    Tdif[ilev, gcol] * flux_dn[ilev+1, gcol] + # Equation 13
                    Rdif[ilev, gcol] * src[ilev, gcol] +
                    src_dn[ilev, gcol]
                ) * denom
            flux_up[ilev, gcol] =
                flux_dn[ilev, gcol] * albedo[ilev, gcol] + # Equation 12
                src[ilev, gcol]
        end
        #--------------------------------------------------------------------------------------
    end

    @inbounds for ilev = 1:nlev
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
    @synchronize
end
