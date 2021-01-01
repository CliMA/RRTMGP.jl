
# Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
# See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13

@kernel function rte_lw_noscat_gray_source_kernel!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D},
    ::Val{nlay},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    nlay,
    ncol,
}
    glay, gcol = @index(Global, NTuple)  # global col & lay ids

    # setting references
    source_up = gray_rrtmgp.src.source_up
    source_dn = gray_rrtmgp.src.source_dn
    lev_source_up = gray_rrtmgp.src.lev_source_inc
    lev_source_dn = gray_rrtmgp.src.lev_source_dec
    lay_source = gray_rrtmgp.src.lay_source
    Ds = gray_rrtmgp.angle_disc.gauss_Ds
    τ = gray_rrtmgp.op.τ
    #---------------------

    τ_thresh = sqrt(eps(FT)) # or abs(eps(FT))?
    τ_loc = τ[glay, gcol] * Ds[1]      # Optical path and transmission,

    trans = exp(-τ_loc)    # used in source function and transport calculations
    # Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
    # is of order epsilon (smallest difference from 1. in working precision)
    # Thanks to Peter Blossey

    fact = (τ_loc > τ_thresh) ? ((FT(1) - trans) / τ_loc - trans) :
        τ_loc * (FT(1 / 2) - FT(1 / 3) * τ_loc)

    # Equations below are developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13

    source_up[glay, gcol] =
        (FT(1) - trans) * lev_source_up[glay, gcol] +
        FT(2) * fact * (lay_source[glay, gcol] - lev_source_up[glay, gcol])

    source_dn[glay, gcol] =
        (FT(1) - trans) * lev_source_dn[glay, gcol] +
        FT(2) * fact * (lay_source[glay, gcol] - lev_source_dn[glay, gcol])

    @synchronize
end


@kernel function rte_lw_noscat_gray_transport_kernel!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D},
    ::Val{nlay},
    ::Val{nlev},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    nlay,
    nlev,
    ncol,
}
    gcol = @index(Global, Linear) # global col ids
    # setting references
    source_up = gray_rrtmgp.src.source_up
    source_dn = gray_rrtmgp.src.source_dn
    sfc_source = gray_rrtmgp.src.sfc_source
    flux_up = gray_rrtmgp.flux.flux_up
    flux_dn = gray_rrtmgp.flux.flux_dn
    flux_net = gray_rrtmgp.flux.flux_net
    sfc_emis = gray_rrtmgp.bcs.sfc_emis
    inc_flux = gray_rrtmgp.bcs.inc_flux
    Ds = gray_rrtmgp.angle_disc.gauss_Ds
    w_μ = gray_rrtmgp.angle_disc.gauss_wts
    n_μ = gray_rrtmgp.angle_disc.n_gauss_angles
    τ = gray_rrtmgp.op.τ
    #--------------------------------------------
    @inbounds for ilev = 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    if inc_flux ≠ nothing
        flux_dn[nlev, gcol] = inc_flux[gcol]
    end
    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    flux_dn[nlev, gcol] /= (FT(2) * FT(π) * w_μ[1])

    #-------Top of domain is index nlev
    # Downward propagation
    @inbounds for ilev = nlay:-1:1
        trans = exp(-τ[ilev, gcol] * Ds[1])
        flux_dn[ilev, gcol] =
            trans * flux_dn[ilev+1, gcol] + source_dn[ilev, gcol]
    end

    # Surface reflection and emission
    flux_up[1, gcol] =
        flux_dn[1, gcol] * (FT(1) - sfc_emis[gcol]) + sfc_source[gcol, 1]

    # Upward propagation
    @inbounds for ilev = 2:nlay+1
        trans = exp(-τ[ilev-1, gcol] * Ds[1])
        flux_up[ilev, gcol] =
            trans * flux_up[ilev-1, gcol] + source_up[ilev-1, gcol]
    end
    #-----------------------------------
    @inbounds for ilev = 1:nlev
        flux_up[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_dn[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
    @synchronize
end

@kernel function rte_lw_2stream_gray_combine_sources_kernel!(
    src::GraySourceLW2Str{FT,FTA2D},
    ::Val{nlev},
    ::Val{ncol},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2},nlev,ncol}
    # RRTMGP provides source functions at each level using the spectral mapping
    # of each adjacent layer. Combine these for two-stream calculations
    # lw combine sources
    glev, gcol = @index(Global, NTuple)  # global col & lev ids

    nlay = nlev - 1

    if glev == 1
        src.lev_source[1, gcol] = src.lev_source_dec[1, gcol]
    elseif glev == nlev
        src.lev_source[nlev, gcol] = src.lev_source_inc[nlay, gcol]
    else
        src.lev_source[glev, gcol] = sqrt(
            src.lev_source_dec[glev, gcol] * src.lev_source_inc[glev-1, gcol],
        )
    end
    @synchronize
end

@kernel function rte_lw_2stream_gray_source_kernel!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D},
    ::Val{nlay},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
    nlay,
    ncol,
}
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
    # setting references
    τ = gray_rrtmgp.op.τ
    ssa = gray_rrtmgp.op.ssa
    g = gray_rrtmgp.op.g
    Rdif = gray_rrtmgp.src.Rdif
    Tdif = gray_rrtmgp.src.Tdif
    lev_source = gray_rrtmgp.src.lev_source
    src_up = gray_rrtmgp.src.src_up
    src_dn = gray_rrtmgp.src.src_dn
    # ------------------
    lw_diff_sec = FT(1.66)

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

    top = glay + 1
    bot = glay

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
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D},
    ::Val{nlev},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
    nlev,
    ncol,
}
    gcol = @index(Global, Linear)  # global col ids
    # setting references
    albedo = gray_rrtmgp.src.albedo
    sfc_source = gray_rrtmgp.src.sfc_source
    Rdif = gray_rrtmgp.src.Rdif
    Tdif = gray_rrtmgp.src.Tdif
    src_up = gray_rrtmgp.src.src_up
    src_dn = gray_rrtmgp.src.src_dn
    src = gray_rrtmgp.src.src
    flux_up = gray_rrtmgp.flux.flux_up
    flux_dn = gray_rrtmgp.flux.flux_dn
    flux_net = gray_rrtmgp.flux.flux_net
    sfc_emis = gray_rrtmgp.bcs.sfc_emis
    # ------------------
    nlay = nlev - 1

    @inbounds for ilev = 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

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

    @inbounds for ilev = 1:nlev
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
    @synchronize
end
