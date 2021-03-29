
# Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
# See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13

@kernel function rte_lw_noscat_source_kernel!(
    slv::Solver{FT,I,FTA1D,FTA2D},
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
    @unpack source_up, source_dn, lev_source_inc, lev_source_dec, lay_source =
        slv.src_lw
    Ds = slv.angle_disc.gauss_Ds
    τ = slv.op.τ

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
        (FT(1) - trans) * lev_source_inc[glay, gcol] +
        FT(2) * fact * (lay_source[glay, gcol] - lev_source_inc[glay, gcol])

    source_dn[glay, gcol] =
        (FT(1) - trans) * lev_source_dec[glay, gcol] +
        FT(2) * fact * (lay_source[glay, gcol] - lev_source_dec[glay, gcol])
end

@kernel function rte_lw_noscat_transport_kernel!(
    slv::Solver{FT,I,FTA1D,FTA2D},
    flux::AbstractFlux{FT,FTA2D},
    igpt::Int,
    ibnd::Int,
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
    @unpack source_up, source_dn, sfc_source = slv.src_lw
    @unpack sfc_emis, inc_flux = slv.bcs_lw
    @unpack flux_up, flux_dn, flux_net = flux

    Ds = slv.angle_disc.gauss_Ds
    w_μ = slv.angle_disc.gauss_wts
    n_μ = slv.angle_disc.n_gauss_angles
    τ = slv.op.τ

    @inbounds for ilev = 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    if inc_flux ≠ nothing
        flux_dn[nlev, gcol] = inc_flux[gcol, igpt]
    end
    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    flux_dn[nlev, gcol] /= (FT(2) * FT(π) * w_μ[1])

    # Top of domain is index nlev
    # Downward propagation
    @inbounds for ilev = nlay:-1:1
        trans = exp(-τ[ilev, gcol] * Ds[1])
        flux_dn[ilev, gcol] =
            trans * flux_dn[ilev+1, gcol] + source_dn[ilev, gcol]
    end

    # Surface reflection and emission
    flux_up[1, gcol] =
        flux_dn[1, gcol] * (FT(1) - sfc_emis[ibnd, gcol]) +
        sfc_emis[ibnd, gcol] * sfc_source[gcol]

    # Upward propagation
    @inbounds for ilev = 2:nlay+1
        trans = exp(-τ[ilev-1, gcol] * Ds[1])
        flux_up[ilev, gcol] =
            trans * flux_up[ilev-1, gcol] + source_up[ilev-1, gcol]
    end

    @inbounds for ilev = 1:nlev
        flux_up[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_dn[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
end

@kernel function rte_lw_2stream_combine_sources_kernel!(
    src::SourceLW2Str{FT,FTA1D,FTA2D},
    ::Val{nlev},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    nlev,
    ncol,
}
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
end

@kernel function rte_lw_2stream_source_kernel!(
    slv::Solver{FT,I,FTA1D,FTA2D},
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
    @unpack τ, ssa, g = slv.op
    @unpack Rdif, Tdif, lev_source, src_up, src_dn = slv.src_lw

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
end

# ---------------------------------------------------------------
#
# Transport of diffuse radiation through a vertically layered atmosphere.
#   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
#   This routine is shared by longwave and shortwave
#
# ---------------------------------------------------------------
@kernel function adding_kernel!(
    slv::Solver{FT,I,FTA1D,FTA2D},
    flux::AbstractFlux{FT,FTA2D},
    islw::Bool,
    igpt::Int,
    ibnd::Int,
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
    nlay = nlev - 1
    # setting references
    @unpack flux_up, flux_dn, flux_net = flux

    if islw && slv.src_lw isa SourceLW2Str
        @unpack albedo, sfc_source, Rdif, Tdif, src_up, src_dn, src = slv.src_lw
        sfc_emis = slv.bcs_lw.sfc_emis

        @inbounds for ilev = 1:nlev
            flux_dn[ilev, gcol] = FT(0)
            flux_up[ilev, gcol] = FT(0)
        end

        if slv.bcs_lw.inc_flux ≠ nothing
            flux_dn[nlev, gcol] = slv.bcs_lw.inc_flux[gcol, igpt]
        end
        # Albedo of lowest level is the surface albedo...
        albedo[1, gcol] = FT(1) - sfc_emis[ibnd, gcol]
        # ... and source of diffuse radiation is surface emission
        src[1, gcol] = FT(π) * sfc_emis[ibnd, gcol] * sfc_source[gcol]
    elseif slv.src_sw isa SourceSW2Str
        @unpack albedo, sfc_source, Rdif, Tdif, src_up, src_dn, src = slv.src_sw

        @inbounds for ilev = 1:nlev
            flux_dn[ilev, gcol] = FT(0)
            flux_up[ilev, gcol] = FT(0)
        end

        # Albedo of lowest level is the surface albedo...
        albedo[1, gcol] = slv.bcs_sw.sfc_alb_diffuse[ibnd, gcol]
        # ... and source of diffuse radiation is surface emission
        src[1, gcol] = sfc_source[gcol]
    else # do nothing
    end
    #--------------------------------------------------------------------------------------
    #    # Albedo of lowest level is the surface albedo...
    #    albedo[1, gcol] = FT(1) - sfc_emis[gcol]
    #    # ... and source of diffuse radiation is surface emission
    #    src[1, gcol] = FT(π) * sfc_emis[gcol] * sfc_source[gcol]
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

    @inbounds for ilev = 1:nlev
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
end
# ---------------------------------------------------------------
@kernel function rte_sw_noscat_solve_kernel!(
    slv::Solver{FT,I,FTA1D,FTA2D},
    flux::AbstractFlux{FT,FTA2D},
    solar_frac::FT,
    ::Val{nlev},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    nlev,
    ncol,
}
    gcol = @index(Global, Linear) # global col ids
    @unpack toa_flux, zenith = slv.bcs_sw
    τ = slv.op.τ
    flux_dn_dir = flux.flux_dn_dir

    flux_dn_dir[nlev, gcol] = toa_flux[gcol] * solar_frac * cos(zenith[gcol])

    for ilev = nlev-1:-1:1
        flux_dn_dir[ilev, gcol] =
            flux_dn_dir[ilev+1, gcol] * exp(-τ[ilev, gcol] / cos(zenith[gcol]))
    end
end
#-------------------------------------------------------------------------------------------------
#
# Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
# with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
# 
# Equations are developed in Meador and Weaver, 1980,
# doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
#
#-------------------------------------------------------------------------------------------------
@kernel function sw_two_stream_kernel!(
    slv::Solver{FT,I,FTA1D,FTA2D},
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
    glay, gcol = @index(Global, NTuple)  # global col & lay ids

    zenith = slv.bcs_sw.zenith
    @unpack τ, ssa, g = slv.op
    @unpack Rdif, Tdif, Rdir, Tdir, Tnoscat = slv.src_sw

    # Zdunkowski Practical Improved Flux Method "PIFM"
    #  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    γ1 = (FT(8) - ssa[glay, gcol] * (FT(5) + FT(3) * g[glay, gcol])) * FT(0.25)
    γ2 = FT(3) * (ssa[glay, gcol] * (FT(1) - g[glay, gcol])) * FT(0.25)
    γ3 = (FT(2) - (FT(3) * cos(zenith[gcol])) * g[glay, gcol]) * FT(0.25)
    γ4 = FT(1) - γ3
    α1 = γ1 * γ4 + γ2 * γ3                          # Eq. 16
    α2 = γ1 * γ3 + γ2 * γ4                          # Eq. 17
    k = sqrt(max((γ1 - γ2) * (γ1 + γ2), FT(1e-12)))

    exp_minusktau = exp(-τ[glay, gcol] * k)
    exp_minus2ktau = exp_minusktau * exp_minusktau

    # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    RT_term =
        FT(1) / (k * (FT(1) + exp_minus2ktau) + γ1 * (FT(1) - exp_minus2ktau))

    Rdif[glay, gcol] = RT_term * γ2 * (FT(1) - exp_minus2ktau) # Eqn. 25
    Tdif[glay, gcol] = RT_term * FT(2) * k * exp_minusktau     # Eqn. 26

    # Transmittance of direct, unscattered beam. Also used below
    Tnoscat[glay, gcol] = exp(-τ[glay, gcol] / cos(zenith[gcol]))

    # Direct reflect and transmission
    k_μ = k * cos(zenith[gcol])
    k_γ3 = k * γ3
    k_γ4 = k * γ4

    # Equation 14, multiplying top and bottom by exp(-k*tau)
    #   and rearranging to avoid div by 0.
    if abs(FT(1) - k_μ * k_μ) ≥ eps(FT)
        RT_term = ssa[glay, gcol] * RT_term / (FT(1) - k_μ * k_μ)
    else
        RT_term = ssa[glay, gcol] * RT_term / eps(FT)
    end

    Rdir[glay, gcol] =
        RT_term * (
            (FT(1) - k_μ) * (α2 + k_γ3) -
            (FT(1) + k_μ) * (α2 - k_γ3) * exp_minus2ktau -
            FT(2) * (k_γ3 - α2 * k_μ) * exp_minusktau * Tnoscat[glay, gcol]
        )
    #
    # Equation 15, multiplying top and bottom by exp(-k*tau),
    #   multiplying through by exp(-tau/mu0) to
    #   prefer underflow to overflow
    # Omitting direct transmittance
    #
    Tdir[glay, gcol] =
        -RT_term * (
            (FT(1) + k_μ) * (α1 + k_γ4) * Tnoscat[glay, gcol] -
            (FT(1) - k_μ) * (α1 - k_γ4) * exp_minus2ktau * Tnoscat[glay, gcol] -
            FT(2) * (k_γ4 + α1 * k_μ) * exp_minusktau
        )
end

@kernel function sw_source_2str_kernel!(
    slv::Solver{FT,I,FTA1D,FTA2D},
    flux::AbstractFlux{FT,FTA2D},
    solar_frac::FT,
    ibnd::I,
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
    gcol = @index(Global, Linear) # global col ids

    @unpack toa_flux, zenith, sfc_alb_direct = slv.bcs_sw
    @unpack Rdir, Tdir, Tnoscat, src_up, src_dn, sfc_source = slv.src_sw
    flux_dn_dir = flux.flux_dn_dir

    # layer index = level index
    # previous level is up (+1)
    flux_dn_dir[nlay+1, gcol] = toa_flux[gcol] * solar_frac * cos(zenith[gcol])
    for ilev = nlay:-1:1
        src_up[ilev, gcol] = Rdir[ilev, gcol] * flux_dn_dir[ilev+1, gcol]
        src_dn[ilev, gcol] = Tdir[ilev, gcol] * flux_dn_dir[ilev+1, gcol]
        flux_dn_dir[ilev, gcol] =
            Tnoscat[ilev, gcol] * flux_dn_dir[ilev+1, gcol]
    end
    sfc_source[gcol] = flux_dn_dir[1, gcol] * sfc_alb_direct[ibnd, gcol]
end
