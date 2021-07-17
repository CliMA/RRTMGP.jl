
"""
    rte_lw_noscat_source!(
        src_lw::SourceLWNoScat{FT},
        op::OneScalar{FT},
        gcol,
        nlay,
    ) where {FT<:AbstractFloat}

Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
"""
function rte_lw_noscat_source!(
    src_lw::SourceLWNoScat{FT},
    op::OneScalar{FT},
    gcol,
    nlay,
) where {FT<:AbstractFloat}
    # setting references
    @unpack src_up, src_dn, lev_source_inc, lev_source_dec, lay_source = src_lw
    #    @unpack nlay = slv.as
    Ds = op.angle_disc.gauss_Ds
    τ = op.τ

    τ_thresh = sqrt(eps(FT)) # or abs(eps(FT))?

    @inbounds for glay = 1:nlay
        τ_loc = τ[glay, gcol] * Ds[1]      # Optical path and transmission,

        trans = exp(-τ_loc)    # used in source function and transport calculations
        # Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
        # is of order epsilon (smallest difference from 1. in working precision)
        # Thanks to Peter Blossey

        fact = (τ_loc > τ_thresh) ? ((FT(1) - trans) / τ_loc - trans) :
            τ_loc * (FT(1 / 2) - FT(1 / 3) * τ_loc)

        # Equations below are developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13

        src_up[glay, gcol] =
            (FT(1) - trans) * lev_source_inc[glay, gcol] +
            FT(2) * fact * (lay_source[glay, gcol] - lev_source_inc[glay, gcol])

        src_dn[glay, gcol] =
            (FT(1) - trans) * lev_source_dec[glay, gcol] +
            FT(2) * fact * (lay_source[glay, gcol] - lev_source_dec[glay, gcol])
    end
end

"""
    rte_lw_noscat_transport!(
        src_lw::SourceLWNoScat{FT},
        bcs_lw::LwBCs{FT},
        op::OneScalar{FT},
        gcol,
        flux::FluxLW{FT},
        igpt,
        ibnd,
        nlay,
        nlev,
    ) where {FT<:AbstractFloat}

Transport for no-scattering longwave problem.
"""
function rte_lw_noscat_transport!(
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
    gcol,
    flux::FluxLW{FT},
    igpt,
    ibnd,
    nlay,
    nlev,
) where {FT<:AbstractFloat}
    # setting references
    @unpack src_up, src_dn, sfc_source = src_lw
    @unpack sfc_emis, inc_flux = bcs_lw
    @unpack flux_up, flux_dn, flux_net = flux

    Ds = op.angle_disc.gauss_Ds
    w_μ = op.angle_disc.gauss_wts
    n_μ = op.angle_disc.n_gauss_angles
    τ = op.τ

    @inbounds for ilev = 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    if inc_flux ≠ nothing
        @inbounds flux_dn[nlev, gcol] = inc_flux[gcol, igpt]
    end
    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    @inbounds flux_dn[nlev, gcol] /= (FT(2) * FT(π) * w_μ[1])

    # Top of domain is index nlev
    # Downward propagation
    @inbounds for ilev = nlay:-1:1
        trans = exp(-τ[ilev, gcol] * Ds[1])
        flux_dn[ilev, gcol] = trans * flux_dn[ilev+1, gcol] + src_dn[ilev, gcol]
    end

    # Surface reflection and emission
    @inbounds flux_up[1, gcol] =
        flux_dn[1, gcol] * (FT(1) - sfc_emis[ibnd, gcol]) +
        sfc_emis[ibnd, gcol] * sfc_source[gcol]

    # Upward propagation
    @inbounds for ilev = 2:nlay+1
        trans = exp(-τ[ilev-1, gcol] * Ds[1])
        flux_up[ilev, gcol] =
            trans * flux_up[ilev-1, gcol] + src_up[ilev-1, gcol]
    end

    @inbounds for ilev = 1:nlev
        flux_up[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_dn[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
end

"""
    rte_lw_2stream_combine_sources!(src::SourceLW2Str, gcol, nlev, ncol)

RRTMGP provides source functions at each level using the spectral mapping
of each adjacent layer. This function combines these for two-stream calculations.
"""
function rte_lw_2stream_combine_sources!(src::SourceLW2Str, gcol, nlev, ncol)
    @inbounds src.lev_source[1, gcol] = src.lev_source_dec[1, gcol] # glev == 1
    @inbounds src.lev_source[nlev, gcol] = src.lev_source_inc[nlev-1, gcol] # glev == nlev
    @inbounds for glev = 2:nlev-1
        src.lev_source[glev, gcol] = sqrt(
            src.lev_source_dec[glev, gcol] * src.lev_source_inc[glev-1, gcol],
        )
    end
end

"""
    rte_lw_2stream_source!(
        op::TwoStream{FT},
        src_lw::SourceLW2Str{FT},
        gcol,
        nlay,
        ncol,
    ) where {FT<:AbstractFloat}

This function combines RRTMGP-specific sources at levels,
computes layer reflectance, transmittance, and
total source function at levels using linear-in-tau approximation.
"""
function rte_lw_2stream_source!(
    op::TwoStream{FT},
    src_lw::SourceLW2Str{FT},
    gcol,
    nlay,
    ncol,
) where {FT<:AbstractFloat}
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
    # setting references
    @unpack τ, ssa, g = op
    @unpack Rdif, Tdif, lev_source, src_up, src_dn = src_lw

    lw_diff_sec = FT(1.66)

    @inbounds for glay = 1:nlay
        γ1 = lw_diff_sec * (1 - FT(0.5) * ssa[glay, gcol] * (1 + g[glay, gcol]))
        γ2 = lw_diff_sec * FT(0.5) * ssa[glay, gcol] * (1 - g[glay, gcol])
        k = sqrt(max((γ1 + γ2) * (γ1 - γ2), FT(1e-12)))

        # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
        RT_term =
            1 / (
                k * (1 + exp(-2 * τ[glay, gcol] * k)) +
                γ1 * (1 - exp(-2 * τ[glay, gcol] * k))
            )

        @inbounds Rdif[glay, gcol] =
            RT_term * γ2 * (1 - exp(-2 * τ[glay, gcol] * k)) # Equation 25
        @inbounds Tdif[glay, gcol] = RT_term * 2 * k * exp(-τ[glay, gcol] * k) # Equation 26

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
            @inbounds src_up[glay, gcol] =
                pi * (
                    Zup_top - Rdif[glay, gcol] * Zdn_top -
                    Tdif[glay, gcol] * Zup_bottom
                )
            @inbounds src_dn[glay, gcol] =
                pi * (
                    Zdn_bottom - Rdif[glay, gcol] * Zup_bottom -
                    Tdif[glay, gcol] * Zdn_top
                )
        else
            @inbounds src_up[glay, gcol] = FT(0)
            @inbounds src_dn[glay, gcol] = FT(0)
        end
    end
end

"""
    adding_lw!(
        flux::FluxLW{FT},
        src_lw::SL,
        bcs_lw::BCL,
        gcol,
        igpt,
        ibnd,
        nlev,
        ncol,
    ) where {FT<:AbstractFloat,SL<:SourceLW2Str{FT},BCL<:LwBCs{FT}}

Transport of diffuse radiation through a vertically layered atmosphere, for the longwave problem.
Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
"""
function adding_lw!(
    flux::FluxLW{FT},
    src_lw::SL,
    bcs_lw::BCL,
    gcol,
    igpt,
    ibnd,
    nlev,
    ncol,
) where {FT<:AbstractFloat,SL<:SourceLW2Str{FT},BCL<:LwBCs{FT}}
    nlay = nlev - 1
    # setting references
    @unpack flux_up, flux_dn, flux_net = flux

    @unpack albedo, sfc_source, Rdif, Tdif, src_up, src_dn, src = src_lw
    @unpack inc_flux, sfc_emis = bcs_lw

    @inbounds for ilev = 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    if inc_flux ≠ nothing
        @inbounds flux_dn[nlev, gcol] = inc_flux[gcol, igpt]
    end
    # Albedo of lowest level is the surface albedo...
    @inbounds albedo[1, gcol] = FT(1) - sfc_emis[ibnd, gcol]
    # ... and source of diffuse radiation is surface emission
    @inbounds src[1, gcol] = FT(π) * sfc_emis[ibnd, gcol] * sfc_source[gcol]

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
    @inbounds flux_up[nlev, gcol] =
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
"""
    rte_sw_noscat_solve_kernel!(
        flux::FluxSW{FT},
        op::OneScalar{FT},
        bcs_sw::SwBCs{FT},
        solar_frac::FT,
        gcol,
        nlev,
    ) where {FT<:AbstractFloat}

No-scattering solver for the shortwave problem.
(Extinction-only i.e. solar direct beam)
"""
function rte_sw_noscat_solve_kernel!(
    flux::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    solar_frac::FT,
    gcol,
    nlev,
) where {FT<:AbstractFloat}
    @unpack toa_flux, zenith = bcs_sw
    τ = op.τ
    flux_dn_dir = flux.flux_dn_dir
    # downward propagation
    flux_dn_dir[nlev, gcol] = toa_flux[gcol] * solar_frac * cos(zenith[gcol])
    for ilev = nlev-1:-1:1
        flux_dn_dir[ilev, gcol] =
            flux_dn_dir[ilev+1, gcol] * exp(-τ[ilev, gcol] / cos(zenith[gcol]))
    end
end

"""
    sw_two_stream!(
        op::TwoStream{FT},
        src_sw::SourceSW2Str{FT},
        bcs_sw::SwBCs{FT},
        gcol,
        nlay,
    ) where {FT<:AbstractFloat}

Computes cell properties (transmittance and reflectance) for direct and diffuse radiation

Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
 
Equations are developed in Meador and Weaver, 1980,
doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
"""
function sw_two_stream!(
    op::TwoStream{FT},
    src_sw::SourceSW2Str{FT},
    bcs_sw::SwBCs{FT},
    gcol,
    nlay,
) where {FT<:AbstractFloat}
    zenith = bcs_sw.zenith
    @unpack τ, ssa, g = op
    @unpack Rdif, Tdif, Rdir, Tdir, Tnoscat = src_sw

    @inbounds for glay = 1:nlay
        # Zdunkowski Practical Improved Flux Method "PIFM"
        #  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
        γ1 =
            (FT(8) - ssa[glay, gcol] * (FT(5) + FT(3) * g[glay, gcol])) *
            FT(0.25)
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
            FT(1) /
            (k * (FT(1) + exp_minus2ktau) + γ1 * (FT(1) - exp_minus2ktau))

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
                (FT(1) - k_μ) *
                (α1 - k_γ4) *
                exp_minus2ktau *
                Tnoscat[glay, gcol] - FT(2) * (k_γ4 + α1 * k_μ) * exp_minusktau
            )
    end
end

"""
    sw_source_2str!(
        src_sw::SourceSW2Str{FT},
        bcs_sw::SwBCs{FT},
        gcol,
        flux::FluxSW{FT},
        solar_frac::FT,
        ibnd,
        nlay,
    ) where {FT<:AbstractFloat}

Direct-beam and source for diffuse radiation
"""
function sw_source_2str!(
    src_sw::SourceSW2Str{FT},
    bcs_sw::SwBCs{FT},
    gcol,
    flux::FluxSW{FT},
    solar_frac::FT,
    ibnd,
    nlay,
) where {FT<:AbstractFloat}
    @unpack toa_flux, zenith, sfc_alb_direct = bcs_sw
    @unpack Rdir, Tdir, Tnoscat, src_up, src_dn, sfc_source = src_sw
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

"""
    adding_sw!(
        src_sw::SourceSW2Str{FT},
        bcs_sw::SwBCs{FT},
        gcol,
        flux::FluxSW{FT},
        igpt,
        ibnd,
        nlev,
    ) where {FT<:AbstractFloat}

Transport for the two stream shortwave problem.

Transport of diffuse radiation through a vertically layered atmosphere.
Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
"""
function adding_sw!(
    src_sw::SourceSW2Str{FT},
    bcs_sw::SwBCs{FT},
    gcol,
    flux::FluxSW{FT},
    igpt,
    ibnd,
    nlev,
) where {FT<:AbstractFloat}
    nlay = nlev - 1
    # setting references
    @unpack flux_up, flux_dn, flux_net = flux


    @unpack albedo, sfc_source, Rdif, Tdif, src_up, src_dn, src = src_sw

    @inbounds for ilev = 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    # Albedo of lowest level is the surface albedo...
    @inbounds albedo[1, gcol] = bcs_sw.sfc_alb_diffuse[ibnd, gcol]
    # ... and source of diffuse radiation is surface emission
    @inbounds src[1, gcol] = sfc_source[gcol]
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
