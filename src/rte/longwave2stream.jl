
function rte_lw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    igpt, ibnd = 1, 1
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            compute_optical_props!(op, as, src_lw, gcol)
            rte_lw_2stream!(
                op,
                flux_lw,
                src_lw,
                bcs_lw,
                gcol,
                igpt,
                ibnd,
                nlev,
                ncol,
            )
            compute_net_flux!(flux_lw, gcol)
        end
    end
    return nothing
end

function rte_lw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    band_flux,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    (; cloud_state, aerosol_state) = as
    bld_cld_mask = cloud_state isa CloudState
    flux_up_lw = flux_lw.flux_up
    flux_dn_lw = flux_lw.flux_dn
    (; flux_up, flux_dn) = flux
    FT = eltype(flux_up)
    # zero the optional per-band accumulator (no-op when spectral fluxes are off)
    set_band_flux_to_zero!(band_flux)
    @inbounds begin
        # initialize LW cloud cover accumulator
        if bld_cld_mask && !isnothing(cloud_state.cld_cover_lw)
            cloud_state.cld_cover_lw .= FT(0)
        end
        if aerosol_state isa AerosolState
            ClimaComms.@threaded device for gcol in 1:ncol
                Optics.compute_aero_mask!(
                    view(aerosol_state.aero_mask, :, gcol),
                    view(aerosol_state.aero_mass, :, :, gcol),
                )
            end
        end
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                ibnd = major_gpt2bnd[igpt]
                bld_cld_mask && Optics.build_cloud_mask!(
                    view(cloud_state.mask_lw, :, gcol),
                    view(cloud_state.cld_frac, :, gcol),
                    cloud_state.mask_type,
                )
                # accumulate LW cloud cover
                if bld_cld_mask && !isnothing(cloud_state.cld_cover_lw)
                    cloud_state.cld_cover_lw[gcol] +=
                        any(view(cloud_state.mask_lw, :, gcol)) ? FT(1) : FT(0)
                end
                compute_optical_props!(
                    op,
                    as,
                    src_lw,
                    gcol,
                    igpt,
                    lookup_lw,
                    lookup_lw_cld,
                    lookup_lw_aero,
                )
                rte_lw_2stream!(
                    op,
                    flux,
                    src_lw,
                    bcs_lw,
                    gcol,
                    igpt,
                    ibnd,
                    nlev,
                    ncol,
                )
                if igpt == 1
                    map!(
                        x -> x,
                        view(flux_up_lw, :, gcol),
                        view(flux_up, :, gcol),
                    )
                    map!(
                        x -> x,
                        view(flux_dn_lw, :, gcol),
                        view(flux_dn, :, gcol),
                    )
                else
                    for ilev in 1:nlev
                        @inbounds flux_up_lw[ilev, gcol] += flux_up[ilev, gcol]
                        @inbounds flux_dn_lw[ilev, gcol] += flux_dn[ilev, gcol]
                    end
                end
                # retain this g-point's contribution in its band (no-op when off)
                accumulate_band_flux!(
                    band_flux,
                    flux_up,
                    flux_dn,
                    gcol,
                    ibnd,
                    nlev,
                )
            end
        end
        # normalize LW cloud cover by number of g-points
        if bld_cld_mask && !isnothing(cloud_state.cld_cover_lw)
            ClimaComms.@threaded device for gcol in 1:ncol
                cloud_state.cld_cover_lw[gcol] /= n_gpt
            end
        end
        ClimaComms.@threaded device for gcol in 1:ncol
            compute_net_flux!(flux_lw, gcol)
        end
    end
    return nothing
end

"""
    lw_2stream_coeffs(τ::FT, ssa::FT, g::FT, lev_src_bot::FT, lev_src_top::FT) where {FT}

This function combines RRTMGP-specific sources at levels,
computes layer reflectance, transmittance, and
total source function at levels using linear-in-tau approximation.
"""
@inline function lw_2stream_coeffs(τ, ssa, g, lev_src_bot, lev_src_top)
    FT = eltype(τ)
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
    k_min = Numerics.k_min(FT)
    lw_diff_sec = FT(1.66) # Fu et al. 1997 diffusivity secant
    τ_thresh = 100 * eps(FT)# tau(icol,ilay) > 1.0e-8_wp used in rte-rrtmgp
    # this is chosen to prevent catastrophic cancellation in src_up and src_dn calculation

    γ1 = lw_diff_sec * (1 - FT(0.5) * ssa * (1 + g))
    γ2 = lw_diff_sec * FT(0.5) * ssa * (1 - g)
    # γ1 − γ2 ≡ lw_diff_sec·(1 − ssa) exactly (Fu et al. Eqs 2.9–2.10); the
    # identity avoids the near-1 cancellation at ssa → 1 (see the shortwave
    # counterpart in sw_2stream_coeffs).
    k = sqrt(max(lw_diff_sec * (FT(1) - ssa) * (γ1 + γ2), k_min))

    e1 = exp(-τ * k)
    coeff = e1 * e1 # = exp(-2τk)
    # 1 − e^{−2kτ} via expm1 and the exact factorization (1 − e)(1 + e);
    # the naive 1 − e² loses ~eps/(2kτ) relative accuracy for thin layers.
    one_minus_e2kt = (-expm1(-τ * k)) * (1 + e1)
    # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    RT_term = 1 / (k * (1 + coeff) + γ1 * one_minus_e2kt)

    Rdif = RT_term * γ2 * one_minus_e2kt # Equation 25
    Tdif = RT_term * 2 * k * e1 # Equation 26

    # Source function for diffuse radiation
    # Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
    # This version straight from ECRAD
    # Source is provided as W/m2-str; factor of pi converts to flux units
    # lw_source_2str

    if τ > τ_thresh
        # Toon et al. (JGR 1989) Eqs 26-27
        Z = (lev_src_bot - lev_src_top) / (τ * (γ1 + γ2))
        Zup_top = Z + lev_src_top
        Zup_bottom = Z + lev_src_bot
        Zdn_top = -Z + lev_src_top
        Zdn_bottom = -Z + lev_src_bot
        src_up = pi * (Zup_top - Rdif * Zdn_top - Tdif * Zup_bottom)
        src_dn = pi * (Zdn_bottom - Rdif * Zup_bottom - Tdif * Zdn_top)
    else
        src_up = FT(0)
        src_dn = FT(0)
    end
    return (Rdif, Tdif, src_up, src_dn)
end


"""
    rte_lw_2stream!(
        op::TwoStream,
        flux::FluxLW,
        src_lw::SL,
        bcs_lw::BCL,
        gcol::Int,
        igpt::Int,
        ibnd::Int,
        nlev::Int,
        ncol::Int,
    ) where {SL, BCL}

Two stream solver for the longwave problem.

Transport of diffuse radiation through a vertically layered atmosphere, for the longwave problem.
Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
"""
@inline function rte_lw_2stream!(
    op::TwoStream,
    flux::FluxLW,
    src_lw::SL,
    bcs_lw::BCL,
    gcol::Int,
    igpt::Int,
    ibnd::Int,
    nlev::Int,
    ncol::Int,
) where {SL, BCL}
    nlay = nlev - 1
    # setting references
    (; τ, ssa, g) = op
    (; flux_up, flux_dn, flux_net) = flux

    (; albedo, lev_source, sfc_source, src) = src_lw
    (; inc_flux, sfc_emis) = bcs_lw
    FT = eltype(τ)
    @inbounds flux_dn_ilevplus1 =
        isnothing(inc_flux) ? FT(0) : inc_flux[gcol, igpt]
    @inbounds flux_dn[nlev, gcol] = flux_dn_ilevplus1
    # Albedo of lowest level is the surface albedo...
    @inbounds albedo_ilev = FT(1) - sfc_emis[ibnd, gcol]
    @inbounds albedo[1, gcol] = albedo_ilev
    # ... and source of diffuse radiation is surface emission
    @inbounds src_ilev = FT(π) * sfc_emis[ibnd, gcol] * sfc_source[gcol]
    @inbounds src[1, gcol] = src_ilev#FT(π) * sfc_emis[ibnd, gcol] * sfc_source[gcol]

    # From bottom to top of atmosphere --
    #   compute albedo and source of upward radiation
    @inbounds lev_src_bot = lev_source[1, gcol]
    @inbounds for ilev in 1:nlay
        lev_src_top = lev_source[ilev + 1, gcol]
        τ_lay, ssa_lay, g_lay =
            τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        Rdif, Tdif, src_up, src_dn = lw_2stream_coeffs(
            τ_lay,
            ssa_lay,
            g_lay,
            lev_src_bot,
            lev_src_top,
        )
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        albedo_ilevplus1 = Rdif + Tdif * Tdif * albedo_ilev * denom # Equation 9
        # 
        # Equation 11 -- source is emitted upward radiation at top of layer plus
        # radiation emitted at bottom of layer,
        # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
        src_ilevplus1 =
            src_up + Tdif * denom * (src_ilev + albedo_ilev * src_dn)
        albedo[ilev + 1, gcol], src[ilev + 1, gcol] =
            albedo_ilevplus1, src_ilevplus1
        lev_src_bot = lev_src_top
        albedo_ilev, src_ilev = albedo_ilevplus1, src_ilevplus1
    end

    # Eq 12, at the top of the domain upwelling diffuse is due to ...
    @inbounds flux_up[nlev, gcol] =
        flux_dn_ilevplus1 * albedo[nlev, gcol] + # ... reflection of incident diffuse and
        src[nlev, gcol]                          # scattering by the direct beam below

    # From the top of the atmosphere downward -- compute fluxes
    @inbounds lev_src_top = lev_source[nlay + 1, gcol]
    ilev = nlay
    @inbounds while ilev ≥ 1
        lev_src_bot, albedo_ilev, src_ilev =
            lev_source[ilev, gcol], albedo[ilev, gcol], src[ilev, gcol]
        τ_lay, ssa_lay, g_lay =
            τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        Rdif, Tdif, _, src_dn = lw_2stream_coeffs(
            τ_lay,
            ssa_lay,
            g_lay,
            lev_src_bot,
            lev_src_top,
        )

        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        flux_dn_ilev =
            (Tdif * flux_dn_ilevplus1 + # Equation 13
             Rdif * src_ilev +
             src_dn) * denom
        flux_up[ilev, gcol] =
            flux_dn_ilev * albedo_ilev + # Equation 12
            src_ilev
        flux_dn[ilev, gcol] = flux_dn_ilev
        flux_dn_ilevplus1 = flux_dn_ilev
        lev_src_top = lev_src_bot
        ilev -= 1
    end
end
