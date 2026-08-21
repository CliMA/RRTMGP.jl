
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
            compute_net_flux!(flux_lw, gcol, nlev)
        end
    end
    return nothing
end

# Device-agnostic per-(g-point, column) body, shared by the CPU driver below
# and the CUDA kernel in ext/cuda. Returns whether this g-point had any cloudy
# layer (for the cloud-cover diagnostic).
@inline function lw_2stream_gpt_col!(
    igpt,
    gcol,
    flux,
    flux_lw,
    band_flux,
    src_lw,
    bcs_lw,
    op,
    as,
    state_cache,
    lookup_lw,
    lookup_lw_cld,
    lookup_lw_aero,
    ibnd,
    nlev,
    ncol,
)
    cloudy = _build_cloud_mask!(as.cloud_state, Val(:mask_lw), gcol)
    compute_optical_props!(
        op,
        as,
        state_cache,
        src_lw,
        gcol,
        igpt,
        lookup_lw,
        lookup_lw_cld,
        lookup_lw_aero,
    )
    rte_lw_2stream!(op, flux, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
    _accumulate_fluxes!(flux_lw, flux, gcol, nlev, igpt)
    # retain this g-point's contribution in its band (no-op when off)
    accumulate_band_flux!(band_flux, flux.flux_up, flux.flux_dn, gcol, ibnd, nlev)
    return cloudy
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
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    (; cloud_state, aerosol_state) = as
    track_cld_cover =
        cloud_state isa CloudState && !isnothing(cloud_state.cld_cover_lw)
    FT = eltype(flux_lw.flux_up)
    # zero the optional per-band accumulator (no-op when spectral fluxes are off)
    set_band_flux_to_zero!(band_flux)
    @inbounds begin
        track_cld_cover && (cloud_state.cld_cover_lw .= FT(0))
        if aerosol_state isa AerosolState
            ClimaComms.@threaded device for gcol in 1:ncol
                _compute_aero_mask!(aerosol_state, gcol)
            end
        end
        for igpt in 1:n_gpt
            ibnd = major_gpt2bnd[igpt]
            ClimaComms.@threaded device for gcol in 1:ncol
                cloudy = lw_2stream_gpt_col!(
                    igpt,
                    gcol,
                    flux,
                    flux_lw,
                    band_flux,
                    src_lw,
                    bcs_lw,
                    op,
                    as,
                    state_cache,
                    lookup_lw,
                    lookup_lw_cld,
                    lookup_lw_aero,
                    ibnd,
                    nlev,
                    ncol,
                )
                track_cld_cover &&
                    (cloud_state.cld_cover_lw[gcol] += FT(cloudy))
            end
        end
        # normalize LW cloud cover by number of g-points
        if track_cld_cover
            ClimaComms.@threaded device for gcol in 1:ncol
                cloud_state.cld_cover_lw[gcol] /= n_gpt
            end
        end
        ClimaComms.@threaded device for gcol in 1:ncol
            compute_net_flux!(flux_lw, gcol, nlev)
        end
    end
    return nothing
end

"""
    lw_2stream_coeffs(τ::FT, ssa::FT, g::FT, lev_src_bot::FT, lev_src_top::FT) where {FT}

Combine RRTMGP-specific sources at levels and compute the layer reflectance,
transmittance, and total source function at levels, using the linear-in-tau
approximation.
"""
@inline function lw_2stream_coeffs(τ, ssa, g, lev_src_bot, lev_src_top)
    FT = eltype(τ)
    # Cell properties: reflection, transmission for diffuse radiation
    # Coupling coefficients needed for source function
    # -------------------------------------------------------------------------------------------------
    #
    # Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
    #    with optical depth tau, single scattering albedo w0, and asymmetry parameter g.
    #
    # Equations are developed in Meador and Weaver, 1980,
    #    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
    #
    # -------------------------------------------------------------------------------------------------    
    # setting references
    k_min = Numerics.k_min(FT)
    lw_diff_sec = FT(1.66) # Fu et al. 1997 diffusivity secant

    γ1 = lw_diff_sec * (1 - FT(0.5) * ssa * (1 + g))
    γ2 = lw_diff_sec * FT(0.5) * ssa * (1 - g)
    # γ1 − γ2 ≡ lw_diff_sec·(1 − ssa) exactly (Fu et al. Eqs 2.9–2.10); the
    # identity avoids the near-1 cancellation at ssa → 1 (see the shortwave
    # counterpart in sw_2stream_coeffs).
    k = sqrt(max(lw_diff_sec * (FT(1) - ssa) * (γ1 + γ2), k_min))

    e1 = exp(-τ * k)
    om1 = -expm1(-τ * k) # = 1 − e1, accurate for small kτ
    coeff = e1 * e1 # = exp(-2τk)
    # 1 − e^{−2kτ} via the exact factorization (1 − e)(1 + e); the naive
    # 1 − e² loses ~eps/(2kτ) relative accuracy for thin layers.
    one_minus_e2kt = om1 * (1 + e1)
    # Formulated to avoid rounding errors when k, gamma1 are of very different magnitudes
    RT_term = 1 / (k * (1 + coeff) + γ1 * one_minus_e2kt)

    Rdif = RT_term * γ2 * one_minus_e2kt # Equation 25
    Tdif = RT_term * 2 * k * e1 # Equation 26

    # Source function for diffuse radiation: linear-in-τ (Toon et al. 1989,
    # JGR, Eqs 26-27, as in ecRad's lw_source_2str), formulated so nothing
    # divides by τ alone. Sources are W/m²-str; π converts to flux units.
    #
    # With B_t/B_b the level sources and Z = (B_b − B_t)/(τ(γ1+γ2)), the
    # Toon sources are algebraically
    #   src_up/π = B_t·(1 − Rdif − Tdif) + Tdif·(B_t − B_b) + Z·(1 + Rdif − Tdif)
    #   src_dn/π = B_b·(1 − Rdif − Tdif) + Tdif·(B_b − B_t) − Z·(1 + Rdif − Tdif)
    # and the coefficient combinations factor exactly (om2 = om1(1 + e1)):
    #   1 ∓ Rdif − Tdif = om1·(k·om1 + (γ1 ∓ γ2)(1 + e1))·RT_term
    # so, with γ1 − γ2 ≡ lw_diff_sec(1 − ssa),
    #   Z·(1 + Rdif − Tdif) = ΔB·(om1/τ)·(k·om1 + (γ1+γ2)(1+e1))·RT_term/(γ1+γ2)
    # where om1/τ → k as τ → 0 and stays accurate via expm1. Every factor is
    # a product/sum of non-cancelling terms, bounding the absolute source
    # error by ~eps·ΔB for ALL τ, preserving the real O(τ) emission of thin
    # layers. Only a genuinely empty layer (τ = 0: Rdif = 0, Tdif = 1, no
    # emission) short-circuits here.
    if τ > FT(0)
        ΔB = lev_src_bot - lev_src_top
        γ_sum = γ1 + γ2
        one_p_e1 = 1 + e1
        # layer emissivity factor, 1 − Rdif − Tdif
        emis_fac =
            om1 *
            (k * om1 + lw_diff_sec * (FT(1) - ssa) * one_p_e1) *
            RT_term
        # Z·(1 + Rdif − Tdif), as ΔB times a well-conditioned factor
        ΔBz =
            ΔB * (om1 / τ) * (k * om1 + γ_sum * one_p_e1) * RT_term /
            max(γ_sum, eps(FT))
        src_up = FT(π) * (lev_src_top * emis_fac - Tdif * ΔB + ΔBz)
        src_dn = FT(π) * (lev_src_bot * emis_fac + Tdif * ΔB - ΔBz)
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
    (; flux_up, flux_dn) = flux

    (; albedo, lev_source, sfc_source, src) = src_lw
    (; inc_flux, sfc_emis) = bcs_lw
    FT = eltype(τ)
    @inbounds flux_dn_ilevplus1 =
        isnothing(inc_flux) ? FT(0) : inc_flux[gcol, igpt]
    @inbounds flux_dn[gcol, nlev] = flux_dn_ilevplus1
    # Albedo of lowest level is the surface albedo...
    @inbounds albedo_ilev = FT(1) - sfc_emis[ibnd, gcol]
    @inbounds albedo[gcol, 1] = albedo_ilev
    # ... and source of diffuse radiation is surface emission
    @inbounds src_ilev = FT(π) * sfc_emis[ibnd, gcol] * sfc_source[gcol]
    @inbounds src[gcol, 1] = src_ilev#FT(π) * sfc_emis[ibnd, gcol] * sfc_source[gcol]

    # From bottom to top of atmosphere --
    #   compute albedo and source of upward radiation
    @inbounds lev_src_bot = lev_source[gcol, 1]
    @inbounds for ilev in 1:nlay
        lev_src_top = lev_source[gcol, ilev + 1]
        τ_lay, ssa_lay, g_lay =
            τ[gcol, ilev], ssa[gcol, ilev], g[gcol, ilev]
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
        albedo[gcol, ilev + 1], src[gcol, ilev + 1] =
            albedo_ilevplus1, src_ilevplus1
        lev_src_bot = lev_src_top
        albedo_ilev, src_ilev = albedo_ilevplus1, src_ilevplus1
    end

    # Eq 12, at the top of the domain upwelling diffuse is due to ...
    @inbounds flux_up[gcol, nlev] =
        flux_dn_ilevplus1 * albedo[gcol, nlev] + # ... reflection of incident diffuse and
        src[gcol, nlev]                          # scattering by the direct beam below

    # From the top of the atmosphere downward -- compute fluxes
    @inbounds lev_src_top = lev_source[gcol, nlay + 1]
    ilev = nlay
    @inbounds while ilev ≥ 1
        lev_src_bot, albedo_ilev, src_ilev =
            lev_source[gcol, ilev], albedo[gcol, ilev], src[gcol, ilev]
        τ_lay, ssa_lay, g_lay =
            τ[gcol, ilev], ssa[gcol, ilev], g[gcol, ilev]
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
        flux_up[gcol, ilev] =
            flux_dn_ilev * albedo_ilev + # Equation 12
            src_ilev
        flux_dn[gcol, ilev] = flux_dn_ilev
        flux_dn_ilevplus1 = flux_dn_ilev
        lev_src_top = lev_src_bot
        ilev -= 1
    end
end
