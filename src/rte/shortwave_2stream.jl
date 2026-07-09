function rte_sw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt, igpt, ibnd = 1, 1, 1
    cos_zenith = bcs_sw.cos_zenith
    FT = eltype(cos_zenith)
    solar_frac = FT(1)
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            compute_optical_props!(op, as, gcol)
            if cos_zenith[gcol] > 0
                # call shortwave rte solver
                rte_sw_2stream!(
                    op,
                    src_sw,
                    bcs_sw,
                    flux_sw,
                    solar_frac,
                    igpt,
                    n_gpt,
                    ibnd,
                    nlev,
                    gcol,
                )
                compute_net_flux!(flux_sw, gcol, nlev)
            else # zero out columns with zenith angle ≥ π/2
                set_flux_to_zero!(flux_sw, gcol, nlev)
            end
        end
    end
    return nothing
end

# Device-agnostic per-(g-point, column) body, shared by the CPU driver below
# and the CUDA kernel in ext/cuda. The cloud mask and optical properties are
# computed for every column (the cloud-cover diagnostic counts night columns
# too); the solve and flux accumulation run only for day columns — night
# columns are zeroed once, after the g-point loop. Returns whether this
# g-point had any cloudy layer.
@inline function sw_2stream_gpt_col!(
    igpt,
    gcol,
    flux,
    flux_sw,
    band_flux,
    op,
    bcs_sw,
    src_sw,
    as,
    lookup_sw,
    lookup_sw_cld,
    lookup_sw_aero,
    μ₀,
    ibnd,
    n_gpt,
    nlev,
)
    cloudy = _build_cloud_mask!(as.cloud_state, Val(:mask_sw), gcol)
    compute_optical_props!(
        op,
        as,
        gcol,
        igpt,
        lookup_sw,
        lookup_sw_cld,
        lookup_sw_aero,
    )
    if μ₀ > 0
        @inbounds solar_frac = lookup_sw.solar_src_scaled[igpt]
        rte_sw_2stream!(
            op,
            src_sw,
            bcs_sw,
            flux,
            solar_frac,
            igpt,
            n_gpt,
            ibnd,
            nlev,
            gcol,
        )
        _accumulate_fluxes!(flux_sw, flux, gcol, nlev, igpt)
        # retain this g-point's contribution in its band (no-op when off)
        accumulate_band_flux!(
            band_flux,
            flux.flux_up,
            flux.flux_dn,
            gcol,
            ibnd,
            nlev,
        )
    end
    return cloudy
end

function rte_sw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    # zero the optional per-band accumulator (no-op when spectral fluxes are off)
    set_band_flux_to_zero!(band_flux)
    @inbounds begin
        (; cloud_state, aerosol_state) = as
        cos_zenith = bcs_sw.cos_zenith
        track_cld_cover =
            cloud_state isa CloudState && !isnothing(cloud_state.cld_cover_sw)
        FT = eltype(flux_sw.flux_up)
        track_cld_cover && (cloud_state.cld_cover_sw .= FT(0))
        if aerosol_state isa AerosolState
            ClimaComms.@threaded device for gcol in 1:ncol
                _compute_aero_mask!(aerosol_state, gcol)
            end
        end
        for igpt in 1:n_gpt
            ibnd = lookup_sw.band_data.major_gpt2bnd[igpt]
            ClimaComms.@threaded device for gcol in 1:ncol
                cloudy = sw_2stream_gpt_col!(
                    igpt,
                    gcol,
                    flux,
                    flux_sw,
                    band_flux,
                    op,
                    bcs_sw,
                    src_sw,
                    as,
                    lookup_sw,
                    lookup_sw_cld,
                    lookup_sw_aero,
                    cos_zenith[gcol],
                    ibnd,
                    n_gpt,
                    nlev,
                )
                track_cld_cover &&
                    (cloud_state.cld_cover_sw[gcol] += FT(cloudy))
            end
        end

        # normalize cloud cover by number of g-points
        if track_cld_cover
            ClimaComms.@threaded device for gcol in 1:ncol
                cloud_state.cld_cover_sw[gcol] /= n_gpt
            end
        end
        ClimaComms.@threaded device for gcol in 1:ncol
            if cos_zenith[gcol] > 0
                compute_net_flux!(flux_sw, gcol, nlev)
            else # zero out columns with zenith angle ≥ π/2
                set_flux_to_zero!(flux_sw, gcol, nlev)
            end
        end
    end
    return nothing
end

"""
    sw_2stream_coeffs(τ::FT, ssa::FT, g::FT, μ₀::FT) where {FT}
Computes cell properties (transmittance and reflectance) for direct and diffuse radiation
Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
 
Equations are developed in Meador and Weaver, 1980,
doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
"""
@inline function sw_2stream_coeffs(τ::FT, ssa::FT, g::FT, μ₀::FT) where {FT}
    k_min = Numerics.k_min(FT)
    # Zdunkowski Practical Improved Flux Method "PIFM"
    #  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    γ1 = (FT(8) - ssa * (FT(5) + FT(3) * g)) * FT(0.25)
    γ2 = FT(3) * (ssa * (FT(1) - g)) * FT(0.25)
    γ3 = (FT(2) - (FT(3) * μ₀) * g) * FT(0.25)
    γ4 = FT(1) - γ3
    α1 = γ1 * γ4 + γ2 * γ3                          # Eq. 16
    α2 = γ1 * γ3 + γ2 * γ4                          # Eq. 17
    # For PIFM, γ1 − γ2 ≡ 2(1 − ssa) exactly; use the identity rather than the
    # difference of two rounded O(1) numbers, which loses ~eps absolute and
    # gives O(1%) errors in k for bright scattering (ssa → 1) at Float32.
    k = sqrt(max(FT(2) * (FT(1) - ssa) * (γ1 + γ2), k_min))

    exp_minusktau = exp(-τ * k)
    exp_minus2ktau = exp_minusktau * exp_minusktau
    # 1 − e^{−kτ} via expm1 (accurate for small kτ), and 1 − e^{−2kτ} from it
    # by the exact factorization (1 − e)(1 + e) — the naive 1 − e² loses
    # ~eps/(2kτ) relative accuracy for optically thin layers.
    om1 = -expm1(-τ * k)
    one_minus_e2kt = om1 * (FT(1) + exp_minusktau)

    # Formulated to avoid rounding errors when k, gamma1 are of very different magnitudes
    RT_term = FT(1) / (k * (FT(1) + exp_minus2ktau) + γ1 * one_minus_e2kt)

    Rdif = RT_term * γ2 * one_minus_e2kt       # Eqn. 25
    Tdif = RT_term * FT(2) * k * exp_minusktau # Eqn. 26

    # Transmittance of direct, unscattered beam. Also used below
    T₀ = Tnoscat = exp(-τ / max(μ₀, Numerics.μ₀_min(FT)))

    # Direct reflect and transmission
    k_μ = k * μ₀
    k_μ2 = k_μ * k_μ
    # Equations 14-15 have a removable singularity at k·μ₀ = 1. Nudge k·μ₀
    # off resonance so numerator and denominator stay mutually consistent.
    # See Numerics.resonance_window for the sqrt(eps) choice.
    if abs(FT(1) - k_μ2) < Numerics.resonance_window(FT)
        k_μ2 = FT(1) - Numerics.resonance_window(FT)
        k_μ = sqrt(k_μ2)
    end
    k_γ3 = k * γ3
    k_γ4 = k * γ4

    # Equation 14, multiplying top and bottom by exp(-k*tau)
    #   and rearranging to avoid div by 0.
    RT_term = ssa * RT_term / (FT(1) - k_μ2)

    Rdir_unconstrained =
        RT_term * (
            (FT(1) - k_μ) * (α2 + k_γ3) -
            (FT(1) + k_μ) * (α2 - k_γ3) * exp_minus2ktau -
            FT(2) * (k_γ3 - α2 * k_μ) * exp_minusktau * T₀
        )
    #
    # Equation 15, multiplying top and bottom by exp(-k*tau),
    #   multiplying through by exp(-tau/mu0) to
    #   prefer underflow to overflow
    # Omitting direct transmittance
    #
    Tdir_unconstrained =
        -RT_term * (
            (FT(1) + k_μ) * (α1 + k_γ4) * T₀ -
            (FT(1) - k_μ) * (α1 - k_γ4) * exp_minus2ktau * T₀ -
            FT(2) * (k_γ4 + α1 * k_μ) * exp_minusktau
        )
    # Final check that energy is not spuriously created, by recognizing that
    # the beam can either be reflected, penetrate unscattered to the base of a layer, 
    # or penetrate through but be scattered on the way - the rest is absorbed
    # Makes the equations safer in single precision. Credit: Robin Hogan, Peter Ukkonen
    # Note: Unlike Eq. 9 and 10 in Hogan et al. (2024), "Improving the Two-Stream 
    # Approximation in RRTMGP" (doi:10.1029/2023MS003932) where the upper bounds are 
    # scaled by μ₀ (the cosine of the solar zenith angle), here they are not scaled. 
    # This is because in RRTMGP, Rdir and Tdir are defined as fractional quantities 
    # relative to the horizontal incident direct flux (which already includes the μ₀ 
    # factor), whereas in ecRad (the model in the paper), they are normalized relative 
    # to the perpendicular incident intensity. Thus, the maximum fractional reflectance 
    # Rdir is bounded by 1 - T₀, and must not be scaled by μ₀.
    Rdir = max(FT(0), min(Rdir_unconstrained, (FT(1) - T₀)))
    Tdir = max(FT(0), min(Tdir_unconstrained, (FT(1) - T₀ - Rdir)))
    return (Rdir, Tdir, Tnoscat, Rdif, Tdif)
end

"""
    rte_sw_2stream!(
        (; τ, ssa, g)::TwoStream,
        (; albedo, src)::SourceSW2Str,
        bcs_sw::SwBCs,
        (; flux_up, flux_dn, flux_dn_dir)::FluxSW,
        solar_frac::FT,
        igpt::Int,
        n_gpt::Int,
        ibnd::Int,
        nlev::Int,
        gcol::Int,
    ) where {FT}

Two stream solver for the shortwave problem.

Transport of diffuse radiation through a vertically layered atmosphere.
Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
"""
@inline function rte_sw_2stream!(
    (; τ, ssa, g)::TwoStream,
    (; albedo, src)::SourceSW2Str,
    bcs_sw::SwBCs,
    (; flux_up, flux_dn, flux_dn_dir)::FluxSW,
    solar_frac::FT,
    igpt::Int,
    n_gpt::Int,
    ibnd::Int,
    nlev::Int,
    gcol::Int,
) where {FT}
    nlay = nlev - 1
    @inbounds begin
        toa_flux = bcs_sw.toa_flux[gcol]
        sfc_alb_direct = bcs_sw.sfc_alb_direct[ibnd, gcol]
        μ₀ = bcs_sw.cos_zenith[gcol]
    end
    # Direct-beam profile, computed top-down with an addition-built
    # cumulative optical depth and stored for the two passes below.
    flux_dn_dir_top = toa_flux * solar_frac * μ₀
    inv_μ₀ = FT(1) / max(μ₀, Numerics.μ₀_min(FT))
    @inbounds flux_dn_dir[nlev, gcol] = flux_dn_dir_top
    τ_cum = FT(0)
    @inbounds for ilev in nlay:-1:1
        τ_cum += τ[ilev, gcol]
        flux_dn_dir[ilev, gcol] = flux_dn_dir_top * exp(-τ_cum * inv_μ₀)
    end
    flux_dn_dir_bot = @inbounds flux_dn_dir[1, gcol] # surface value
    sfc_source = flux_dn_dir_bot * sfc_alb_direct

    @inbounds flux_dn[nlev, gcol] = FT(0) # set to incoming flux when provided?
    # Albedo of lowest level is the surface albedo...
    @inbounds surface_albedo =
        albedo[1, gcol] = bcs_sw.sfc_alb_diffuse[ibnd, gcol]
    # ... and source of diffuse radiation is surface emission
    @inbounds src[1, gcol] = sfc_source
    # From bottom to top of atmosphere --
    #   compute albedo and source of upward radiation
    albedo_ilev, src_ilev = surface_albedo, sfc_source
    @inbounds for ilev in 1:nlay
        τ_ilev, ssa_ilev, g_ilev =
            τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        (Rdir, Tdir, _, Rdif, Tdif) =
            sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀)
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        albedo_ilevplus1 = Rdif + Tdif * Tdif * albedo_ilev * denom # Equation 9
        #
        # Equation 11 -- source is emitted upward radiation at top of layer plus
        # radiation emitted at bottom of layer,
        # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
        flux_dn_dir_ilevplus1 = flux_dn_dir[ilev + 1, gcol]
        src_up_ilev = Rdir * flux_dn_dir_ilevplus1
        src_dn_ilev = Tdir * flux_dn_dir_ilevplus1
        src_ilevplus1 =
            src_up_ilev +
            Tdif * denom * (src_ilev + albedo_ilev * src_dn_ilev)
        albedo[ilev + 1, gcol], src[ilev + 1, gcol] =
            albedo_ilevplus1, src_ilevplus1
        albedo_ilev = albedo_ilevplus1
        src_ilev = src_ilevplus1
    end
    # Eq 12, at the top of the domain upwelling diffuse is due to ...
    @inbounds flux_up[nlev, gcol] =
        flux_dn[nlev, gcol] * albedo[nlev, gcol] + # ... reflection of incident diffuse and
        src[nlev, gcol]                          # scattering by the direct beam below

    # From the top of the atmosphere downward -- compute fluxes
    @inbounds flux_dn_ilevplus1 = flux_dn[nlev, gcol]
    @inbounds flux_dn[nlev, gcol] += flux_dn_dir_top

    ilev = nlay
    @inbounds while ilev ≥ 1
        τ_ilev, ssa_ilev, g_ilev =
            τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        albedo_ilev, src_ilev = albedo[ilev, gcol], src[ilev, gcol]
        (_, Tdir, _, Rdif, Tdif) =
            sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀)
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        src_dn_ilev = Tdir * flux_dn_dir[ilev + 1, gcol]
        flux_dn_ilev =
            (Tdif * flux_dn_ilevplus1 + # Equation 13
             Rdif * src_ilev +
             src_dn_ilev) * denom
        flux_up[ilev, gcol] =
            flux_dn_ilev * albedo_ilev + # Equation 12
            src_ilev
        flux_dn[ilev, gcol] = flux_dn_ilev + flux_dn_dir[ilev, gcol]
        flux_dn_ilevplus1 = flux_dn_ilev
        ilev -= 1
    end
    return nothing
end
