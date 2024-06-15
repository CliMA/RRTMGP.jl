function rte_sw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    max_threads,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt, igpt, ibnd = 1, 1, 1
    cos_zenith = bcs_sw.cos_zenith
    FT = eltype(cos_zenith)
    solar_frac = FT(1)
    @inbounds begin
        flux_up_sw = flux_sw.flux_up
        flux_dn_sw = flux_sw.flux_dn
        flux_net_sw = flux_sw.flux_net

        ClimaComms.@threaded device for gcol in 1:ncol
            if cos_zenith[gcol] > 0 # zero out columns with zenith angle ≥ π/2
                compute_optical_props!(op, as, gcol)
                # call shortwave rte solver
                rte_sw_2stream!(op, src_sw, bcs_sw, flux_sw, solar_frac, igpt, n_gpt, ibnd, nlev, gcol)
                for ilev in 1:nlev
                    flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
                end
            else
                set_flux_to_zero!(flux_sw, gcol)
            end
        end
    end
    return nothing
end

function rte_sw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    max_threads,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    @inbounds begin
        cloud_state = as.cloud_state
        bld_cld_mask = cloud_state isa CloudState
        flux_up_sw = flux_sw.flux_up
        flux_dn_sw = flux_sw.flux_dn
        flux_net_sw = flux_sw.flux_net
        (; flux_up, flux_dn) = flux
        cos_zenith = bcs_sw.cos_zenith
        FT = eltype(flux_up)
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                if cos_zenith[gcol] > 0
                    bld_cld_mask && Optics.build_cloud_mask!(
                        view(cloud_state.mask_sw, :, gcol),
                        view(cloud_state.cld_frac, :, gcol),
                        cloud_state.mask_type,
                    )
                    # compute optical properties
                    compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
                    solar_frac = lookup_sw.solar_src_scaled[igpt]
                    ibnd = lookup_sw.band_data.major_gpt2bnd[igpt]
                    # call rte shortwave solver
                    rte_sw_2stream!(op, src_sw, bcs_sw, flux, solar_frac, igpt, n_gpt, ibnd, nlev, gcol)
                    if igpt == 1
                        map!(x -> x, view(flux_up_sw, :, gcol), view(flux_up, :, gcol))
                        map!(x -> x, view(flux_dn_sw, :, gcol), view(flux_dn, :, gcol))
                    else
                        for ilev in 1:nlev
                            @inbounds flux_up_sw[ilev, gcol] += flux_up[ilev, gcol]
                            @inbounds flux_dn_sw[ilev, gcol] += flux_dn[ilev, gcol]
                        end
                    end
                else
                    set_flux_to_zero!(flux_sw, gcol)
                end
            end
        end

        ClimaComms.@threaded device for gcol in 1:ncol
            if cos_zenith[gcol] > 0
                for ilev in 1:nlev
                    flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
                end
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
function sw_2stream_coeffs(τ::FT, ssa::FT, g::FT, μ₀::FT) where {FT}
    k_min = sqrt(eps(FT)) #FT(1e4 * eps(FT)) # Suggestion from Chiel van Heerwaarden
    # Zdunkowski Practical Improved Flux Method "PIFM"
    #  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    γ1 = (FT(8) - ssa * (FT(5) + FT(3) * g)) * FT(0.25)
    γ2 = FT(3) * (ssa * (FT(1) - g)) * FT(0.25)
    γ3 = (FT(2) - (FT(3) * μ₀) * g) * FT(0.25)
    γ4 = FT(1) - γ3
    α1 = γ1 * γ4 + γ2 * γ3                          # Eq. 16
    α2 = γ1 * γ3 + γ2 * γ4                          # Eq. 17
    k = sqrt(max((γ1 - γ2) * (γ1 + γ2), k_min))

    #τ = max(τ, FT(0))
    exp_minusktau = exp(-τ * k)
    exp_minus2ktau = exp_minusktau * exp_minusktau

    # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    RT_term = FT(1) / (k * (FT(1) + exp_minus2ktau) + γ1 * (FT(1) - exp_minus2ktau))

    Rdif = RT_term * γ2 * (FT(1) - exp_minus2ktau) # Eqn. 25
    Tdif = RT_term * FT(2) * k * exp_minusktau     # Eqn. 26

    # Transmittance of direct, unscattered beam. Also used below
    T₀ = Tnoscat = exp(-τ / μ₀)

    # Direct reflect and transmission
    k_μ = k * μ₀
    k_γ3 = k * γ3
    k_γ4 = k * γ4

    # Equation 14, multiplying top and bottom by exp(-k*tau)
    #   and rearranging to avoid div by 0.
    if abs(FT(1) - k_μ * k_μ) ≥ eps(FT)
        RT_term = ssa * RT_term / (FT(1) - k_μ * k_μ)
    else
        RT_term = ssa * RT_term / eps(FT)
    end

    Rdir_unconstrained =
        RT_term * (
            (FT(1) - k_μ) * (α2 + k_γ3) - (FT(1) + k_μ) * (α2 - k_γ3) * exp_minus2ktau -
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
            (FT(1) + k_μ) * (α1 + k_γ4) * T₀ - (FT(1) - k_μ) * (α1 - k_γ4) * exp_minus2ktau * T₀ -
            FT(2) * (k_γ4 + α1 * k_μ) * exp_minusktau
        )
    # Final check that energy is not spuriously created, by recognizing that
    # the beam can either be reflected, penetrate unscattered to the base of a layer, 
    # or penetrate through but be scattered on the way - the rest is absorbed
    # Makes the equations safer in single precision. Credit: Robin Hogan, Peter Ukkonen
    Rdir = max(FT(0), min(Rdir_unconstrained, (FT(1) - T₀)))
    Tdir = max(FT(0), min(Tdir_unconstrained, (FT(1) - T₀ - Rdir)))
    return (Rdir, Tdir, Tnoscat, Rdif, Tdif)
end

# Direct-beam and source for diffuse radiation
function get_flux_dn_dir(τ, μ₀, flux_dn_dir_top, lev)
    nlay = length(τ)
    τ_sum = zero(eltype(τ))
    for ilev in nlay:-1:lev
        τ_sum += τ[ilev]
    end
    #τ_sum = max(τ_sum, zero(eltype(τ)))
    return flux_dn_dir_top * exp(-τ_sum / μ₀)
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
function rte_sw_2stream!(
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
    τ_sum = FT(0)
    for ilay in 1:nlay
        @inbounds τ_sum += τ[ilay, gcol]
    end
    #τ_sum = max(τ_sum, FT(0))
    # Direct-beam and source for diffuse radiation
    flux_dn_dir_top = toa_flux * solar_frac * μ₀
    flux_dn_dir_bot = flux_dn_dir_top * exp(-τ_sum / μ₀) # store value at surface
    @inbounds flux_dn_dir[1, gcol] = flux_dn_dir_bot # store value at surface
    sfc_source = flux_dn_dir_bot * sfc_alb_direct

    @inbounds flux_dn[nlev, gcol] = FT(0) # set to incoming flux when provided?
    # Albedo of lowest level is the surface albedo...
    @inbounds surface_albedo = albedo[1, gcol] = bcs_sw.sfc_alb_diffuse[ibnd, gcol]
    # ... and source of diffuse radiation is surface emission
    @inbounds src[1, gcol] = sfc_source
    # From bottom to top of atmosphere --
    #   compute albedo and source of upward radiation
    τ_cum = τ_sum
    albedo_ilev, src_ilev = surface_albedo, sfc_source
    @inbounds for ilev in 1:nlay
        τ_ilev, ssa_ilev, g_ilev = τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        (Rdir, Tdir, _, Rdif, Tdif) = sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀)
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        albedo_ilevplus1 = Rdif + Tdif * Tdif * albedo_ilev * denom # Equation 9
        # 
        # Equation 11 -- source is emitted upward radiation at top of layer plus
        # radiation emitted at bottom of layer,
        # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
        τ_cum -= τ_ilev
        τ_cum = max(τ_cum, FT(0))
        flux_dn_dir_ilevplus1 = flux_dn_dir_top * exp(-τ_cum / μ₀)
        src_up_ilev = Rdir * flux_dn_dir_ilevplus1 #flux_dn_dir[ilev + 1]
        src_dn_ilev = Tdir * flux_dn_dir_ilevplus1 #flux_dn_dir[ilev + 1]
        src_ilevplus1 = src_up_ilev + Tdif * denom * (src_ilev + albedo_ilev * src_dn_ilev)
        albedo[ilev + 1, gcol], src[ilev + 1, gcol] = albedo_ilevplus1, src_ilevplus1
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
    τ_cum = FT(0)

    ilev = nlay
    @inbounds while ilev ≥ 1
        τ_ilev, ssa_ilev, g_ilev = τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        albedo_ilev, src_ilev = albedo[ilev, gcol], src[ilev, gcol]
        (_, Tdir, _, Rdif, Tdif) = sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀)
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        src_dn_ilev = Tdir * flux_dn_dir_top * exp(-τ_cum / μ₀)
        τ_cum += τ_ilev
        #τ_cum = max(τ_cum, FT(0))
        flux_dn_ilev = (Tdif * flux_dn_ilevplus1 + # Equation 13
                        Rdif * src_ilev +
                        src_dn_ilev) * denom
        flux_up[ilev, gcol] =
            flux_dn_ilev * albedo_ilev + # Equation 12
            src_ilev
        flux_dn[ilev, gcol] = flux_dn_ilev + flux_dn_dir_top * exp(-τ_cum / μ₀)
        flux_dn_ilevplus1 = flux_dn_ilev
        ilev -= 1
        if isnan(flux_dn[ilev, gcol])
            my_print(flux_dn[ilev, gcol],
            ilev,
            gcol,
            τ_cum,
            μ₀,
            exp(-τ_cum / μ₀),
            denom,
            flux_dn_ilev,
            Tdif,
            Rdif,
            src_ilev,
            src_dn_ilev,
            flux_dn[nlev, gcol])
        end
    end
    return nothing
end

function my_print end
