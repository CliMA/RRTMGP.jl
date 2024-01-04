function rte_sw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    max_threads,
    as::GrayAtmosphericState{FT},
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt, igpt, ibnd = 1, 1, UInt8(1)
    solar_frac = FT(1)
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol, igpt, nothing, nothing)
            # call shortwave rte solver
            rte_sw_2stream!(op, src_sw, bcs_sw, flux_sw, solar_frac, igpt, n_gpt, ibnd, nlev, gcol)
        end
    end
    return nothing
end

function rte_sw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    max_threads,
    as::GrayAtmosphericState{FT},
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux_sw, op, bcs_sw, src_sw, nlay, ncol, as)
    @cuda threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(args...)
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    nlay,
    ncol,
    as::GrayAtmosphericState{FT},
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt, ibnd = 1, 1, UInt8(1)
    solar_frac = FT(1)
    if gcol ≤ ncol
        @inbounds begin
            set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol, igpt, nothing, nothing)
            # call shortwave rte solver
            rte_sw_2stream!(op, src_sw, bcs_sw, flux_sw, solar_frac, igpt, n_gpt, ibnd, nlev, gcol)
        end
    end
    return nothing
end

function rte_sw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    max_threads,
    as::AtmosphericState{FT},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    @inbounds begin
        bld_cld_mask = as.cld_mask_type isa AbstractCloudMask
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                bld_cld_mask && Optics.build_cloud_mask!(
                    view(as.cld_mask_sw, :, gcol),
                    view(as.cld_frac, :, gcol),
                    as.cld_mask_type,
                )
                igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                # compute optical properties
                compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
                solar_frac = lookup_sw.solar_src_scaled[igpt]
                ibnd = lookup_sw.major_gpt2bnd[igpt]
                # call rte shortwave solver
                rte_sw_2stream!(op, src_sw, bcs_sw, flux, solar_frac, igpt, n_gpt, ibnd, nlev, gcol)
                n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
            end
        end
    end
    return nothing
end

function rte_sw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    max_threads,
    as::AtmosphericState{FT},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux, flux_sw, op, bcs_sw, src_sw, nlay, ncol, as, lookup_sw, lookup_sw_cld)
    @cuda threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(args...)
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    nlay,
    ncol,
    as::AtmosphericState{FT},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.major_gpt2bnd)
    if gcol ≤ ncol
        flux_up_sw = flux_sw.flux_up
        flux_dn_sw = flux_sw.flux_dn
        flux_net_sw = flux_sw.flux_net
        flux_up = flux.flux_up
        flux_dn = flux.flux_dn
        @inbounds for ilev in 1:nlev
            flux_up_sw[ilev, gcol] = FT(0)
            flux_dn_sw[ilev, gcol] = FT(0)
        end
        @inbounds for igpt in 1:n_gpt
            # set cloud mask, if applicable
            if as.cld_mask_type isa AbstractCloudMask
                Optics.build_cloud_mask!(view(as.cld_mask_sw, :, gcol), view(as.cld_frac, :, gcol), as.cld_mask_type)
            end
            # compute optical properties
            compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
            solar_frac = lookup_sw.solar_src_scaled[igpt]
            ibnd = lookup_sw.major_gpt2bnd[igpt]
            # rte shortwave solver
            rte_sw_2stream!(op, src_sw, bcs_sw, flux, solar_frac, igpt, n_gpt, ibnd, nlev, gcol)
            if n_gpt > 1
                for ilev in 1:nlev
                    @inbounds flux_up_sw[ilev, gcol] += flux_up[ilev, gcol]
                    @inbounds flux_dn_sw[ilev, gcol] += flux_dn[ilev, gcol]
                end
            end
        end
        for ilev in 1:nlev
            @inbounds flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
        end
    end
    return nothing
end

function sw_2stream_coeffs(τ::FT, ssa::FT, g::FT, μ₀::FT) where {FT}
    k_min = FT(1e4 * eps(FT)) # Suggestion from Chiel van Heerwaarden
    # Zdunkowski Practical Improved Flux Method "PIFM"
    #  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    γ1 = (FT(8) - ssa * (FT(5) + FT(3) * g)) * FT(0.25)
    γ2 = FT(3) * (ssa * (FT(1) - g)) * FT(0.25)
    γ3 = (FT(2) - (FT(3) * μ₀) * g) * FT(0.25)
    γ4 = FT(1) - γ3
    α1 = γ1 * γ4 + γ2 * γ3                          # Eq. 16
    α2 = γ1 * γ3 + γ2 * γ4                          # Eq. 17
    k = sqrt(max((γ1 - γ2) * (γ1 + γ2), k_min))

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
        ibnd::UInt8,
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
    ibnd::UInt8,
    nlev::Int,
    gcol::Int,
) where {FT}
    nlay = nlev - 1
    toa_flux = bcs_sw.toa_flux[gcol]
    sfc_alb_direct = bcs_sw.sfc_alb_direct[ibnd, gcol]
    μ₀ = cos(bcs_sw.zenith[gcol])
    τ_sum = FT(0)
    for ilay in 1:nlay
        @inbounds τ_sum += τ[ilay, gcol]
    end
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
        τ_ilev = τ[ilev, gcol]
        (Rdir, Tdir, _, Rdif, Tdif) = sw_2stream_coeffs(τ_ilev, ssa[ilev, gcol], g[ilev, gcol], μ₀)
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        albedo_ilevplus1 = Rdif + Tdif * Tdif * albedo_ilev * denom # Equation 9
        # 
        # Equation 11 -- source is emitted upward radiation at top of layer plus
        # radiation emitted at bottom of layer,
        # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
        τ_cum -= τ_ilev
        flux_dn_dir_ilevplus1 = flux_dn_dir_top * exp(-τ_cum / μ₀)
        src_up_ilev = Rdir * flux_dn_dir_ilevplus1 #flux_dn_dir[ilev + 1]
        src_dn_ilev = Tdir * flux_dn_dir_ilevplus1 #flux_dn_dir[ilev + 1]
        src_ilev = src_up_ilev + Tdif * denom * (src_ilev + albedo_ilev * src_dn_ilev)
        src[ilev + 1, gcol] = src_ilev
        albedo[ilev + 1, gcol] = albedo_ilev = albedo_ilevplus1
    end
    # Eq 12, at the top of the domain upwelling diffuse is due to ...
    @inbounds flux_up[nlev, gcol] =
        flux_dn[nlev, gcol] * albedo[nlev, gcol] + # ... reflection of incident diffuse and
        src[nlev, gcol]                          # scattering by the direct beam below

    # From the top of the atmosphere downward -- compute fluxes
    @inbounds flux_dn_ilevplus1 = flux_dn[nlev, gcol]
    @inbounds flux_dn[nlev, gcol] += flux_dn_dir_top
    τ_cum = FT(0)

    @inbounds for ilev in nlay:-1:1
        τ_ilev, albedo_ilev, src_ilev = τ[ilev, gcol], albedo[ilev, gcol], src[ilev, gcol]
        (_, Tdir, _, Rdif, Tdif) = sw_2stream_coeffs(τ_ilev, ssa[ilev, gcol], g[ilev, gcol], μ₀)
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        src_dn_ilev = Tdir * flux_dn_dir_top * exp(-τ_cum / μ₀)
        τ_cum += τ_ilev
        flux_dn_ilev = (Tdif * flux_dn_ilevplus1 + # Equation 13
                        Rdif * src_ilev +
                        src_dn_ilev) * denom
        flux_up[ilev, gcol] =
            flux_dn_ilev * albedo_ilev + # Equation 12
            src_ilev
        flux_dn[ilev, gcol] = flux_dn_ilev + flux_dn_dir_top * exp(-τ_cum / μ₀)
        flux_dn_ilevplus1 = flux_dn_ilev
    end
    return nothing
end
