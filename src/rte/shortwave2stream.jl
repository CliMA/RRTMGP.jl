"""
    rte_sw_2stream_solve!(
        context,
        flux::FluxSW{FT},
        flux_sw::FluxSW{FT},
        op::TwoStream{FT},
        bcs_sw::SwBCs{FT},
        src_sw::SourceSW2Str{FT},
        major_gpt2bnd::AbstractArray{UInt8,1},
        solar_src_scaled::AbstractArray{FT,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lookup_sw::Union{LookUpSW, Nothing} = nothing,
        lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    ) where {FT<:AbstractFloat}

Two stream solver for the shortwave problem.
"""
function rte_sw_2stream_solve!(
    context,
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    major_gpt2bnd::AbstractArray{UInt8, 1},
    solar_src_scaled::AbstractArray{FT, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (
            flux,
            flux_sw,
            op,
            bcs_sw,
            src_sw,
            nlay,
            ncol,
            major_gpt2bnd,
            solar_src_scaled,
            as,
            lookup_sw,
            lookup_sw_cld,
        )
        @cuda threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(args...)
    else
        @inbounds begin
            bld_cld_mask = as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
            # setting references for flux_sw
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                    bld_cld_mask && Optics.build_cloud_mask!(
                        as.cld_mask_sw,
                        as.cld_frac,
                        as.random_sw,
                        gcol,
                        igpt,
                        as.cld_mask_type,
                    )
                    compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
                    sw_two_stream!(op, src_sw, bcs_sw, gcol, nlay) # Cell properties: transmittance and reflectance for direct and diffuse radiation
                    sw_source_2str!(src_sw, bcs_sw, gcol, flux, igpt, solar_src_scaled, major_gpt2bnd, nlay) # Direct-beam and source for diffuse radiation
                    adding_sw!(src_sw, bcs_sw, gcol, flux, igpt, major_gpt2bnd, nlev) # Transport
                    n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
                end
            end
        end
    end
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
    major_gpt2bnd::AbstractArray{UInt8, 1},
    solar_src_scaled::AbstractArray{FT, 1},
    as::AbstractAtmosphericState{FT},
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        @inbounds for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
            if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
                Optics.build_cloud_mask!(as.cld_mask_sw, as.cld_frac, as.random_sw, gcol, igpt, as.cld_mask_type)
            end
            compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
            sw_two_stream!(op, src_sw, bcs_sw, gcol, nlay) # Cell properties: transmittance and reflectance for direct and diffuse radiation
            sw_source_2str!(src_sw, bcs_sw, gcol, flux, igpt, solar_src_scaled, major_gpt2bnd, nlay) # Direct-beam and source for diffuse radiation
            adding_sw!(src_sw, bcs_sw, gcol, flux, igpt, major_gpt2bnd, nlev) # Transport
            n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
        end
    end
    return nothing
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
    gcol::Int,
    nlay::Int,
) where {FT <: AbstractFloat}
    zenith = bcs_sw.zenith
    (; τ, ssa, g) = op
    (; Rdif, Tdif, Rdir, Tdir, Tnoscat) = src_sw
    k_min = FT(1e4 * eps(FT)) # Suggestion from Chiel van Heerwaarden
    μ₀ = FT(cos(zenith[gcol]))

    @inbounds for glay in 1:nlay
        # Zdunkowski Practical Improved Flux Method "PIFM"
        #  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
        γ1 = (FT(8) - ssa[glay, gcol] * (FT(5) + FT(3) * g[glay, gcol])) * FT(0.25)
        γ2 = FT(3) * (ssa[glay, gcol] * (FT(1) - g[glay, gcol])) * FT(0.25)
        γ3 = (FT(2) - (FT(3) * cos(zenith[gcol])) * g[glay, gcol]) * FT(0.25)
        γ4 = FT(1) - γ3
        α1 = γ1 * γ4 + γ2 * γ3                          # Eq. 16
        α2 = γ1 * γ3 + γ2 * γ4                          # Eq. 17
        k = sqrt(max((γ1 - γ2) * (γ1 + γ2), k_min))

        exp_minusktau = exp(-τ[glay, gcol] * k)
        exp_minus2ktau = exp_minusktau * exp_minusktau

        # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
        RT_term = FT(1) / (k * (FT(1) + exp_minus2ktau) + γ1 * (FT(1) - exp_minus2ktau))

        @inbounds Rdif[glay, gcol] = RT_term * γ2 * (FT(1) - exp_minus2ktau) # Eqn. 25
        @inbounds Tdif[glay, gcol] = RT_term * FT(2) * k * exp_minusktau     # Eqn. 26

        # Transmittance of direct, unscattered beam. Also used below
        @inbounds T₀ = Tnoscat[glay, gcol] = exp(-τ[glay, gcol] / cos(zenith[gcol]))

        # Direct reflect and transmission
        k_μ = k * μ₀
        k_γ3 = k * γ3
        k_γ4 = k * γ4

        # Equation 14, multiplying top and bottom by exp(-k*tau)
        #   and rearranging to avoid div by 0.
        if abs(FT(1) - k_μ * k_μ) ≥ eps(FT)
            RT_term = ssa[glay, gcol] * RT_term / (FT(1) - k_μ * k_μ)
        else
            RT_term = ssa[glay, gcol] * RT_term / eps(FT)
        end

        @inbounds Rdir_unconstrained =
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
        @inbounds Tdir_unconstrained =
            -RT_term * (
                (FT(1) + k_μ) * (α1 + k_γ4) * T₀ - (FT(1) - k_μ) * (α1 - k_γ4) * exp_minus2ktau * T₀ -
                FT(2) * (k_γ4 + α1 * k_μ) * exp_minusktau
            )
        # Final check that energy is not spuriously created, by recognizing that
        # the beam can either be reflected, penetrate unscattered to the base of a layer, 
        # or penetrate through but be scattered on the way - the rest is absorbed
        # Makes the equations safer in single precision. Credit: Robin Hogan, Peter Ukkonen
        Rdir[glay, gcol] = max(FT(0), min(Rdir_unconstrained, (FT(1) - T₀)))
        Tdir[glay, gcol] = max(FT(0), min(Tdir_unconstrained, (FT(1) - T₀ - Rdir[glay, gcol])))
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
    gcol::Int,
    flux::FluxSW{FT},
    igpt::Int,
    solar_src_scaled::AbstractArray{FT, 1},
    major_gpt2bnd::AbstractArray{UInt8, 1},
    nlay::Int,
) where {FT <: AbstractFloat}
    ibnd = major_gpt2bnd[igpt]
    solar_frac = solar_src_scaled[igpt]
    (; toa_flux, zenith, sfc_alb_direct) = bcs_sw
    (; Rdir, Tdir, Tnoscat, src_up, src_dn, sfc_source) = src_sw
    flux_dn_dir = flux.flux_dn_dir

    # layer index = level index
    # previous level is up (+1)
    @inbounds flux_dn_dir[nlay + 1, gcol] = toa_flux[gcol] * solar_frac * cos(zenith[gcol])
    @inbounds for ilev in nlay:-1:1
        src_up[ilev, gcol] = Rdir[ilev, gcol] * flux_dn_dir[ilev + 1, gcol]
        src_dn[ilev, gcol] = Tdir[ilev, gcol] * flux_dn_dir[ilev + 1, gcol]
        flux_dn_dir[ilev, gcol] = Tnoscat[ilev, gcol] * flux_dn_dir[ilev + 1, gcol]
    end
    @inbounds sfc_source[gcol] = flux_dn_dir[1, gcol] * sfc_alb_direct[ibnd, gcol]
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
    gcol::Int,
    flux::FluxSW{FT},
    igpt::Int,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    nlev::Int,
) where {FT <: AbstractFloat}
    ibnd = major_gpt2bnd[igpt]
    nlay = nlev - 1
    n_gpt = length(major_gpt2bnd)
    # destructuring
    (; flux_up, flux_dn, flux_dn_dir, flux_net) = flux
    (; albedo, sfc_source, Rdif, Tdif, src_up, src_dn, src) = src_sw

    @inbounds for ilev in 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    # Albedo of lowest level is the surface albedo...
    @inbounds albedo[1, gcol] = bcs_sw.sfc_alb_diffuse[ibnd, gcol]
    # ... and source of diffuse radiation is surface emission
    @inbounds src[1, gcol] = sfc_source[gcol]
    # From bottom to top of atmosphere --
    #   compute albedo and source of upward radiation
    @inbounds for ilev in 1:nlay
        denom = FT(1) / (FT(1) - Rdif[ilev, gcol] * albedo[ilev, gcol])  # Eq 10
        albedo[ilev + 1, gcol] = Rdif[ilev, gcol] + Tdif[ilev, gcol] * Tdif[ilev, gcol] * albedo[ilev, gcol] * denom # Equation 9
        # 
        # Equation 11 -- source is emitted upward radiation at top of layer plus
        # radiation emitted at bottom of layer,
        # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
        src[ilev + 1, gcol] =
            src_up[ilev, gcol] + Tdif[ilev, gcol] * denom * (src[ilev, gcol] + albedo[ilev, gcol] * src_dn[ilev, gcol])
    end

    # Eq 12, at the top of the domain upwelling diffuse is due to ...
    @inbounds flux_up[nlev, gcol] =
        flux_dn[nlev, gcol] * albedo[nlev, gcol] + # ... reflection of incident diffuse and
        src[nlev, gcol]                          # scattering by the direct beam below

    # From the top of the atmosphere downward -- compute fluxes
    @inbounds for ilev in nlay:-1:1
        denom = FT(1) / (FT(1) - Rdif[ilev, gcol] * albedo[ilev, gcol])  # Eq 10
        flux_dn[ilev, gcol] =
            (Tdif[ilev, gcol] * flux_dn[ilev + 1, gcol] + # Equation 13
             Rdif[ilev, gcol] * src[ilev, gcol] +
             src_dn[ilev, gcol]) * denom
        flux_up[ilev, gcol] =
            flux_dn[ilev, gcol] * albedo[ilev, gcol] + # Equation 12
            src[ilev, gcol]
    end

    @inbounds for ilev in 1:nlev
        flux_dn[ilev, gcol] += flux_dn_dir[ilev, gcol]
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
end
# ---------------------------------------------------------------
