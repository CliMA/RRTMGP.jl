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
            # get column views
            flux_sw_col = FluxSWc(flux_sw, gcol)
            op_col = TwoStreamc(op, gcol)
            src_sw_col = SourceSW2Strc(src_sw, gcol)
            bcs_sw_col = SwBCsc(bcs_sw, gcol)
            #set_flux_to_zero!(flux_sw, gcol)
            set_flux_to_zero!(flux_sw_col)
            compute_optical_props!(op, as, gcol, igpt, nothing, nothing)
            # Cell properties: transmittance and reflectance for direct and diffuse radiation
            sw_two_stream!(op_col, src_sw_col, bcs_sw_col, nlay)
            # Direct-beam and source for diffuse radiation
            sw_source_2str!(src_sw_col, bcs_sw_col, flux_sw_col, igpt, ibnd, solar_frac, nlay)
            # Transport
            adding_sw!(src_sw_col, bcs_sw_col, flux_sw_col, igpt, n_gpt, ibnd, nlev)
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
    as::AbstractAtmosphericState{FT},
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt, ibnd = 1, 1, UInt8(1)
    solar_frac = FT(1)
    if gcol ≤ ncol
        @inbounds begin
            # get column views
            flux_sw_col = FluxSWc(flux_sw, gcol)
            op_col = TwoStreamc(op, gcol)
            src_sw_col = SourceSW2Strc(src_sw, gcol)
            bcs_sw_col = SwBCsc(bcs_sw, gcol)
            #set_flux_to_zero!(flux_sw, gcol)
            set_flux_to_zero!(flux_sw_col)
            compute_optical_props!(op, as, gcol, igpt, nothing, nothing)
            # Cell properties: transmittance and reflectance for direct and diffuse radiation
            sw_two_stream!(op_col, src_sw_col, bcs_sw_col, nlay)
            # Direct-beam and source for diffuse radiation
            sw_source_2str!(src_sw_col, bcs_sw_col, flux_sw_col, igpt, ibnd, solar_frac, nlay)
            # Transport
            adding_sw!(src_sw_col, bcs_sw_col, flux_sw_col, igpt, n_gpt, ibnd, nlev)
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
        bld_cld_mask = as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
        # setting references for flux_sw
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                # get column views
                flux_col = FluxSWc(flux, gcol)
                flux_sw_col = FluxSWc(flux_sw, gcol)
                op_col = TwoStreamc(op, gcol)
                src_sw_col = SourceSW2Strc(src_sw, gcol)
                bcs_sw_col = SwBCsc(bcs_sw, gcol)
                bld_cld_mask &&
                    Optics.build_cloud_mask!(as.cld_mask_sw, as.cld_frac, as.random_sw, gcol, igpt, as.cld_mask_type)
                #igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                igpt == 1 && set_flux_to_zero!(flux_sw_col)
                compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
                # Cell properties: transmittance and reflectance for direct and diffuse radiation
                sw_two_stream!(op_col, src_sw_col, bcs_sw_col, nlay)
                solar_frac = lookup_sw.solar_src_scaled[igpt]
                ibnd = lookup_sw.major_gpt2bnd[igpt]
                # Direct-beam and source for diffuse radiation
                sw_source_2str!(src_sw_col, bcs_sw_col, flux_col, igpt, ibnd, solar_frac, nlay)
                # Transport
                adding_sw!(src_sw_col, bcs_sw_col, flux_col, igpt, n_gpt, ibnd, nlev)
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
        # get column views
        flux_col = FluxSWc(flux, gcol)
        flux_sw_col = FluxSWc(flux_sw, gcol)
        op_col = TwoStreamc(op, gcol)
        src_sw_col = SourceSW2Strc(src_sw, gcol)
        bcs_sw_col = SwBCsc(bcs_sw, gcol)
        #set_flux_to_zero!(flux_sw, gcol)
        set_flux_to_zero!(flux_sw_col)
        @inbounds for igpt in 1:n_gpt
            if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
                Optics.build_cloud_mask!(as.cld_mask_sw, as.cld_frac, as.random_sw, gcol, igpt, as.cld_mask_type)
            end
            compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
            # Cell properties: transmittance and reflectance for direct and diffuse radiation
            sw_two_stream!(op_col, src_sw_col, bcs_sw_col, nlay)
            solar_frac = lookup_sw.solar_src_scaled[igpt]
            ibnd = lookup_sw.major_gpt2bnd[igpt]
            # Direct-beam and source for diffuse radiation
            sw_source_2str!(src_sw_col, bcs_sw_col, flux_col, igpt, ibnd, solar_frac, nlay)
            # Transport
            adding_sw!(src_sw_col, bcs_sw_col, flux_col, igpt, n_gpt, ibnd, nlev)
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
function sw_two_stream!(op::TwoStreamc, src_sw::SourceSW2Strc, bcs_sw::SwBCsc, nlay::Int)
    zenith = bcs_sw.zenith[1]
    FT = eltype(zenith)
    (; τ, ssa, g) = op
    (; Rdif, Tdif, Rdir, Tdir, Tnoscat) = src_sw
    k_min = FT(1e4 * eps(FT)) # Suggestion from Chiel van Heerwaarden
    μ₀ = FT(cos(zenith))

    @inbounds for glay in 1:nlay
        τ_glay, ssa_glay, g_glay = τ[glay], ssa[glay], g[glay]
        # Zdunkowski Practical Improved Flux Method "PIFM"
        #  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
        γ1 = (FT(8) - ssa_glay * (FT(5) + FT(3) * g_glay)) * FT(0.25)
        γ2 = FT(3) * (ssa_glay * (FT(1) - g_glay)) * FT(0.25)
        γ3 = (FT(2) - (FT(3) * μ₀) * g_glay) * FT(0.25)
        γ4 = FT(1) - γ3
        α1 = γ1 * γ4 + γ2 * γ3                          # Eq. 16
        α2 = γ1 * γ3 + γ2 * γ4                          # Eq. 17
        k = sqrt(max((γ1 - γ2) * (γ1 + γ2), k_min))

        exp_minusktau = exp(-τ_glay * k)
        exp_minus2ktau = exp_minusktau * exp_minusktau

        # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
        RT_term = FT(1) / (k * (FT(1) + exp_minus2ktau) + γ1 * (FT(1) - exp_minus2ktau))

        @inbounds Rdif[glay] = RT_term * γ2 * (FT(1) - exp_minus2ktau) # Eqn. 25
        @inbounds Tdif[glay] = RT_term * FT(2) * k * exp_minusktau     # Eqn. 26

        # Transmittance of direct, unscattered beam. Also used below
        @inbounds T₀ = Tnoscat[glay] = exp(-τ_glay / μ₀)

        # Direct reflect and transmission
        k_μ = k * μ₀
        k_γ3 = k * γ3
        k_γ4 = k * γ4

        # Equation 14, multiplying top and bottom by exp(-k*tau)
        #   and rearranging to avoid div by 0.
        if abs(FT(1) - k_μ * k_μ) ≥ eps(FT)
            RT_term = ssa_glay * RT_term / (FT(1) - k_μ * k_μ)
        else
            RT_term = ssa_glay * RT_term / eps(FT)
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
        Rdir[glay] = max(FT(0), min(Rdir_unconstrained, (FT(1) - T₀)))
        Tdir[glay] = max(FT(0), min(Tdir_unconstrained, (FT(1) - T₀ - Rdir[glay])))
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
    src_sw::SourceSW2Strc,
    bcs_sw::SwBCsc,
    flux::FluxSWc,
    igpt::Int,
    ibnd::UInt8,
    solar_frac::FT,
    nlay::Int,
) where {FT <: AbstractFloat}
    zenith = bcs_sw.zenith[1]
    toa_flux = bcs_sw.toa_flux[1]
    (; sfc_alb_direct) = bcs_sw
    (; Rdir, Tdir, Tnoscat, src_up, src_dn, sfc_source) = src_sw
    flux_dn_dir = flux.flux_dn_dir

    # layer index = level index
    # previous level is up (+1)
    @inbounds flux_dn_dir[nlay + 1] = toa_flux * solar_frac * cos(zenith)
    @inbounds for ilev in nlay:-1:1
        src_up[ilev] = Rdir[ilev] * flux_dn_dir[ilev + 1]
        src_dn[ilev] = Tdir[ilev] * flux_dn_dir[ilev + 1]
        flux_dn_dir[ilev] = Tnoscat[ilev] * flux_dn_dir[ilev + 1]
    end
    @inbounds sfc_source[1] = flux_dn_dir[1] * sfc_alb_direct[ibnd]
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
function adding_sw!(src_sw::SourceSW2Strc, bcs_sw::SwBCsc, flux::FluxSWc, igpt::Int, n_gpt::Int, ibnd::UInt8, nlev::Int)
    nlay = nlev - 1
    # destructuring
    (; flux_up, flux_dn, flux_dn_dir, flux_net) = flux
    (; albedo, Rdif, Tdif, src_up, src_dn, src) = src_sw
    sfc_source = src_sw.sfc_source[1]
    FT = eltype(sfc_source)
    @inbounds for ilev in 1:nlev
        flux_dn[ilev] = FT(0)
        flux_up[ilev] = FT(0)
    end

    # Albedo of lowest level is the surface albedo...
    @inbounds albedo[1] = bcs_sw.sfc_alb_diffuse[ibnd]
    # ... and source of diffuse radiation is surface emission
    @inbounds src[1] = sfc_source
    # From bottom to top of atmosphere --
    #   compute albedo and source of upward radiation
    @inbounds for ilev in 1:nlay
        denom = FT(1) / (FT(1) - Rdif[ilev] * albedo[ilev])  # Eq 10
        albedo[ilev + 1] = Rdif[ilev] + Tdif[ilev] * Tdif[ilev] * albedo[ilev] * denom # Equation 9
        # 
        # Equation 11 -- source is emitted upward radiation at top of layer plus
        # radiation emitted at bottom of layer,
        # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
        src[ilev + 1] = src_up[ilev] + Tdif[ilev] * denom * (src[ilev] + albedo[ilev] * src_dn[ilev])
    end

    # Eq 12, at the top of the domain upwelling diffuse is due to ...
    @inbounds flux_up[nlev] =
        flux_dn[nlev] * albedo[nlev] + # ... reflection of incident diffuse and
        src[nlev]                          # scattering by the direct beam below

    # From the top of the atmosphere downward -- compute fluxes
    @inbounds for ilev in nlay:-1:1
        denom = FT(1) / (FT(1) - Rdif[ilev] * albedo[ilev])  # Eq 10
        flux_dn[ilev] = (Tdif[ilev] * flux_dn[ilev + 1] + # Equation 13
                         Rdif[ilev] * src[ilev] +
                         src_dn[ilev]) * denom
        flux_up[ilev] =
            flux_dn[ilev] * albedo[ilev] + # Equation 12
            src[ilev]
    end

    @inbounds for ilev in 1:nlev
        flux_dn[ilev] += flux_dn_dir[ilev]
        flux_net[ilev] = flux_up[ilev] - flux_dn[ilev]
    end
end
# ---------------------------------------------------------------
