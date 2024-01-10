"""
    rte_lw_solve!(
        context,
        flux::FluxLW{FT},
        flux_lw::FluxLW{FT},
        src_lw::SourceLW2Str{FT},
        bcs_lw::LwBCs{FT},
        op::TwoStream{FT},
        major_gpt2bnd::AbstractArray{UInt8,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    ) where {FT<:AbstractFloat}

Two stream solver for the longwave problem.

RRTMGP provides source functions at each level using the spectral mapping
of each adjacent layer. Combine these for two-stream calculations
lw combine sources (rte_lw_2stream_combine_sources!)
combine RRTMGP-specific sources at levels
compute layer reflectance, transmittance
compute total source function at levels using linear-in-tau (rte_lw_2stream_source!)
transport (adding_lw!)
"""
function rte_lw_solve!(
    context,
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLW2Str{FT},
    bcs_lw::LwBCs{FT},
    op::TwoStream{FT},
    major_gpt2bnd::AbstractArray{UInt8, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lookup_lw::Union{LookUpLW, Nothing} = nothing,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, flux_lw, src_lw, bcs_lw, op, nlay, ncol, major_gpt2bnd, as, lookup_lw, lookup_lw_cld)
        @cuda threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(args...)
    else
        bld_cld_mask = as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
        @inbounds begin
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
                    bld_cld_mask && Optics.build_cloud_mask!(
                        view(as.cld_mask_lw, :, gcol),
                        view(as.cld_frac, :, gcol),
                        as.cld_mask_type,
                    )
                    compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
                    rte_lw_2stream_combine_sources!(src_lw, gcol, nlev, ncol)
                    rte_lw_2stream_source!(op, src_lw, gcol, nlay, ncol)
                    adding_lw!(flux, src_lw, bcs_lw, gcol, igpt, major_gpt2bnd, nlev, ncol)
                    n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
                end
            end
        end
    end
    return nothing
end

function rte_lw_2stream_solve_CUDA!(
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLW2Str{FT},
    bcs_lw::LwBCs{FT},
    op::TwoStream{FT},
    nlay,
    ncol,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    as::AbstractAtmosphericState{FT},
    lookup_lw::Union{LookUpLW, Nothing} = nothing,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        @inbounds for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
            if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
                Optics.build_cloud_mask!(view(as.cld_mask_lw, :, gcol), view(as.cld_frac, :, gcol), as.cld_mask_type)
            end
            compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
            rte_lw_2stream_combine_sources!(src_lw, gcol, nlev, ncol)
            rte_lw_2stream_source!(op, src_lw, gcol, nlay, ncol)
            adding_lw!(flux, src_lw, bcs_lw, gcol, igpt, major_gpt2bnd, nlev, ncol)
            n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
        end
    end
    return nothing
end


"""
    rte_lw_2stream_combine_sources!(src::SourceLW2Str, gcol, nlev, ncol)

RRTMGP provides source functions at each level using the spectral mapping
of each adjacent layer. This function combines these for two-stream calculations.
"""
function rte_lw_2stream_combine_sources!(src::SourceLW2Str, gcol, nlev, ncol)
    @inbounds src.lev_source[1, gcol] = src.lev_source_dec[1, gcol] # glev == 1
    @inbounds src.lev_source[nlev, gcol] = src.lev_source_inc[nlev - 1, gcol] # glev == nlev
    @inbounds for glev in 2:(nlev - 1)
        src.lev_source[glev, gcol] = sqrt(src.lev_source_dec[glev, gcol] * src.lev_source_inc[glev - 1, gcol])
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
) where {FT <: AbstractFloat}
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
    (; τ, ssa, g) = op
    (; Rdif, Tdif, lev_source, src_up, src_dn) = src_lw
    #k_min = FT === Float64 ? FT(1e-12) : FT(1e-4) used in RRTMGP-RTE FORTRAN code
    k_min = sqrt(eps(FT)) #FT(1e4 * eps(FT))
    lw_diff_sec = FT(1.66)
    τ_thresh = 100 * eps(FT)# tau(icol,ilay) > 1.0e-8_wp used in rte-rrtmgp
    # this is chosen to prevent catastrophic cancellation in src_up and src_dn calculation

    @inbounds for glay in 1:nlay
        γ1 = lw_diff_sec * (1 - FT(0.5) * ssa[glay, gcol] * (1 + g[glay, gcol]))
        γ2 = lw_diff_sec * FT(0.5) * ssa[glay, gcol] * (1 - g[glay, gcol])
        k = sqrt(max((γ1 + γ2) * (γ1 - γ2), k_min))
        τ_lay = τ[glay, gcol]

        coeff = exp(-2 * τ_lay * k)
        # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
        RT_term = 1 / (k * (1 + coeff) + γ1 * (1 - coeff))

        @inbounds Rdif[glay, gcol] = RT_term * γ2 * (1 - coeff) # Equation 25
        @inbounds Tdif[glay, gcol] = RT_term * 2 * k * exp(-τ_lay * k) # Equation 26

        # Source function for diffuse radiation
        # Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
        # This version straight from ECRAD
        # Source is provided as W/m2-str; factor of pi converts to flux units
        # lw_source_2str

        lev_src_bot = lev_source[glay, gcol]
        lev_src_top = lev_source[glay + 1, gcol]
        Rdif_lay = Rdif[glay, gcol]
        Tdif_lay = Tdif[glay, gcol]

        if τ_lay > τ_thresh
            # Toon et al. (JGR 1989) Eqs 26-27
            Z = (lev_src_bot - lev_src_top) / (τ_lay * (γ1 + γ2))
            Zup_top = Z + lev_src_top
            Zup_bottom = Z + lev_src_bot
            Zdn_top = -Z + lev_src_top
            Zdn_bottom = -Z + lev_src_bot
            @inbounds src_up[glay, gcol] = pi * (Zup_top - Rdif_lay * Zdn_top - Tdif_lay * Zup_bottom)
            @inbounds src_dn[glay, gcol] = pi * (Zdn_bottom - Rdif_lay * Zup_bottom - Tdif_lay * Zdn_top)
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
    gcol::Int,
    igpt::Int,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    nlev::Int,
    ncol::Int,
) where {FT <: AbstractFloat, SL <: SourceLW2Str{FT}, BCL <: LwBCs{FT}}
    nlay = nlev - 1
    ibnd = major_gpt2bnd[igpt]
    # setting references
    (; flux_up, flux_dn, flux_net) = flux

    (; albedo, sfc_source, Rdif, Tdif, src_up, src_dn, src) = src_lw
    (; inc_flux, sfc_emis) = bcs_lw

    @inbounds for ilev in 1:nlev
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
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
end
