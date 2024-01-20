

function rte_lw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    max_threads,
    as::GrayAtmosphericState,
)
    (; nlay, ncol) = as
    nlev = nlay + 1
    igpt, ibnd = 1, UInt8(1)
    (; flux_up, flux_dn, flux_net) = flux_lw
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            compute_optical_props!(op, as, src_lw, gcol)
            rte_lw_2stream!(op, flux_lw, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
            for ilev in 1:nlev
                flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
            end
        end
    end
    return nothing
end

function rte_lw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    max_threads,
    as::GrayAtmosphericState,
)
    (; nlay, ncol) = as
    nlev = nlay + 1
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux_lw, src_lw, bcs_lw, op, nlay, ncol, as)
    @cuda threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(args...)
    return nothing
end

function rte_lw_2stream_solve_CUDA!(
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    nlay,
    ncol,
    as::GrayAtmosphericState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    igpt, ibnd = 1, UInt8(1)
    if gcol ≤ ncol
        (; flux_up, flux_dn, flux_net) = flux_lw
        compute_optical_props!(op, as, src_lw, gcol)
        rte_lw_2stream!(op, flux_lw, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
        @inbounds begin
            for ilev in 1:nlev
                flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
            end
        end
    end
    return nothing
end

function rte_lw_2stream_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    max_threads,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    (; nlay, ncol) = as
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw
    n_gpt = length(major_gpt2bnd)
    bld_cld_mask = as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
    flux_up_lw = flux_lw.flux_up
    flux_dn_lw = flux_lw.flux_dn
    flux_net_lw = flux_lw.flux_net
    (; flux_up, flux_dn) = flux
    @inbounds begin
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                ibnd = major_gpt2bnd[igpt]
                bld_cld_mask && Optics.build_cloud_mask!(
                    view(as.cld_mask_lw, :, gcol),
                    view(as.cld_frac, :, gcol),
                    as.cld_mask_type,
                )
                compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
                rte_lw_2stream!(op, flux, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
                if igpt == 1
                    map!(x -> x, view(flux_up_lw, :, gcol), view(flux_up, :, gcol))
                    map!(x -> x, view(flux_dn_lw, :, gcol), view(flux_dn, :, gcol))
                else
                    for ilev in 1:nlev
                        @inbounds flux_up_lw[ilev, gcol] += flux_up[ilev, gcol]
                        @inbounds flux_dn_lw[ilev, gcol] += flux_dn[ilev, gcol]
                    end
                end
            end
        end
        ClimaComms.@threaded device for gcol in 1:ncol
            for ilev in 1:nlev
                flux_net_lw[ilev, gcol] = flux_up_lw[ilev, gcol] - flux_dn_lw[ilev, gcol]
            end
        end
    end
    return nothing
end

function rte_lw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    max_threads,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    (; nlay, ncol) = as
    nlev = nlay + 1
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux, flux_lw, src_lw, bcs_lw, op, nlay, ncol, as, lookup_lw, lookup_lw_cld)
    @cuda threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(args...)
    return nothing
end

function rte_lw_2stream_solve_CUDA!(
    flux::FluxLW,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    nlay,
    ncol,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        flux_up_lw = flux_lw.flux_up
        flux_dn_lw = flux_lw.flux_dn
        flux_net_lw = flux_lw.flux_net
        (; flux_up, flux_dn) = flux
        @inbounds for igpt in 1:n_gpt
            ibnd = major_gpt2bnd[igpt]
            if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
                Optics.build_cloud_mask!(view(as.cld_mask_lw, :, gcol), view(as.cld_frac, :, gcol), as.cld_mask_type)
            end
            compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
            rte_lw_2stream!(op, flux, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
            if igpt == 1
                map!(x -> x, view(flux_up_lw, :, gcol), view(flux_up, :, gcol))
                map!(x -> x, view(flux_dn_lw, :, gcol), view(flux_dn, :, gcol))
            else
                for ilev in 1:nlev
                    @inbounds flux_up_lw[ilev, gcol] += flux_up[ilev, gcol]
                    @inbounds flux_dn_lw[ilev, gcol] += flux_dn[ilev, gcol]
                end
            end
        end
        @inbounds begin
            for ilev in 1:nlev
                flux_net_lw[ilev, gcol] = flux_up_lw[ilev, gcol] - flux_dn_lw[ilev, gcol]
            end
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
    #k_min = FT === Float64 ? FT(1e-12) : FT(1e-4) used in RRTMGP-RTE FORTRAN code
    k_min = sqrt(eps(FT)) #FT(1e4 * eps(FT))
    lw_diff_sec = FT(1.66)
    τ_thresh = 100 * eps(FT)# tau(icol,ilay) > 1.0e-8_wp used in rte-rrtmgp
    # this is chosen to prevent catastrophic cancellation in src_up and src_dn calculation

    γ1 = lw_diff_sec * (1 - FT(0.5) * ssa * (1 + g))
    γ2 = lw_diff_sec * FT(0.5) * ssa * (1 - g)
    k = sqrt(max((γ1 + γ2) * (γ1 - γ2), k_min))

    coeff = exp(-2 * τ * k)
    # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    RT_term = 1 / (k * (1 + coeff) + γ1 * (1 - coeff))

    Rdif = RT_term * γ2 * (1 - coeff) # Equation 25
    Tdif = RT_term * 2 * k * exp(-τ * k) # Equation 26

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
        ibnd::UInt8,
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
    ibnd::UInt8,
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
    @inbounds flux_dn_ilevplus1 = isnothing(inc_flux) ? FT(0) : inc_flux[gcol, igpt]
    @inbounds flux_dn[nlev, gcol] = flux_dn_ilevplus1
    # Albedo of lowest level is the surface albedo...
    @inbounds albedo_ilev = FT(1) - sfc_emis[ibnd, gcol]
    @inbounds albedo[1, gcol] = albedo_ilev
    # ... and source of diffuse radiation is surface emission
    @inbounds src[1, gcol] = FT(π) * sfc_emis[ibnd, gcol] * sfc_source[gcol]

    # From bottom to top of atmosphere --
    #   compute albedo and source of upward radiation
    @inbounds lev_src_bot = lev_source[1, gcol]
    @inbounds for ilev in 1:nlay
        lev_src_top = lev_source[ilev + 1, gcol]
        τ_lay, ssa_lay, g_lay = τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        Rdif, Tdif, src_up, src_dn = lw_2stream_coeffs(τ_lay, ssa_lay, g_lay, lev_src_bot, lev_src_top)
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)  # Eq 10
        albedo_ilevplus1 = Rdif + Tdif * Tdif * albedo_ilev * denom # Equation 9
        albedo[ilev + 1, gcol] = albedo_ilevplus1
        # 
        # Equation 11 -- source is emitted upward radiation at top of layer plus
        # radiation emitted at bottom of layer,
        # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
        src[ilev + 1, gcol] = src_up + Tdif * denom * (src[ilev, gcol] + albedo_ilev * src_dn)
        lev_src_bot = lev_src_top
        albedo_ilev = albedo_ilevplus1
    end

    # Eq 12, at the top of the domain upwelling diffuse is due to ...
    @inbounds flux_up[nlev, gcol] =
        flux_dn_ilevplus1 * albedo[nlev, gcol] + # ... reflection of incident diffuse and
        src[nlev, gcol]                          # scattering by the direct beam below

    # From the top of the atmosphere downward -- compute fluxes
    @inbounds lev_src_top = lev_source[nlay + 1, gcol]
    @inbounds for ilev in nlay:-1:1
        lev_src_bot = lev_source[ilev, gcol]
        src_ilev = src[ilev, gcol]
        τ_lay, ssa_lay, g_lay = τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        Rdif, Tdif, _, src_dn = lw_2stream_coeffs(τ_lay, ssa_lay, g_lay, lev_src_bot, lev_src_top)

        denom = FT(1) / (FT(1) - Rdif * albedo[ilev, gcol])  # Eq 10
        #flux_dn[ilev, gcol] = (Tdif * flux_dn[ilev + 1, gcol] + # Equation 13
        flux_dn_ilev = (Tdif * flux_dn_ilevplus1 + # Equation 13
                        Rdif * src_ilev +
                        src_dn) * denom
        flux_up[ilev, gcol] =
            flux_dn_ilev * albedo[ilev, gcol] + # Equation 12
            src_ilev
        flux_dn[ilev, gcol] = flux_dn_ilev
        flux_dn_ilevplus1 = flux_dn_ilev
        lev_src_top = lev_src_bot
    end
end
