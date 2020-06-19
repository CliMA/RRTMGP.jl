

struct GrayRTE{FT}
    fluxes::FluxesBroadBand{FT}
    hr_lay::Array{FT,2}

    function GrayRTE(
        ncol::Int,
        nlay::Int,
        nlev::Int,
        ::Type{FT},
    ) where {FT<:AbstractFloat}
        return new{FT}(FluxesBroadBand(FT, (ncol, nlev)), zeros(FT, ncol, nlay))
    end
end

function rte_lw_gray_solve!(
    fluxes::FluxesBroadBand{FT},
    op::AbstractOpticalPropsArry{FT,I},
    mesh_orientation::MeshOrientation{I},
    bcs::LongwaveBCs{FT},
    source::SourceFuncLongWave{FT,I},
    angle_disc::Union{GaussQuadrature{FT,I},Nothing} = nothing,
) where {FT<:AbstractFloat,I<:Int}
    # Upper boundary condition
    sfc_emis = bcs.sfc_emis
    flux_up, flux_dn, flux_net = fluxes.flux_up, fluxes.flux_dn, fluxes.flux_net
    inc_flux = bcs.inc_flux
    i_lev_top = ilev_top(mesh_orientation)
    top_at_1 = mesh_orientation.top_at_1
    n̂ = nhat(mesh_orientation)
    b = binary(mesh_orientation)
    n_sfc = ilev_bot(mesh_orientation)

    fill!(flux_dn, 0)
    fill!(flux_up, 0)
    if inc_flux ≠ nothing
        flux_dn[1, i_lev_top] = inc_flux[1]
    end

    ncol = size(flux_up, 1)
    nlev = size(flux_up, 2)
    nlay = nlev - 1
    ngpt = 1

    n_μ = angle_disc.n_gauss_angles
    Ds = angle_disc.gauss_Ds
    w_μ = angle_disc.gauss_wts

    lev_source_up = top_at_1 ? source.lev_source_dec : source.lev_source_inc
    lev_source_dn = top_at_1 ? source.lev_source_inc : source.lev_source_dec
    lay_source = source.lay_source
    sfc_source = source.sfc_source

    τ = op.τ
    τ_thresh = sqrt(eps(FT)) # or abs(eps(FT))?

    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    for icol = 1:ncol
        flux_dn[icol, i_lev_top] /= 2 * FT(π) * w_μ[1]
    end

    # Downward propagation
    @inbounds for ilev in lev_range(mesh_orientation), icol = 1:ncol
        τ_loc = τ[icol, ilev-b] * Ds[1] # τ at layer ilev-b
        trans = exp(-τ_loc)
        fact = τ_loc > τ_thresh ? (1 - trans) / τ_loc - trans :
            τ_loc * (FT(1 / 2) - FT(1 / 3) * τ_loc)
        source_dn =
            (1 - trans) * lev_source_dn[icol, ilev-b] +
            2 * fact * (lay_source[icol, ilev-b] - lev_source_dn[icol, ilev-b])
        flux_dn[icol, ilev] = trans * flux_dn[icol, ilev-n̂] + source_dn
    end

    # Surface reflection and emission
    for icol = 1:ncol
        flux_up[icol, n_sfc] =
            flux_dn[icol, n_sfc] * (1 - sfc_emis[icol]) +
            sfc_emis[icol] * sfc_source[icol]
    end

    # Upward propagation
    @inbounds for ilev in lev_range_reversed(mesh_orientation), icol = 1:ncol
        τ_loc = τ[icol, ilev-1+b] * Ds[1]
        trans = exp(-τ_loc)
        fact = τ_loc > τ_thresh ? (1 - trans) / τ_loc - trans :
            τ_loc * (FT(1 / 2) - FT(1 / 3) * τ_loc)
        source_up =
            (1 - trans) * lev_source_up[icol, ilev-1+b] +
            2 *
            fact *
            (lay_source[icol, ilev-1+b] - lev_source_up[icol, ilev-1+b])
        flux_up[icol, ilev] = trans * flux_up[icol, ilev+n̂] + source_up
    end

    # Convert intensity to flux assuming azimuthal isotropy and quadrature weight
    flux_up .*= 2 * FT(π) * w_μ[1]
    flux_dn .*= 2 * FT(π) * w_μ[1]
    flux_net .= flux_up .- flux_dn

    return nothing
end

function rte_lw_2stream_gray_solve!(
    fluxes::FluxesBroadBand{FT},
    op::AbstractOpticalPropsArry{FT,I},
    mesh_orientation::MeshOrientation{I},
    bcs::LongwaveBCs{FT},
    source::SourceFuncLongWave{FT,I},
    angle_disc::Union{GaussQuadrature{FT,I},Nothing} = nothing,
) where {FT<:AbstractFloat,I<:Int}
    sfc_emis = bcs.sfc_emis
    flux_up, flux_dn, flux_net = fluxes.flux_up, fluxes.flux_dn, fluxes.flux_net
    lev_src_inc, lev_src_dec = source.lev_source_inc, source.lev_source_dec
    lay_source, sfc_source = source.lay_source, source.sfc_source
    inc_flux = bcs.inc_flux
    i_lev_top = ilev_top(mesh_orientation)
    top_at_1 = mesh_orientation.top_at_1
    n̂ = nhat(mesh_orientation)
    b = binary(mesh_orientation)
    n_sfc = ilev_bot(mesh_orientation)

    τ = op.τ
    g = op.g
    ssa = op.ssa

    fill!(flux_dn, 0)
    fill!(flux_up, 0)
    if inc_flux ≠ nothing
        flux_dn[1, i_lev_top] = inc_flux[1]
    end

    ncol = size(flux_up, 1)
    nlev = size(flux_up, 2)
    nlay = nlev - 1
    ngpt = 1
    #-----------------------------------------------------------------------------------
    lev_source = Array{FT}(undef, ncol, nlev)

    γ1 = Array{FT}(undef, ncol, nlay)
    γ2 = similar(γ1)
    Rdif = similar(γ1)
    Tdif = similar(γ1)
    lw_diff_sec = FT(1.66)
    #-----------------------------------------------------------------------------------

    for icol = 1:ncol
        # lw combine sources
        lev_source[icol, 1] = lev_src_dec[icol, 1]
        lev_source[icol, nlev] = lev_src_inc[icol, nlay]

        for ilay = 2:nlay
            lev_source[icol, ilay] =
                sqrt(lev_src_dec[icol, ilay] * lev_src_inc[icol, ilay-1])
        end
        #-----------------------------------------------
    end
    # Cell properties: reflection, transmission for diffuse radiation
    # Coupling coefficients needed for source function
    for icol = 1:ncol
        for ilay = 1:nlay
            γ1[icol, ilay] =
                lw_diff_sec *
                (1.0 - 0.5 * ssa[icol, ilay] * (1.0 + g[icol, ilay]))
            γ2[icol, ilay] =
                lw_diff_sec * 0.5 * ssa[icol, ilay] * (1.0 - g[icol, ilay])
            k = sqrt(max(
                (γ1[icol, ilay] + γ2[icol, ilay]) *
                (γ1[icol, ilay] - γ2[icol, ilay]),
                1e-12,
            ))

            # Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
            RT_term =
                1.0 / (
                    k * (1.0 + exp(-2.0 * τ[icol, ilay] * k)) +
                    γ1[icol, ilay] * (1.0 - exp(-2 * τ[icol, ilay] * k))
                )

            Rdif[icol, ilay] =
                RT_term * γ2[icol, ilay] * (1.0 - exp(-2 * τ[icol, ilay] * k)) # Equation 25
            Tdif[icol, ilay] = RT_term * 2.0 * k * exp(-τ[icol, ilay] * k) # Equation 26
        end
    end

    # Source function for diffuse radiation
    # lw_source_2str
    source_up = Array{FT}(undef, ncol, nlay)
    source_dn = Array{FT}(undef, ncol, nlay)

    for ilay = 1:nlay
        if top_at_1
            top = ilay
            bot = ilay + 1
        else
            top = ilay + 1
            bot = ilay
        end

        for icol = 1:ncol
            if τ[icol, ilay] > 1.0e-8
                # Toon et al. (JGR 1989) Eqs 26-27
                Z =
                    (lev_source[icol, bot] - lev_source[icol, top]) /
                    (τ[icol, ilay] * (γ1[icol, ilay] + γ2[icol, ilay]))
                Zup_top = Z + lev_source[icol, top]
                Zup_bottom = Z + lev_source[icol, bot]
                Zdn_top = -Z + lev_source[icol, top]
                Zdn_bottom = -Z + lev_source[icol, bot]
                source_up[icol, ilay] =
                    pi * (
                        Zup_top - Rdif[icol, ilay] * Zdn_top -
                        Tdif[icol, ilay] * Zup_bottom
                    )
                source_dn[icol, ilay] =
                    pi * (
                        Zdn_bottom - Rdif[icol, ilay] * Zup_bottom -
                        Tdif[icol, ilay] * Zdn_top
                    )
            else
                source_up[icol, ilay] = 0.0
                source_dn[icol, ilay] = 0.0
            end
        end

    end
    # Transport
    adding_gray!(
        top_at_1,
        sfc_emis,
        Rdif,
        Tdif,
        source_up,
        source_dn,
        sfc_source,
        flux_up,
        flux_dn,
    )
    flux_net .= flux_up .- flux_dn

    return nothing
end

function adding_gray!(
    top_at_1::Bool,
    sfc_emis::Array{FT,2},
    Rdif::Array{FT,2},
    Tdif::Array{FT,2},
    src_up::Array{FT,2},
    src_dn::Array{FT,2},
    sfc_source::Array{FT,2},
    flux_up::Array{FT,2},
    flux_dn::Array{FT,2},
) where {FT,I}
    ncol = size(flux_up, 1)
    nlev = size(flux_up, 2)
    nlay = nlev - 1

    albedo = zeros(FT, nlev)
    src = zeros(FT, nlev)
    denom = zeros(FT, nlev)
    if top_at_1
        #--------------------------------------------------------------------------------------
        for icol = 1:ncol
            # Albedo of lowest level is the surface albedo...
            albedo[nlev] = 1.0 - sfc_emis[icol]
            # ... and source of diffuse radiation is surface emission
            src[nlev] = π * sfc_emis[icol] * sfc_source[icol]
            # From bottom to top of atmosphere --
            # compute albedo and source of upward radiation
            for ilev = nlay:-1:1
                denom[ilev] = 1.0 / (1.0 - Rdif[icol, ilev] * albedo[ilev+1])
                albedo[ilev] =
                    Rdif[icol, ilev] +
                    Tdif[ilev] * Tdif[ilev] * albedo[ilev+1] * denom[ilev]
                # Equation 11 -- source is emitted upward radiation at top of layer plus
                # radiation emitted at bottom of layer,
                # transmitted through the layer and reflected from layers below (tdiff*src*albedo)
                src[ilev] =
                    src_up[icol, ilev] +
                    Tdif[icol, ilev] *
                    denom[ilev] *
                    (src[ilev+1] + albedo[ilev+1] * src_dn[icol, ilev])
            end

            # Eq 12, at the top of the domain upwelling diffuse is due to ...
            flux_up[icol, 1] = flux_dn[icol, 1] * albedo[1] + src[1] # reflection of incident diffuse and
            # emission from below

            # From the top of the atmosphere downward -- compute fluxes
            for ilev = 2:nlay+1
                flux_dn[icol, ilev] =
                    (
                        Tdif[icol, ilev-1] * flux_dn[icol, ilev-1] + # Equation 13
                        Rdif[icol, ilev-1] * src[ilev] +
                        src_dn[icol, ilev-1]
                    ) * denom[ilev-1]
                flux_up[icol, ilev] =
                    flux_dn[icol, ilev] * albedo[ilev] + # Equation 12
                    src[ilev]
            end

        end
        #--------------------------------------------------------------------------------------
    else
        #--------------------------------------------------------------------------------------
        for icol = 1:ncol
            # Albedo of lowest level is the surface albedo...
            albedo[1] = albedo_sfc[icol]
            # ... and source of diffuse radiation is surface emission
            src[1] = src_sfc[icol]
            #--------------------------------------------------------
            # From bottom to top of atmosphere --
            #   compute albedo and source of upward radiation
            for ilev = 1:nlay
                denom[ilev] = 1.0 / (1.0 - Rdif[icol, ilev] * albedo[ilev])  # Eq 10
                albedo[ilev+1] =
                    Rdif[icol, ilev] +
                    Tdif[icol, ilev] *
                    Tdif[icol, ilev] *
                    albedo[ilev] *
                    denom[ilev] # Equation 9
                # 
                # Equation 11 -- source is emitted upward radiation at top of layer plus
                # radiation emitted at bottom of layer,
                # transmitted through the layer and reflected from layers below (Tdiff*src*albedo)
                src[ilev+1] =
                    src_up[icol, ilev] +
                    Tdif[icol, ilev] *
                    denom[ilev] *
                    (src[ilev] + albedo[ilev] * src_dn[icol, ilev])
            end

            # Eq 12, at the top of the domain upwelling diffuse is due to ...
            flux_up[icol, nlev] =
                flux_dn[icol, nlev] * albedo[nlev] + # ... reflection of incident diffuse and
                src[nlev]                          # scattering by the direct beam below

            # From the top of the atmosphere downward -- compute fluxes
            for ilev = nlay:-1:1
                flux_dn[icol, ilev] =
                    (
                        Tdif[icol, ilev] * flux_dn[icol, ilev+1] + # Equation 13
                        Rdif[icol, ilev] * src[ilev] +
                        src_dn[icol, ilev]
                    ) * denom[ilev]
                flux_up[icol, ilev] =
                    flux_dn[icol, ilev] * albedo[ilev] + # Equation 12
                    src[ilev]

            end
            #--------------------------------------------------------


        end
        #--------------------------------------------------------------------------------------
    end

    return nothing
end
