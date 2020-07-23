

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
) where {FT,I}
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
    Ds = @view angle_disc.gauss_Ds[1:n_μ, n_μ]
    w_μ = @view angle_disc.gauss_wts[1:n_μ, n_μ]

    lev_source_up = top_at_1 ? source.lev_source_dec : source.lev_source_inc
    lev_source_dn = top_at_1 ? source.lev_source_inc : source.lev_source_dec
    lay_source = source.lay_source
    sfc_source = source.sfc_source

    τ = op.τ
    τ_thresh = sqrt(eps(FT)) # or abs(eps(FT))?

    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    flux_dn_top = @view(flux_dn[:, i_lev_top])
    flux_dn_top .*= 1/(2 * FT(π) * w_μ[1])

    # Downward propagation
    @inbounds for ilev in lev_range(mesh_orientation), icol = 1:ncol
        τ_loc = τ[icol, ilev-b] * Ds[1] # τ at layer ilev-b
        trans = exp(-τ_loc)
        fact = τ_loc > τ_thresh ? (1 - trans) / τ_loc - trans :
            τ_loc * (FT(1 / 2) - 1 / 3 * τ_loc)
        source_dn =
            (1 - trans) * lev_source_dn[icol, ilev-b] +
            2 * fact * (lay_source[icol, ilev-b] - lev_source_dn[icol, ilev-b])
        flux_dn[icol, ilev] = trans * flux_dn[icol, ilev-n̂] + source_dn
    end

    # Surface reflection and emission
    for icol = 1:ncol
        flux_up[icol, n_sfc] =
            flux_dn[icol, n_sfc] * (1.0 - sfc_emis[icol]) +
            sfc_emis[icol] * sfc_source[icol]
    end

    # Upward propagation
    @inbounds for ilev in lev_range_reversed(mesh_orientation), icol = 1:ncol
        τ_loc = τ[icol, ilev-1+b] * Ds[1]
        trans = exp(-τ_loc)
        fact = τ_loc > τ_thresh ? (1 - trans) / τ_loc - trans :
            τ_loc * (FT(1 / 2) - 1 / 3 * τ_loc)
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
