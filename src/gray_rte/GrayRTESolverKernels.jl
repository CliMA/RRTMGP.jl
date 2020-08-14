

@kernel function rte_lw_noscat_gray_source_kernel!(
    lay_source::FTA3D,
    lev_source_up::FTA3D,
    lev_source_dn::FTA3D,
    τ::FTA2D,
    trans::FTA2D,
    source_up::FTA2D,
    source_dn::FTA2D,
    Ds::FTA1D,
    ::Val{ncol},
    ::Val{nlay},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    ncol,
    nlay,
}
    gcol, glay = @index(Global, NTuple)  # global col & lay ids

    τ_thresh = sqrt(eps(FT)) # or abs(eps(FT))?
    τ_loc = τ[gcol, glay] * Ds[1]
    trans[gcol, glay] = exp(-τ_loc)
    fact = (τ_loc > τ_thresh) ?
        ((FT(1) - trans[gcol, glay]) / τ_loc - trans[gcol, glay]) :
        τ_loc * (FT(1 / 2) - FT(1 / 3) * τ_loc)

    source_up[gcol, glay] =
        (FT(1) - trans[gcol, glay]) * lev_source_up[gcol, glay, 1] +
        FT(2) *
        fact *
        (lay_source[gcol, glay, 1] - lev_source_up[gcol, glay, 1])

    source_dn[gcol, glay] =
        (FT(1) - trans[gcol, glay]) * lev_source_dn[gcol, glay, 1] +
        FT(2) *
        fact *
        (lay_source[gcol, glay, 1] - lev_source_dn[gcol, glay, 1])

    @synchronize
end

@kernel function rte_lw_noscat_gray_transport_kernel!(
    flux_up::FTA2D,
    flux_dn::FTA2D,
    flux_net::FTA2D,
    trans::FTA2D,
    source_up::FTA2D,
    source_dn::FTA2D,
    sfc_source::FTA2D,
    sfc_emis::FTA1D,
    inc_flux::Union{Nothing,FTA1D},
    τ::FTA2D,
    n_μ::Int,
    Ds::FTA1D,
    w_μ::FTA1D,
    top_at_1::Bool,
    ::Val{ncol},
    ::Val{nlev},
    ::Val{nlay},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    ncol,
    nlev,
    nlay,
}
    gcol = @index(Global, Linear) # global col ids
    i_lev_top = top_at_1 ? 1 : nlev # index for top level of column

    for ilev = 1:nlev
        flux_dn[gcol, ilev] = FT(0)
        flux_up[gcol, ilev] = FT(0)
    end

    if inc_flux ≠ nothing
        flux_dn[gcol, i_lev_top] = inc_flux[gcol]
    end
    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    flux_dn[gcol, i_lev_top] /= (FT(2) * FT(π) * w_μ[1])

    if top_at_1
        #-------Top of domain is index 1
        # Downward propagation
        for ilev = 2:nlev
            flux_dn[gcol, ilev] =
                trans[gcol, ilev-1] * flux_dn[gcol, ilev-1] +
                source_dn[gcol, ilev-1]
        end
        # Surface reflection and emission
        flux_up[gcol, nlev] =
            flux_dn[gcol, nlev] * (FT(1) - sfc_emis[gcol]) + sfc_source[gcol, 1]

        # Upward propagation
        for ilev = nlay:-1:1
            flux_up[gcol, ilev] =
                trans[gcol, ilev] * flux_up[gcol, ilev+1] +
                source_up[gcol, ilev]
        end
        #-----------------------------------
    else
        #-------Top of domain is index nlev
        # Downward propagation
        for ilev = nlay:-1:1
            flux_dn[gcol, ilev] =
                trans[gcol, ilev] * flux_dn[gcol, ilev+1] +
                source_dn[gcol, ilev]
        end

        # Surface reflection and emission
        flux_up[gcol, 1] =
            flux_dn[gcol, 1] * (FT(1) - sfc_emis[gcol]) + sfc_source[gcol]

        # Upward propagation
        for ilev = 2:nlay+1
            flux_up[gcol, ilev] =
                trans[gcol, ilev-1] * flux_up[gcol, ilev-1] +
                source_up[gcol, ilev-1]
        end
        #-----------------------------------
    end

    for ilev = 1:nlev
        flux_up[gcol, ilev] *= (FT(2) * FT(π) * w_μ[1])
        flux_dn[gcol, ilev] *= (FT(2) * FT(π) * w_μ[1])
        flux_net[gcol, ilev] = flux_up[gcol, ilev] - flux_dn[gcol, ilev]
    end
    @synchronize
end
