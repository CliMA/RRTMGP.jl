"""
    rte_lw_solve!(
        context,
        flux::FluxLW{FT},
        flux_lw::FluxLW{FT},
        src_lw::SourceLWNoScat{FT},
        bcs_lw::LwBCs{FT},
        op::OneScalar{FT},
        major_gpt2bnd::AbstractArray{UInt8,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    ) where {FT<:AbstractFloat}

No scattering solver for the longwave problem.
LW fluxes, no scattering, mu (cosine of integration angle) specified by column
Does radiation calculation at user-supplied angles; converts radiances to flux
using user-supplied weights. Currently, only n_gauss_angles = 1 supported.
"""
function rte_lw_solve!(
    context,
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
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
        @cuda threads = (tx) blocks = (bx) rte_lw_noscat_solve_CUDA!(args...)
    else
        @inbounds begin
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
                    compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
                    rte_lw_noscat_source!(src_lw, op, gcol, nlay)
                    rte_lw_noscat_transport!(src_lw, bcs_lw, op, gcol, flux, igpt, major_gpt2bnd, nlay, nlev)
                    n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
                end
            end
        end
    end
    return nothing
end

function rte_lw_noscat_solve_CUDA!(
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
    nlay,
    ncol,
    major_gpt2bnd,
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
            compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
            rte_lw_noscat_source!(src_lw, op, gcol, nlay)
            rte_lw_noscat_transport!(src_lw, bcs_lw, op, gcol, flux, igpt, major_gpt2bnd, nlay, nlev)
            n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
        end
    end
    return nothing
end

"""
    rte_lw_noscat_source!(
        src_lw::SourceLWNoScat{FT},
        op::OneScalar{FT},
        gcol,
        nlay,
    ) where {FT<:AbstractFloat}

Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
"""
function rte_lw_noscat_source!(src_lw::SourceLWNoScat{FT}, op::OneScalar{FT}, gcol, nlay) where {FT <: AbstractFloat}
    # setting references
    (; src_up, src_dn, lev_source_inc, lev_source_dec, lay_source) = src_lw
    Ds = op.angle_disc.gauss_Ds
    τ = op.τ

    τ_thresh = sqrt(sqrt(eps(FT))) # or abs(eps(FT))?

    @inbounds for glay in 1:nlay
        τ_loc = τ[glay, gcol] * Ds[1]      # Optical path and transmission,

        trans = exp(-τ_loc)    # used in source function and transport calculations
        # Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
        # is of order epsilon (smallest difference from 1. in working precision)
        # Thanks to Peter Blossey
        # Updated to 3rd order series and lower threshold based on suggestion from Dmitry Alexeev (Nvidia)

        fact =
            (τ_loc > τ_thresh) ? ((FT(1) - trans) / τ_loc - trans) :
            τ_loc * (FT(1 / 2) + τ_loc * (-FT(1 / 3) + τ_loc * FT(1 / 8)))

        # Equations below are developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
        src_up[glay, gcol] =
            (FT(1) - trans) * lev_source_inc[glay, gcol] +
            FT(2) * fact * (lay_source[glay, gcol] - lev_source_inc[glay, gcol])

        src_dn[glay, gcol] =
            (FT(1) - trans) * lev_source_dec[glay, gcol] +
            FT(2) * fact * (lay_source[glay, gcol] - lev_source_dec[glay, gcol])
    end
end

"""
    rte_lw_noscat_transport!(
        src_lw::SourceLWNoScat{FT},
        bcs_lw::LwBCs{FT},
        op::OneScalar{FT},
        gcol,
        flux::FluxLW{FT},
        igpt,
        ibnd,
        nlay,
        nlev,
    ) where {FT<:AbstractFloat}

Transport for no-scattering longwave problem.
"""
function rte_lw_noscat_transport!(
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
    gcol,
    flux::FluxLW{FT},
    igpt,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    nlay,
    nlev,
) where {FT <: AbstractFloat}
    ibnd = major_gpt2bnd[igpt]
    # setting references
    (; src_up, src_dn, sfc_source) = src_lw
    (; sfc_emis, inc_flux) = bcs_lw
    (; flux_up, flux_dn, flux_net) = flux

    Ds = op.angle_disc.gauss_Ds
    w_μ = op.angle_disc.gauss_wts
    n_μ = op.angle_disc.n_gauss_angles
    τ = op.τ

    @inbounds for ilev in 1:nlev
        flux_dn[ilev, gcol] = FT(0)
        flux_up[ilev, gcol] = FT(0)
    end

    if inc_flux ≠ nothing
        @inbounds flux_dn[nlev, gcol] = inc_flux[gcol, igpt]
    end
    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    @inbounds flux_dn[nlev, gcol] /= (FT(2) * FT(π) * w_μ[1])

    # Top of domain is index nlev
    # Downward propagation
    @inbounds for ilev in nlay:-1:1
        trans = exp(-τ[ilev, gcol] * Ds[1])
        flux_dn[ilev, gcol] = trans * flux_dn[ilev + 1, gcol] + src_dn[ilev, gcol]
    end

    # Surface reflection and emission
    @inbounds flux_up[1, gcol] =
        flux_dn[1, gcol] * (FT(1) - sfc_emis[ibnd, gcol]) + sfc_emis[ibnd, gcol] * sfc_source[gcol]

    # Upward propagation
    @inbounds for ilev in 2:(nlay + 1)
        trans = exp(-τ[ilev - 1, gcol] * Ds[1])
        flux_up[ilev, gcol] = trans * flux_up[ilev - 1, gcol] + src_up[ilev - 1, gcol]
    end

    @inbounds for ilev in 1:nlev
        flux_up[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_dn[ilev, gcol] *= (FT(2) * FT(π) * w_μ[1])
        flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
    end
end
