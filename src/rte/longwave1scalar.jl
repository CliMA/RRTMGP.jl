function rte_lw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_lw::FluxLW,
    src_lw::SourceLWNoScat,
    bcs_lw::LwBCs,
    op::OneScalar,
    angle_disc::AngularDiscretization,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    igpt, ibnd = 1, 1
    τ = op.τ
    Ds = angle_disc.gauss_Ds
    (; flux_up, flux_dn, flux_net) = flux_lw
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            compute_optical_props!(op, as, src_lw, gcol)
            rte_lw_noscat!(src_lw, bcs_lw, op, angle_disc, gcol, flux_lw, igpt, ibnd, nlay, nlev)
            for ilev in 1:nlev
                flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
            end
        end
    end
    return nothing
end

function rte_lw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    src_lw::SourceLWNoScat,
    bcs_lw::LwBCs,
    op::OneScalar,
    angle_disc::AngularDiscretization,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    τ = op.τ
    Ds = angle_disc.gauss_Ds
    flux_up_lw = flux_lw.flux_up
    flux_dn_lw = flux_lw.flux_dn
    flux_net_lw = flux_lw.flux_net
    @inbounds begin
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                ibnd = major_gpt2bnd[igpt]
                igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
                compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
                rte_lw_noscat!(src_lw, bcs_lw, op, angle_disc, gcol, flux, igpt, ibnd, nlay, nlev)
                add_to_flux!(flux_lw, flux, gcol)
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

"""
    lw_noscat_source_up(lev_source_inc::FT, lay_source::FT, τ_loc::FT, trans::FT, τ_thresh) where {FT}

Compute LW source function for upward emission at levels using linear-in-tau assumption
See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
"""
@inline function lw_noscat_source_up(lev_source_inc::FT, lay_source::FT, τ_loc::FT, trans::FT, τ_thresh) where {FT}
    fact =
        (τ_loc > τ_thresh) ? ((FT(1) - trans) / τ_loc - trans) :
        τ_loc * (FT(1 / 2) + τ_loc * (-FT(1 / 3) + τ_loc * FT(1 / 8)))
    # Equations below are developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
    return (FT(1) - trans) * lev_source_inc + FT(2) * fact * (lay_source - lev_source_inc)
end

"""
    lw_noscat_source_dn(lev_source_dec::FT, lay_source::FT, τ_loc::FT, trans::FT, τ_thresh) where {FT}

Compute LW source function for downward emission at levels using linear-in-tau assumption
See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
"""
@inline function lw_noscat_source_dn(lev_source_dec::FT, lay_source::FT, τ_loc::FT, trans::FT, τ_thresh) where {FT}
    fact =
        (τ_loc > τ_thresh) ? ((FT(1) - trans) / τ_loc - trans) :
        τ_loc * (FT(1 / 2) + τ_loc * (-FT(1 / 3) + τ_loc * FT(1 / 8)))
    # Equations below are developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
    return (FT(1) - trans) * lev_source_dec + FT(2) * fact * (lay_source - lev_source_dec)
end

"""
    rte_lw_noscat!(
        src_lw::SourceLWNoScat{FT},
        bcs_lw::LwBCs{FT},
        op::OneScalar{FT},
        gcol,
        flux::FluxLW,
        igpt,
        ibnd,
        nlay,
        nlev,
    ) where {FT<:AbstractFloat}

Transport for no-scattering longwave problem.
"""
@inline function rte_lw_noscat!(
    src_lw::SourceLWNoScat,
    bcs_lw::LwBCs,
    op::OneScalar,
    angle_disc::AngularDiscretization,
    gcol,
    flux::FluxLW,
    igpt,
    ibnd,
    nlay,
    nlev,
)
    # setting references
    (; sfc_source) = src_lw
    (; lay_source, lev_source_inc, lev_source_dec) = src_lw
    (; sfc_emis, inc_flux) = bcs_lw
    (; flux_up, flux_dn) = flux

    Ds = angle_disc.gauss_Ds
    w_μ = angle_disc.gauss_wts
    n_μ = angle_disc.n_gauss_angles
    τ = op.τ
    FT = eltype(τ)
    τ_thresh = 100 * eps(FT) # or abs(eps(FT))?

    intensity_to_flux = FT(π) * w_μ[1]
    flux_to_intensity = FT(1) / intensity_to_flux

    # Transport is for intensity
    #   convert flux at top of domain to intensity assuming azimuthal isotropy
    intensity_dn_ilevplus1 = isnothing(inc_flux) ? FT(0) : inc_flux[gcol, igpt] * flux_to_intensity
    @inbounds flux_dn[nlev, gcol] = intensity_dn_ilevplus1 * intensity_to_flux

    # Top of domain is index nlev
    # Downward propagation
    ilev = nlay
    @inbounds while ilev ≥ 1
        τ_loc = τ[ilev, gcol] * Ds[1]
        trans = exp(-τ_loc)
        lay_src, lev_src_dec = lay_source[ilev, gcol], lev_source_dec[ilev, gcol]
        intensity_dn_ilev =
            trans * intensity_dn_ilevplus1 + lw_noscat_source_dn(lev_src_dec, lay_src, τ_loc, trans, τ_thresh)
        intensity_dn_ilevplus1 = intensity_dn_ilev
        flux_dn[ilev, gcol] = intensity_dn_ilev * intensity_to_flux
        ilev -= 1
    end

    # Surface reflection and emission
    @inbounds intensity_up_ilevminus1 =
        intensity_dn_ilevplus1 * (FT(1) - sfc_emis[ibnd, gcol]) + sfc_emis[ibnd, gcol] * sfc_source[gcol]
    #flux_dn[1, gcol] * (FT(1) - sfc_emis[ibnd, gcol]) + sfc_emis[ibnd, gcol] * sfc_source[gcol]
    @inbounds flux_up[1, gcol] = intensity_up_ilevminus1 * intensity_to_flux

    # Upward propagation
    @inbounds for ilev in 2:(nlay + 1)
        τ_loc = τ[ilev - 1, gcol] * Ds[1]
        trans = exp(-τ_loc)
        lay_src, lev_src_inc = lay_source[ilev - 1, gcol], lev_source_inc[ilev - 1, gcol]
        intensity_up_ilev =
            trans * intensity_up_ilevminus1 + lw_noscat_source_up(lev_src_inc, lay_src, τ_loc, trans, τ_thresh)
        intensity_up_ilevminus1 = intensity_up_ilev
        flux_up[ilev, gcol] = intensity_up_ilev * intensity_to_flux
    end
    return nothing
end
