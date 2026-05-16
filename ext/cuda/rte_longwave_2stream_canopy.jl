function rte_lw_2stream_canopy_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    op_atm::TwoStream,
    op_exp::TwoStream,
    bcs_lw::LwBCs,
    src_atm::SourceLW2Str,
    src_exp::SourceLW2Str,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
    canopy::CanopyState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    args = (
        flux, flux_lw, op_atm, op_exp, bcs_lw, src_atm, src_exp,
        nlay, ncol, as, lookup_lw, lookup_lw_cld, lookup_lw_aero, canopy,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_canopy_solve_CUDA!(args...)
    return nothing
end

function rte_lw_2stream_canopy_solve_CUDA!(
    flux::FluxLW,
    flux_lw::FluxLW,
    op_atm::TwoStream,
    op_exp::TwoStream,
    bcs_lw::LwBCs,
    src_atm::SourceLW2Str,
    src_exp::SourceLW2Str,
    nlay_atm,
    ncol,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
    canopy::CanopyState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    ncanlay = canopy.ncanlay
    nlev_total = nlay_atm + ncanlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        flux_up_lw = flux_lw.flux_up
        flux_dn_lw = flux_lw.flux_dn
        flux_net_lw = flux_lw.flux_net
        flux_up = flux.flux_up
        flux_dn = flux.flux_dn
        FT = eltype(flux_up)
        (; cloud_state, aerosol_state) = as
        # Planck tables are on-device (via Adapt)
        t_planck = lookup_lw.planck.t_planck
        tot_planck = lookup_lw.planck.tot_planck
        @inbounds begin
            if aerosol_state isa AerosolState
                Optics.compute_aero_mask!(
                    view(aerosol_state.aero_mask, :, gcol),
                    view(aerosol_state.aero_mass, :, :, gcol),
                )
            end
            for igpt in 1:n_gpt
                ibnd = major_gpt2bnd[igpt]
                if cloud_state isa CloudState
                    Optics.build_cloud_mask!(
                        view(cloud_state.mask_lw, :, gcol),
                        view(cloud_state.cld_frac, :, gcol),
                        cloud_state.mask_type,
                    )
                end
                compute_optical_props!(op_atm, as, src_atm, gcol, igpt, lookup_lw, lookup_lw_cld, lookup_lw_aero)
                fill_canopy_lw_optics!(op_exp, op_atm, canopy, ibnd, gcol, nlay_atm)
                totplnk_band = view(tot_planck, :, ibnd)
                fill_canopy_lw_sources!(src_exp, src_atm, canopy, ibnd, gcol, nlay_atm, t_planck, totplnk_band)
                rte_lw_2stream!(op_exp, flux, src_exp, bcs_lw, gcol, igpt, ibnd, nlev_total, ncol)
                if igpt == 1
                    map!(x -> x, view(flux_up_lw, :, gcol), view(flux_up, :, gcol))
                    map!(x -> x, view(flux_dn_lw, :, gcol), view(flux_dn, :, gcol))
                else
                    for ilev in 1:nlev_total
                        flux_up_lw[ilev, gcol] += flux_up[ilev, gcol]
                        flux_dn_lw[ilev, gcol] += flux_dn[ilev, gcol]
                    end
                end
            end
            for ilev in 1:nlev_total
                flux_net_lw[ilev, gcol] = flux_up_lw[ilev, gcol] - flux_dn_lw[ilev, gcol]
            end
        end
    end
    return nothing
end
