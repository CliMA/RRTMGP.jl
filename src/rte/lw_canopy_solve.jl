using ..Canopy: CanopyState, fill_canopy_lw_optics!, fill_canopy_lw_sources!, planck_interp

"""
    rte_lw_2stream_canopy_solve!(device, flux, flux_lw, op_atm, op_exp, bcs_lw, src_atm, src_exp, as, lookup_lw, lookup_lw_cld, lookup_lw_aero, canopy)

CPU dispatch for the coupled atmosphere-canopy longwave 2-stream solver.

For each g-point and column:
1. Computes atmospheric optical properties and Planck sources
2. Fills expanded column optics and temperature-scaled Planck sources
3. Solves using the standard `rte_lw_2stream!` (no coefficient changes for LW)
"""
function rte_lw_2stream_canopy_solve!(
    device::ClimaComms.AbstractCPUDevice,
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
    nlay_atm, ncol = AtmosphericStates.get_dims(as)
    ncanlay = canopy.ncanlay
    nlay_total = nlay_atm + ncanlay
    nlev_total = nlay_total + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    (; cloud_state, aerosol_state) = as
    bld_cld_mask = cloud_state isa CloudState
    flux_up_lw = flux_lw.flux_up
    flux_dn_lw = flux_lw.flux_dn
    flux_net_lw = flux_lw.flux_net
    (; flux_up, flux_dn) = flux

    # Extract Planck lookup tables as CPU arrays for canopy T-profile scaling
    t_planck_cpu = Array(lookup_lw.planck.t_planck)
    tot_planck_cpu = Array(lookup_lw.planck.tot_planck)

    @inbounds begin
        if aerosol_state isa AerosolState
            ClimaComms.@threaded device for gcol in 1:ncol
                Optics.compute_aero_mask!(
                    view(aerosol_state.aero_mask, :, gcol),
                    view(aerosol_state.aero_mass, :, :, gcol),
                )
            end
        end
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                ibnd = major_gpt2bnd[igpt]
                bld_cld_mask && Optics.build_cloud_mask!(
                    view(cloud_state.mask_lw, :, gcol),
                    view(cloud_state.cld_frac, :, gcol),
                    cloud_state.mask_type,
                )
                # Compute atmospheric optical properties + Planck sources
                compute_optical_props!(op_atm, as, src_atm, gcol, igpt, lookup_lw, lookup_lw_cld, lookup_lw_aero)

                # Fill expanded column optics
                fill_canopy_lw_optics!(op_exp, op_atm, canopy, ibnd, gcol, nlay_atm)

                # Fill expanded column Planck sources with canopy temperature profile
                totplnk_band = view(tot_planck_cpu, :, ibnd)
                fill_canopy_lw_sources!(src_exp, src_atm, canopy, ibnd, gcol, nlay_atm, t_planck_cpu, totplnk_band)

                # Solve using standard LW 2-stream (no coefficient changes needed)
                rte_lw_2stream!(op_exp, flux, src_exp, bcs_lw, gcol, igpt, ibnd, nlev_total, ncol)

                # Accumulate fluxes
                if igpt == 1
                    map!(x -> x, view(flux_up_lw, :, gcol), view(flux_up, :, gcol))
                    map!(x -> x, view(flux_dn_lw, :, gcol), view(flux_dn, :, gcol))
                else
                    for ilev in 1:nlev_total
                        @inbounds flux_up_lw[ilev, gcol] += flux_up[ilev, gcol]
                        @inbounds flux_dn_lw[ilev, gcol] += flux_dn[ilev, gcol]
                    end
                end
            end
        end
        ClimaComms.@threaded device for gcol in 1:ncol
            for ilev in 1:nlev_total
                flux_net_lw[ilev, gcol] = flux_up_lw[ilev, gcol] - flux_dn_lw[ilev, gcol]
            end
        end
    end
    return nothing
end
