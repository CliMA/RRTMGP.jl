function rte_sw_2stream_canopy_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op_atm::TwoStream,
    op_exp::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
    canopy::CanopyState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    args = (
        flux, flux_sw, op_atm, op_exp, bcs_sw, src_sw,
        nlay, ncol, as, lookup_sw, lookup_sw_cld, lookup_sw_aero, canopy,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_canopy_solve_CUDA!(args...)
    return nothing
end

function rte_sw_2stream_canopy_solve_CUDA!(
    flux::FluxSW,
    flux_sw::FluxSW,
    op_atm::TwoStream,
    op_exp::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay_atm,
    ncol,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
    canopy::CanopyState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    ncanlay = canopy.ncanlay
    nlev_total = nlay_atm + ncanlay + 1
    n_gpt = length(lookup_sw.band_data.major_gpt2bnd)
    if gcol ≤ ncol
        flux_up_sw = flux_sw.flux_up
        flux_dn_sw = flux_sw.flux_dn
        flux_dn_dir_sw = flux_sw.flux_dn_dir
        flux_net_sw = flux_sw.flux_net
        flux_up = flux.flux_up
        flux_dn = flux.flux_dn
        flux_dn_dir = flux.flux_dn_dir
        FT = eltype(flux_up)
        (; cloud_state, aerosol_state) = as
        μ₀ = bcs_sw.cos_zenith[gcol]
        @inbounds begin
            if aerosol_state isa AerosolState
                Optics.compute_aero_mask!(
                    view(aerosol_state.aero_mask, :, gcol),
                    view(aerosol_state.aero_mass, :, :, gcol),
                )
            end
            for igpt in 1:n_gpt
                if cloud_state isa CloudState
                    Optics.build_cloud_mask!(
                        view(cloud_state.mask_sw, :, gcol),
                        view(cloud_state.cld_frac, :, gcol),
                        cloud_state.mask_type,
                    )
                end
                compute_optical_props!(op_atm, as, gcol, igpt, lookup_sw, lookup_sw_cld, lookup_sw_aero)
                solar_frac = lookup_sw.solar_src_scaled[igpt]
                ibnd = lookup_sw.band_data.major_gpt2bnd[igpt]
                fill_canopy_sw_optics!(op_exp, op_atm, canopy, ibnd, gcol, nlay_atm)
                rte_sw_2stream!(op_exp, src_sw, bcs_sw, flux, solar_frac, igpt, n_gpt, ibnd, nlev_total, gcol, canopy)
                if igpt == 1
                    map!(x -> x, view(flux_up_sw, :, gcol), view(flux_up, :, gcol))
                    map!(x -> x, view(flux_dn_sw, :, gcol), view(flux_dn, :, gcol))
                    map!(x -> x, view(flux_dn_dir_sw, :, gcol), view(flux_dn_dir, :, gcol))
                else
                    for ilev in 1:nlev_total
                        flux_up_sw[ilev, gcol] += flux_up[ilev, gcol]
                        flux_dn_sw[ilev, gcol] += flux_dn[ilev, gcol]
                    end
                    flux_dn_dir_sw[1, gcol] += flux_dn_dir[1, gcol]
                end
            end
            if μ₀ ≤ 0
                for ilev in 1:nlev_total
                    flux_up_sw[ilev, gcol] = FT(0)
                end
                for ilev in 1:nlev_total
                    flux_dn_sw[ilev, gcol] = FT(0)
                end
            end
            for ilev in 1:nlev_total
                flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
            end
        end
    end
    return nothing
end
