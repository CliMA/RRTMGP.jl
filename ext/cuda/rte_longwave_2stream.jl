function rte_lw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx, bx = _configure_threadblock(ncol)
    args = (flux_lw, src_lw, bcs_lw, op, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(
        args...,
    )
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
    igpt, ibnd = 1, 1
    if gcol ≤ ncol
        compute_optical_props!(op, as, src_lw, gcol)
        rte_lw_2stream!(
            op,
            flux_lw,
            src_lw,
            bcs_lw,
            gcol,
            igpt,
            ibnd,
            nlev,
            ncol,
        )
        compute_net_flux!(flux_lw, gcol)
    end
    return nothing
end

function rte_lw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    band_flux,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    set_band_flux_to_zero!(band_flux)
    tx, bx = _configure_threadblock(ncol)
    args = (
        flux,
        flux_lw,
        band_flux,
        src_lw,
        bcs_lw,
        op,
        nlay,
        ncol,
        as,
        lookup_lw,
        lookup_lw_cld,
        lookup_lw_aero,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_lw_2stream_solve_CUDA!(
    flux::FluxLW,
    flux_lw::FluxLW,
    band_flux,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    nlay,
    ncol,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        flux_up_lw = flux_lw.flux_up
        flux_dn_lw = flux_lw.flux_dn
        (; flux_up, flux_dn) = flux
        FT = eltype(flux_up)
        (; cloud_state, aerosol_state) = as
        if aerosol_state isa AerosolState
            Optics.compute_aero_mask!(
                view(aerosol_state.aero_mask, :, gcol),
                view(aerosol_state.aero_mass, :, :, gcol),
            )
        end
        n_cloudy_gpts = 0  # thread-local counter for LW cloud cover
        @inbounds for igpt in 1:n_gpt
            ibnd = major_gpt2bnd[igpt]
            if cloud_state isa CloudState
                Optics.build_cloud_mask!(
                    view(cloud_state.mask_lw, :, gcol),
                    view(cloud_state.cld_frac, :, gcol),
                    cloud_state.mask_type,
                )
                # count g-points with any cloudy layer
                n_cloudy_gpts += any(view(cloud_state.mask_lw, :, gcol)) ? 1 : 0
            end
            compute_optical_props!(
                op,
                as,
                src_lw,
                gcol,
                igpt,
                lookup_lw,
                lookup_lw_cld,
                lookup_lw_aero,
            )
            rte_lw_2stream!(
                op,
                flux,
                src_lw,
                bcs_lw,
                gcol,
                igpt,
                ibnd,
                nlev,
                ncol,
            )
            if igpt == 1
                for ilev in 1:nlev
                    flux_up_lw[ilev, gcol] = flux_up[ilev, gcol]
                    flux_dn_lw[ilev, gcol] = flux_dn[ilev, gcol]
                end
            else
                for ilev in 1:nlev
                    @inbounds flux_up_lw[ilev, gcol] += flux_up[ilev, gcol]
                    @inbounds flux_dn_lw[ilev, gcol] += flux_dn[ilev, gcol]
                end
            end
            # retain this g-point's contribution in its band (no-op when off)
            accumulate_band_flux!(band_flux, flux_up, flux_dn, gcol, ibnd, nlev)
        end
        @inbounds begin
            compute_net_flux!(flux_lw, gcol)
            # write out LW cloud cover
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_lw)
                cloud_state.cld_cover_lw[gcol] = FT(n_cloudy_gpts) / n_gpt
            end
        end
    end
    return nothing
end
