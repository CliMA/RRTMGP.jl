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
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(args...)
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
        (; flux_up, flux_dn, flux_net) = flux_lw
        compute_optical_props!(op, as, src_lw, gcol)
        rte_lw_2stream!(op, flux_lw, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
        @inbounds begin
            for ilev in 1:nlev
                flux_net[gcol, ilev] = flux_up[gcol, ilev] - flux_dn[gcol, ilev]
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
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx, bx = _configure_threadblock(ncol)
    args = (flux, flux_lw, src_lw, bcs_lw, op, nlay, ncol, as, lookup_lw, lookup_lw_cld, lookup_lw_aero)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(args...)
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
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        flux_up_lw = flux_lw.flux_up
        flux_dn_lw = flux_lw.flux_dn
        flux_net_lw = flux_lw.flux_net
        (; flux_up, flux_dn) = flux
        (; cloud_state, aerosol_state) = as
        if aerosol_state isa AerosolState
            Optics.compute_aero_mask!(view(aerosol_state.aero_mask, :, gcol), view(aerosol_state.aero_mass, :, :, gcol))
        end
        @inbounds for igpt in 1:n_gpt
            ibnd = major_gpt2bnd[igpt]
            if cloud_state isa CloudState
                Optics.build_cloud_mask!(
                    view(cloud_state.mask_lw, :, gcol),
                    view(cloud_state.cld_frac, :, gcol),
                    cloud_state.mask_type,
                )
            end
            compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld, lookup_lw_aero)
            rte_lw_2stream!(op, flux, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
            if igpt == 1
                map!(x -> x, view(flux_up_lw, gcol, :), view(flux_up, gcol, :))
                map!(x -> x, view(flux_dn_lw, gcol, :), view(flux_dn, gcol, :))
            else
                for ilev in 1:nlev
                    @inbounds flux_up_lw[gcol, ilev] += flux_up[gcol, ilev]
                    @inbounds flux_dn_lw[gcol, ilev] += flux_dn[gcol, ilev]
                end
            end
        end
        @inbounds begin
            for ilev in 1:nlev
                flux_net_lw[gcol, ilev] = flux_up_lw[gcol, ilev] - flux_dn_lw[gcol, ilev]
            end
        end
    end
    return nothing
end
