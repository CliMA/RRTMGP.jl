function rte_sw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx, bx = _configure_threadblock(ncol)
    args = (flux_sw, op, bcs_sw, src_sw, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(args...)
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay,
    ncol,
    as::GrayAtmosphericState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt, ibnd = 1, 1, 1
    FT = eltype(bcs_sw.cos_zenith)
    solar_frac = FT(1)
    if gcol ≤ ncol
        flux_up_sw = flux_sw.flux_up
        flux_dn_sw = flux_sw.flux_dn
        flux_net_sw = flux_sw.flux_net
        μ₀ = bcs_sw.cos_zenith[gcol]
        @inbounds begin
            compute_optical_props!(op, as, gcol)
            # call shortwave rte solver
            rte_sw_2stream!(op, src_sw, bcs_sw, flux_sw, solar_frac, igpt, n_gpt, ibnd, nlev, gcol)
            for ilev in 1:nlev
                flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
            end
        end
        if μ₀ ≤ 0 # zero out columns with zenith angle ≥ π/2
            for ilev in 1:nlev
                flux_up_sw[ilev, gcol] = FT(0)
            end
            for ilev in 1:nlev
                flux_dn_sw[ilev, gcol] = FT(0)
            end
            for ilev in 1:nlev
                flux_net_sw[ilev, gcol] = FT(0)
            end
        end
    end
    return nothing
end

function rte_sw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    tx, bx = _configure_threadblock(ncol)
    args = (flux, flux_sw, op, bcs_sw, src_sw, nlay, ncol, as, lookup_sw, lookup_sw_cld, lookup_sw_aero)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(args...)
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux::FluxSW,
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay,
    ncol,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
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
                # set cloud mask, if applicable
                if cloud_state isa CloudState
                    Optics.build_cloud_mask!(
                        view(cloud_state.mask_sw, :, gcol),
                        view(cloud_state.cld_frac, :, gcol),
                        cloud_state.mask_type,
                    )
                end
                # compute optical properties
                compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld, lookup_sw_aero)
                solar_frac = lookup_sw.solar_src_scaled[igpt]
                ibnd = lookup_sw.band_data.major_gpt2bnd[igpt]
                # rte shortwave solver
                rte_sw_2stream!(op, src_sw, bcs_sw, flux, solar_frac, igpt, n_gpt, ibnd, nlev, gcol)
                if igpt == 1
                    map!(x -> x, view(flux_up_sw, :, gcol), view(flux_up, :, gcol))
                    map!(x -> x, view(flux_dn_sw, :, gcol), view(flux_dn, :, gcol))
                    map!(x -> x, view(flux_dn_dir_sw, :, gcol), view(flux_dn_dir, :, gcol))
                else
                    for ilev in 1:nlev
                        flux_up_sw[ilev, gcol] += flux_up[ilev, gcol]
                        flux_dn_sw[ilev, gcol] += flux_dn[ilev, gcol]
                        flux_dn_dir_sw[ilev, gcol] += flux_dn_dir[ilev, gcol]
                    end
                end
            end
            if μ₀ ≤ 0 # zero out columns with zenith angle ≥ π/2
                for ilev in 1:nlev
                    flux_up_sw[ilev, gcol] = FT(0)
                end
                for ilev in 1:nlev
                    flux_dn_sw[ilev, gcol] = FT(0)
                end
            end
            for ilev in 1:nlev
                flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
            end
        end
    end
    return nothing
end
