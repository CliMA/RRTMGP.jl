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
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(
        args...,
    )
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
            rte_sw_2stream!(
                op,
                src_sw,
                bcs_sw,
                flux_sw,
                solar_frac,
                igpt,
                n_gpt,
                ibnd,
                nlev,
                gcol,
            )
            for ilev in 1:nlev
                flux_net_sw[ilev, gcol] =
                    flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
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
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    set_band_flux_to_zero!(band_flux)
    tx, bx = _configure_threadblock(ncol)
    args = (
        flux,
        flux_sw,
        band_flux,
        op,
        bcs_sw,
        src_sw,
        nlay,
        ncol,
        as,
        lookup_sw,
        lookup_sw_cld,
        lookup_sw_aero,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux::FluxSW,
    flux_sw::FluxSW,
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay,
    ncol,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.band_data.major_gpt2bnd)
    if gcol ≤ ncol
        FT = eltype(flux_sw.flux_up)
        (; cloud_state, aerosol_state) = as
        μ₀ = @inbounds bcs_sw.cos_zenith[gcol]
        n_cloudy_gpts = 0  # thread-local counter for cloud cover
        @inbounds begin
            _compute_aero_mask!(aerosol_state, gcol)
            for igpt in 1:n_gpt
                cloudy = sw_2stream_gpt_col!(
                    igpt,
                    gcol,
                    flux,
                    flux_sw,
                    band_flux,
                    op,
                    bcs_sw,
                    src_sw,
                    as,
                    lookup_sw,
                    lookup_sw_cld,
                    lookup_sw_aero,
                    μ₀,
                    lookup_sw.band_data.major_gpt2bnd[igpt],
                    n_gpt,
                    nlev,
                )
                n_cloudy_gpts += cloudy ? 1 : 0
            end
            if μ₀ ≤ 0 # zero out columns with zenith angle ≥ π/2
                set_flux_to_zero!(flux_sw, gcol)
            else
                compute_net_flux!(flux_sw, gcol)

            end
            # write out SW cloud cover
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_sw)
                cloud_state.cld_cover_sw[gcol] = FT(n_cloudy_gpts) / n_gpt
            end
        end
    end
    return nothing
end
