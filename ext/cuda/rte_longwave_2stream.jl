function rte_lw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux_lw::FluxLW,
    src_lw::SourceLW2Str,
    bcs_lw::LwBCs,
    op::TwoStream,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
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
        compute_net_flux!(flux_lw, gcol, nlev)
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
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
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
        state_cache,
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
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, Nothing},
    lookup_lw_aero::Union{LookUpAerosolMerra, Nothing},
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        (; cloud_state, aerosol_state) = as
        FT = eltype(flux_lw.flux_up)
        _compute_aero_mask!(aerosol_state, gcol)
        n_cloudy_gpts = 0  # thread-local counter for LW cloud cover
        @inbounds for igpt in 1:n_gpt
            cloudy = lw_2stream_gpt_col!(
                igpt,
                gcol,
                flux,
                flux_lw,
                band_flux,
                src_lw,
                bcs_lw,
                op,
                as,
                state_cache,
                lookup_lw,
                lookup_lw_cld,
                lookup_lw_aero,
                major_gpt2bnd[igpt],
                nlev,
                ncol,
            )
            n_cloudy_gpts += cloudy ? 1 : 0
        end
        @inbounds begin
            compute_net_flux!(flux_lw, gcol, nlev)
            # write out LW cloud cover
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_lw)
                cloud_state.cld_cover_lw[gcol] = FT(n_cloudy_gpts) / n_gpt
            end
        end
    end
    return nothing
end
