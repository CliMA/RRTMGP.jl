function rte_sw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    args = (flux_sw, op, bcs_sw, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_sw_noscat_solve_CUDA!(
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    nlay,
    ncol,
    as::GrayAtmosphericState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt = 1, 1
    FT = eltype(bcs_sw.cos_zenith)
    solar_frac = FT(1)
    if gcol ≤ ncol
        @inbounds begin
            μ₀ = bcs_sw.cos_zenith[gcol]
            if μ₀ > 0
                compute_optical_props!(op, as, gcol)
                rte_sw_noscat!(
                    flux_sw,
                    op,
                    bcs_sw,
                    igpt,
                    n_gpt,
                    solar_frac,
                    gcol,
                    nlev,
                )
                compute_net_flux!(flux_sw, gcol, nlev)
            else
                set_flux_to_zero!(flux_sw, gcol, nlev)
            end
        end
    end
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    args = (flux, flux_sw, op, bcs_sw, nlay, ncol, as, state_cache, lookup_sw)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_sw_noscat_solve_CUDA!(
    flux::FluxSW,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    nlay,
    ncol,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    if gcol ≤ ncol
        μ₀ = @inbounds bcs_sw.cos_zenith[gcol]
        @inbounds begin
            if μ₀ > 0
                for igpt in 1:n_gpt
                    sw_noscat_gpt_col!(
                        igpt,
                        gcol,
                        flux,
                        flux_sw,
                        op,
                        bcs_sw,
                        as,
                        state_cache,
                        lookup_sw,
                        n_gpt,
                        nlev,
                    )
                end
                compute_net_flux!(flux_sw, gcol, nlev)
            else
                set_flux_to_zero!(flux_sw, gcol, nlev)
            end
        end
    end
    return nothing
end
