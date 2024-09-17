function rte_sw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx, bx = _configure_threadblock(ncol)
    args = (flux_sw, op, bcs_sw, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    return nothing
end

function rte_sw_noscat_solve_CUDA!(flux_sw::FluxSW, op::OneScalar, bcs_sw::SwBCs, nlay, ncol, as::GrayAtmosphericState)
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
                rte_sw_noscat!(flux_sw, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
            else
                set_flux!(flux_sw, FT(0), gcol)
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
    lookup_sw::LookUpSW,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx, bx = _configure_threadblock(ncol)
    args = (flux, flux_sw, op, bcs_sw, nlay, ncol, as, lookup_sw)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
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
    lookup_sw::LookUpSW,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    if gcol ≤ ncol
        μ₀ = bcs_sw.cos_zenith[gcol]
        @inbounds begin
            for igpt in 1:n_gpt
                compute_optical_props!(op, as, gcol, igpt, lookup_sw, nothing)
                solar_frac = lookup_sw.solar_src_scaled[igpt]
                rte_sw_noscat!(flux, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
                igpt == 1 ? set_flux!(flux_sw, flux, gcol) : add_to_flux!(flux_sw, flux, gcol)
            end
            if μ₀ <= 0
                set_flux!(flux_sw, FT(0), gcol)
                set_net_flux!(flux_sw, FT(0), gcol)
            else
                compute_net_flux!(flux_sw, gcol)
            end
        end
    end
    return nothing
end
