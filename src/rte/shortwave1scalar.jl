function rte_sw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    max_threads,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt, igpt = 1, 1
    FT = eltype(op.angle_disc.gauss_Ds)
    solar_frac = FT(1)
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol)
            rte_sw_noscat!(flux_sw, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
        end
    end
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    max_threads,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux_sw, op, bcs_sw, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    return nothing
end

function rte_sw_noscat_solve_CUDA!(flux_sw::FluxSW, op::OneScalar, bcs_sw::SwBCs, nlay, ncol, as::GrayAtmosphericState)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt = 1, 1
    FT = eltype(op.angle_disc.gauss_Ds)
    solar_frac = FT(1)
    # setting references for flux_sw
    if gcol ≤ ncol
        @inbounds begin
            set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol)
            rte_sw_noscat!(flux_sw, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
        end
    end
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    max_threads,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    # setting references for flux_sw
    @inbounds begin
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
                solar_frac = lookup_sw.solar_src_scaled[igpt]
                rte_sw_noscat!(flux, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
                n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
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
    max_threads,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    # setting references for flux_sw
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux, flux_sw, op, bcs_sw, nlay, ncol, as, lookup_sw, lookup_sw_cld)
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
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    # setting references for flux_sw
    if gcol ≤ ncol
        set_flux_to_zero!(flux_sw, gcol)
        @inbounds for igpt in 1:n_gpt
            compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
            solar_frac = lookup_sw.solar_src_scaled[igpt]
            rte_sw_noscat!(flux, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
            n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
        end
    end
    return nothing
end

"""
    rte_sw_noscat!(
        flux::FluxSW,
        op::OneScalar,
        bcs_sw::SwBCs,
        solar_frac::AbstractFloat,
        gcol,
        nlev,
    )

No-scattering solver for the shortwave problem.
(Extinction-only i.e. solar direct beam)
"""
function rte_sw_noscat!(
    flux::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    igpt::Int,
    n_gpt::Int,
    solar_frac::AbstractFloat,
    gcol::Int,
    nlev::Int,
)
    (; toa_flux, cos_zenith) = bcs_sw
    τ = op.τ
    (; flux_dn_dir, flux_net) = flux
    # downward propagation
    @inbounds flux_dn_dir[nlev, gcol] = toa_flux[gcol] * solar_frac * cos_zenith[gcol]
    ilev = nlev - 1
    @inbounds while ilev ≥ 1
        flux_dn_dir[ilev, gcol] = flux_dn_dir[ilev + 1, gcol] * exp(-τ[ilev, gcol] / cos_zenith[gcol])
        flux_net[ilev, gcol] = -flux_dn_dir[ilev, gcol]
        ilev -= 1
    end
end
