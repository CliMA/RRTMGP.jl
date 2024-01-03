function rte_sw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_sw::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    max_threads,
    as::GrayAtmosphericState{FT},
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt, igpt = 1, 1
    solar_frac = FT(1)
    # setting references for flux_sw
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol, igpt, nothing, nothing)
            rte_sw_noscat_solve_kernel!(flux_sw, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
        end
    end
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    max_threads,
    as::GrayAtmosphericState{FT},
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    # setting references for flux_sw
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux_sw, op, bcs_sw, nlay, ncol, as)
    @cuda threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    return nothing
end

function rte_sw_noscat_solve_CUDA!(
    flux_sw::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    nlay,
    ncol,
    as::GrayAtmosphericState{FT},
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt = 1, 1
    solar_frac = FT(1)
    # setting references for flux_sw
    if gcol ≤ ncol
        @inbounds begin
            set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol, igpt, nothing, nothing)
            rte_sw_noscat_solve_kernel!(flux_sw, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
        end
    end
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    max_threads,
    as::AtmosphericState{FT},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    τ = op.τ
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    # setting references for flux_sw
    @inbounds begin
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                for glay in 1:nlay
                    τ[glay, gcol], _, _ = compute_optical_props(as, gcol, glay, igpt, lookup_sw, lookup_sw_cld)
                end
                solar_frac = as isa GrayAtmosphericState ? FT(1) : lookup_sw.solar_src_scaled[igpt]
                rte_sw_noscat_solve_kernel!(flux, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
                n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
            end
        end
    end
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    max_threads,
    as::AtmosphericState{FT},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    # setting references for flux_sw
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux, flux_sw, op, bcs_sw, nlay, ncol, as, lookup_sw, lookup_sw_cld)
    @cuda threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    return nothing
end

function rte_sw_noscat_solve_CUDA!(
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    nlay,
    ncol,
    as::AtmosphericState{FT},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    # setting references for flux_sw
    if gcol ≤ ncol
        τ = op.τ
        set_flux_to_zero!(flux_sw, gcol)
        @inbounds for igpt in 1:n_gpt
            for glay in 1:nlay
                τ[glay, gcol] = compute_optical_props(as, gcol, glay, igpt, lookup_sw, lookup_sw_cld)
            end
            solar_frac = as isa GrayAtmosphericState ? FT(1) : lookup_sw.solar_src_scaled[igpt]
            rte_sw_noscat_solve_kernel!(flux, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
            n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
        end
    end
    return nothing
end

"""
    rte_sw_noscat_solve_kernel!(
        flux::FluxSW{FT},
        op::OneScalar{FT},
        bcs_sw::SwBCs{FT},
        solar_frac::FT,
        gcol,
        nlev,
    ) where {FT<:AbstractFloat}

No-scattering solver for the shortwave problem.
(Extinction-only i.e. solar direct beam)
"""
function rte_sw_noscat_solve_kernel!(
    flux::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    igpt::Int,
    n_gpt::Int,
    solar_frac::FT,
    gcol::Int,
    nlev::Int,
) where {FT <: AbstractFloat}
    #solar_frac = solar_src_scaled[igpt]
    (; toa_flux, cos_zenith) = bcs_sw
    #n_gpt = length(solar_src_scaled)
    τ = op.τ
    (; flux_dn_dir, flux_net) = flux
    # downward propagation
    @inbounds flux_dn_dir[nlev, gcol] = toa_flux[gcol] * solar_frac * cos_zenith[gcol]
    @inbounds for ilev in (nlev - 1):-1:1
        flux_dn_dir[ilev, gcol] = flux_dn_dir[ilev + 1, gcol] * exp(-τ[ilev, gcol] / cos_zenith[gcol])
        flux_net[ilev, gcol] = -flux_dn_dir[ilev, gcol]
    end
end
