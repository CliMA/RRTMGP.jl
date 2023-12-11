"""
    rte_sw_noscat_solve!(
        context,
        flux::FluxSW{FT},
        flux_sw::FluxSW{FT},
        op::OneScalar{FT},
        bcs_sw::SwBCs{FT},
        solar_src_scaled::AbstractArray{FT,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lookup_sw::Union{LookUpSW, Nothing} = nothing,
        lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    ) where {FT<:AbstractFloat}

No-scattering solver for the shortwave problem.
(Extinction-only i.e. solar direct beam)
"""
function rte_sw_noscat_solve!(
    context,
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    solar_src_scaled::AbstractArray{FT, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt = length(solar_src_scaled)
    device = ClimaComms.device(context)
    # setting references for flux_sw
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, flux_sw, op, bcs_sw, nlay, ncol, solar_src_scaled, as, lookup_sw, lookup_sw_cld)
        @cuda threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    else
        @inbounds begin
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                    compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
                    rte_sw_noscat_solve_kernel!(flux, op, bcs_sw, igpt, solar_src_scaled, gcol, nlev)
                    n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
                end
            end
        end
    end
    return nothing
end

function rte_sw_noscat_solve_CUDA!(
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    nlay,
    ncol,
    solar_src_scaled::AbstractArray{FT, 1},
    as::AbstractAtmosphericState{FT},
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(solar_src_scaled)
    # setting references for flux_sw
    if gcol ≤ ncol
        @inbounds for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
            rte_sw_noscat_solve_kernel!(flux, op, bcs_sw, igpt, solar_src_scaled, gcol, nlev)
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
    solar_src_scaled::AbstractArray{FT, 1},
    gcol::Int,
    nlev::Int,
) where {FT <: AbstractFloat}
    solar_frac = solar_src_scaled[igpt]
    (; toa_flux, zenith) = bcs_sw
    n_gpt = length(solar_src_scaled)
    τ = op.τ
    (; flux_dn_dir, flux_net) = flux
    # downward propagation
    @inbounds flux_dn_dir[nlev, gcol] = toa_flux[gcol] * solar_frac * cos(zenith[gcol])
    @inbounds for ilev in (nlev - 1):-1:1
        flux_dn_dir[ilev, gcol] = flux_dn_dir[ilev + 1, gcol] * exp(-τ[ilev, gcol] / cos(zenith[gcol]))
        flux_net[ilev, gcol] = -flux_dn_dir[ilev, gcol]
    end
end
