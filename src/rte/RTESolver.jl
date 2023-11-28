module RTESolver

using StaticArrays
import ClimaComms
using CUDA
using Adapt

using ..AngularDiscretizations
using ..Vmrs
using ..AtmosphericStates
using ..Fluxes
using ..BCs
using ..Sources
using ..RTE
using ..Optics
using ..LookUpTables

export solve_lw!, solve_sw!

include("RTESolverKernels.jl")

"""
    solve_lw!(
        slv::Solver,
        max_threads::Int,
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    )

Solver for the longwave radiation problem
"""
function solve_lw!(
    slv::Solver,
    max_threads::Int,
    lookup_lw::Union{LookUpLW, Nothing} = nothing,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
)
    (; as, op, bcs_lw, src_lw, flux_lw, fluxb_lw) = slv
    context = RTE.context(slv)
    DA = RTE.array_type(slv)
    no_args = lookup_lw isa Nothing && lookup_lw_cld isa Nothing
    if no_args
        major_gpt2bnd = DA(UInt8[1])
        flux = flux_lw
    else
        (; major_gpt2bnd) = lookup_lw
        flux = fluxb_lw
    end

    rte_lw_solve!(context, flux, flux_lw, src_lw, bcs_lw, op, major_gpt2bnd, max_threads, as, lookup_lw, lookup_lw_cld)
    return nothing
end

"""
    rte_lw_solve!(
        context,
        flux::FluxLW{FT},
        flux_lw::FluxLW{FT},
        src_lw::SourceLWNoScat{FT},
        bcs_lw::LwBCs{FT},
        op::OneScalar{FT},
        major_gpt2bnd::AbstractArray{UInt8,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    ) where {FT<:AbstractFloat}

No scattering solver for the longwave problem.
LW fluxes, no scattering, mu (cosine of integration angle) specified by column
Does radiation calculation at user-supplied angles; converts radiances to flux
using user-supplied weights. Currently, only n_gauss_angles = 1 supported.
"""
function rte_lw_solve!(
    context,
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
    major_gpt2bnd::AbstractArray{UInt8, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lookup_lw::Union{LookUpLW, Nothing} = nothing,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, flux_lw, src_lw, bcs_lw, op, nlay, ncol, major_gpt2bnd, as, lookup_lw, lookup_lw_cld)
        @cuda threads = (tx) blocks = (bx) rte_lw_noscat_solve_CUDA!(args...)
    else
        @inbounds begin
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
                    compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
                    rte_lw_noscat_source!(src_lw, op, gcol, nlay)
                    rte_lw_noscat_transport!(src_lw, bcs_lw, op, gcol, flux, igpt, major_gpt2bnd, nlay, nlev)
                    n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
                end
            end
        end
    end
    #------------------------------------------------
    return nothing
end
# ---------------------------------------------------------------    
function rte_lw_noscat_solve_CUDA!(
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
    nlay,
    ncol,
    major_gpt2bnd,
    as::AbstractAtmosphericState{FT},
    lookup_lw::Union{LookUpLW, Nothing} = nothing,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        @inbounds for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
            compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
            rte_lw_noscat_source!(src_lw, op, gcol, nlay)
            rte_lw_noscat_transport!(src_lw, bcs_lw, op, gcol, flux, igpt, major_gpt2bnd, nlay, nlev)
            n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
        end
    end
    return nothing
end

"""
    rte_lw_solve!(
        context,
        flux::FluxLW{FT},
        flux_lw::FluxLW{FT},
        src_lw::SourceLW2Str{FT},
        bcs_lw::LwBCs{FT},
        op::TwoStream{FT},
        major_gpt2bnd::AbstractArray{UInt8,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    ) where {FT<:AbstractFloat}

Two stream solver for the longwave problem.

RRTMGP provides source functions at each level using the spectral mapping
of each adjacent layer. Combine these for two-stream calculations
lw combine sources (rte_lw_2stream_combine_sources!)
combine RRTMGP-specific sources at levels
compute layer reflectance, transmittance
compute total source function at levels using linear-in-tau (rte_lw_2stream_source!)
transport (adding_lw!)
"""
function rte_lw_solve!(
    context,
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLW2Str{FT},
    bcs_lw::LwBCs{FT},
    op::TwoStream{FT},
    major_gpt2bnd::AbstractArray{UInt8, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lookup_lw::Union{LookUpLW, Nothing} = nothing,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, flux_lw, src_lw, bcs_lw, op, nlay, ncol, major_gpt2bnd, as, lookup_lw, lookup_lw_cld)
        @cuda threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(args...)
    else
        bld_cld_mask = as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
        @inbounds begin
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
                    bld_cld_mask && Optics.build_cloud_mask!(
                        as.cld_mask_lw,
                        as.cld_frac,
                        as.random_lw,
                        gcol,
                        igpt,
                        as.cld_mask_type,
                    )
                    compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
                    rte_lw_2stream_combine_sources!(src_lw, gcol, nlev, ncol)
                    rte_lw_2stream_source!(op, src_lw, gcol, nlay, ncol)
                    adding_lw!(flux, src_lw, bcs_lw, gcol, igpt, major_gpt2bnd, nlev, ncol)
                    n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
                end
            end
        end
    end
    return nothing
end
# -------------------------------------------------------------------------------------------------
function rte_lw_2stream_solve_CUDA!(
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLW2Str{FT},
    bcs_lw::LwBCs{FT},
    op::TwoStream{FT},
    nlay,
    ncol,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    as::AbstractAtmosphericState{FT},
    lookup_lw::Union{LookUpLW, Nothing} = nothing,
    lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        @inbounds for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
            if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
                Optics.build_cloud_mask!(as.cld_mask_lw, as.cld_frac, as.random_lw, gcol, igpt, as.cld_mask_type)
            end
            compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
            rte_lw_2stream_combine_sources!(src_lw, gcol, nlev, ncol)
            rte_lw_2stream_source!(op, src_lw, gcol, nlay, ncol)
            adding_lw!(flux, src_lw, bcs_lw, gcol, igpt, major_gpt2bnd, nlev, ncol)
            n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
        end
    end
    return nothing
end

"""
    solve_sw!(
        slv::Solver,
        max_threads::Int,
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, Nothing} = nothing,
    )

Solver for the shortwave radiation problem
"""
function solve_sw!(
    slv::Solver,
    max_threads::Int,
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
)
    (; as, op, bcs_sw, src_sw, flux_sw, fluxb_sw) = slv
    context = RTE.context(slv)

    no_args = lookup_sw isa Nothing && lookup_sw_cld isa Nothing
    if no_args
        DA = RTE.array_type(slv)
        FT = RTE.float_type(slv)
        major_gpt2bnd = DA(UInt8[1])
        solar_src_scaled = DA(FT[1])
        flux = flux_sw
    else
        (; major_gpt2bnd, solar_src_scaled) = lookup_sw
        flux = fluxb_sw
    end

    # solving radiative transfer equation
    if slv.op isa OneScalar
        rte_sw_noscat_solve!(
            context,
            flux,
            flux_sw,
            op,
            bcs_sw,
            solar_src_scaled,
            max_threads,
            as,
            lookup_sw,
            lookup_sw_cld,
        ) # no-scattering solver
    else
        rte_sw_2stream_solve!(
            context,
            flux,
            flux_sw,
            op,
            bcs_sw,
            src_sw,
            major_gpt2bnd,
            solar_src_scaled,
            max_threads,
            as,
            lookup_sw,
            lookup_sw_cld,
        ) # 2-stream solver
    end
    return nothing
end

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
        lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
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
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
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
    #-----------------------------------------------------------------------------------
    return nothing
end
#--------------------------------------------------------------------------------------------------
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
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
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
#--------------------------------------------------------------------------------------------------
"""
    rte_sw_2stream_solve!(
        context,
        flux::FluxSW{FT},
        flux_sw::FluxSW{FT},
        op::TwoStream{FT},
        bcs_sw::SwBCs{FT},
        src_sw::SourceSW2Str{FT},
        major_gpt2bnd::AbstractArray{UInt8,1},
        solar_src_scaled::AbstractArray{FT,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lookup_sw::Union{LookUpSW, Nothing} = nothing,
        lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    ) where {FT<:AbstractFloat}

Two stream solver for the shortwave problem.
"""
function rte_sw_2stream_solve!(
    context,
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    major_gpt2bnd::AbstractArray{UInt8, 1},
    solar_src_scaled::AbstractArray{FT, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    (; nlay, ncol) = as
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (
            flux,
            flux_sw,
            op,
            bcs_sw,
            src_sw,
            nlay,
            ncol,
            major_gpt2bnd,
            solar_src_scaled,
            as,
            lookup_sw,
            lookup_sw_cld,
        )
        @info("metadata size",
            flux = sizeof(CUDA.cudaconvert(flux)),
            flux_sw = sizeof(CUDA.cudaconvert(flux_sw)),
            op = sizeof(CUDA.cudaconvert(op)),
            bcs_sw = sizeof(CUDA.cudaconvert(bcs_sw)),
            src_sw = sizeof(CUDA.cudaconvert(src_sw)),
            nlay = sizeof(CUDA.cudaconvert(nlay)),
            ncol = sizeof(CUDA.cudaconvert(ncol)),
            major_gpt2bnd = sizeof(CUDA.cudaconvert(major_gpt2bnd)),
            solar_src_scaled = sizeof(CUDA.cudaconvert(solar_src_scaled)),
            as = sizeof(CUDA.cudaconvert(as)),
            lookup_sw = sizeof(CUDA.cudaconvert(lookup_sw)),
            lookup_sw_cld = sizeof(CUDA.cudaconvert(lookup_sw_cld)),
        )
        @cuda threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(args...)
    else
        @inbounds begin
            bld_cld_mask = as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
            # setting references for flux_sw
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                    bld_cld_mask && Optics.build_cloud_mask!(
                        as.cld_mask_sw,
                        as.cld_frac,
                        as.random_sw,
                        gcol,
                        igpt,
                        as.cld_mask_type,
                    )
                    compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
                    sw_two_stream!(op, src_sw, bcs_sw, gcol, nlay) # Cell properties: transmittance and reflectance for direct and diffuse radiation
                    sw_source_2str!(src_sw, bcs_sw, gcol, flux, igpt, solar_src_scaled, major_gpt2bnd, nlay) # Direct-beam and source for diffuse radiation
                    adding_sw!(src_sw, bcs_sw, gcol, flux, igpt, major_gpt2bnd, nlev) # Transport
                    n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
                end
            end
        end
    end
    return nothing
end
#--------------------------------------------------------------------------------------------------
function rte_sw_2stream_solve_CUDA!(
    flux::FluxSW{FT},
    flux_sw::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    nlay,
    ncol,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    solar_src_scaled::AbstractArray{FT, 1},
    as::AbstractAtmosphericState{FT},
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        @inbounds for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
            if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
                Optics.build_cloud_mask!(as.cld_mask_sw, as.cld_frac, as.random_sw, gcol, igpt, as.cld_mask_type)
            end
            compute_optical_props!(op, as, gcol, igpt, lookup_sw, lookup_sw_cld)
            sw_two_stream!(op, src_sw, bcs_sw, gcol, nlay) # Cell properties: transmittance and reflectance for direct and diffuse radiation
            sw_source_2str!(src_sw, bcs_sw, gcol, flux, igpt, solar_src_scaled, major_gpt2bnd, nlay) # Direct-beam and source for diffuse radiation
            adding_sw!(src_sw, bcs_sw, gcol, flux, igpt, major_gpt2bnd, nlev) # Transport
            n_gpt > 1 && add_to_flux!(flux_sw, flux, gcol)
        end
    end
    return nothing
end
#--------------------------------------------------------------------------------------------------
end
