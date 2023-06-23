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
        lkp_args...
    )

Solver for the longwave radiation problem
"""
function solve_lw!(slv::Solver, max_threads::Int, lkp_args...)
    (; as, op, bcs_lw, src_lw, flux_lw, fluxb_lw) = slv
    (; nlay, ncol) = as
    context = RTE.context(slv)
    DA = RTE.array_type(slv)
    nargs = length(lkp_args)
    @assert nargs < 3
    if nargs > 0
        @assert lkp_args[1] isa LookUpLW
        if nargs == 2
            @assert lkp_args[2] isa LookUpCld
        end
        (; n_gpt, major_gpt2bnd) = lkp_args[1]
        flux = fluxb_lw
    else
        n_gpt = 1
        major_gpt2bnd = DA(UInt8[1])
        flux = flux_lw
    end

    if op isa OneScalar
        rte_lw_noscat_solve!(
            context,
            flux,
            flux_lw,
            src_lw,
            bcs_lw,
            op,
            nlay,
            ncol,
            major_gpt2bnd,
            max_threads,
            as,
            lkp_args...,
        ) # no-scattering solver
    else
        rte_lw_2stream_solve!(
            context,
            flux,
            flux_lw,
            src_lw,
            bcs_lw,
            op,
            nlay,
            ncol,
            major_gpt2bnd,
            max_threads,
            as,
            lkp_args...,
        ) # 2-stream solver
    end
    return nothing
end

"""
    rte_lw_noscat_solve!(
        context,
        flux::FluxLW{FT},
        flux_lw::FluxLW{FT},
        src_lw::SourceLWNoScat{FT},
        bcs_lw::LwBCs{FT},
        op::OneScalar{FT},
        nlay,
        ncol,
        major_gpt2bnd::AbstractArray{UInt8,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lkp_args...,
    ) where {FT<:AbstractFloat}

No scattering solver for the longwave problem.
LW fluxes, no scattering, mu (cosine of integration angle) specified by column
Does radiation calculation at user-supplied angles; converts radiances to flux
using user-supplied weights. Currently, only n_gauss_angles = 1 supported.
"""
function rte_lw_noscat_solve!(
    context,
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
    nlay,
    ncol,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lkp_args...,
) where {FT <: AbstractFloat}
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, flux_lw, src_lw, bcs_lw, op, nlay, ncol, major_gpt2bnd, as, lkp_args...)
        @cuda threads = (tx) blocks = (bx) rte_lw_noscat_solve_CUDA!(args...)
    else
        @inbounds begin
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
                    compute_optical_props!(op, as, src_lw, gcol, igpt, lkp_args...)
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
    lkp_args...,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
            compute_optical_props!(op, as, src_lw, gcol, igpt, lkp_args...)
            rte_lw_noscat_source!(src_lw, op, gcol, nlay)
            rte_lw_noscat_transport!(src_lw, bcs_lw, op, gcol, flux, igpt, major_gpt2bnd, nlay, nlev)
            n_gpt > 1 && add_to_flux!(flux_lw, flux, gcol)
        end
    end
    return nothing
end

"""
    rte_lw_2stream_solve!(
        context,
        flux::FluxLW{FT},
        flux_lw::FluxLW{FT},
        src_lw::SourceLW2Str{FT},
        bcs_lw::LwBCs{FT},
        op::TwoStream{FT},
        nlay,
        ncol,
        major_gpt2bnd::AbstractArray{UInt8,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lkp_args...,
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
function rte_lw_2stream_solve!(
    context,
    flux::FluxLW{FT},
    flux_lw::FluxLW{FT},
    src_lw::SourceLW2Str{FT},
    bcs_lw::LwBCs{FT},
    op::TwoStream{FT},
    nlay,
    ncol,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lkp_args...,
) where {FT <: AbstractFloat}
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, flux_lw, src_lw, bcs_lw, op, nlay, ncol, major_gpt2bnd, as, lkp_args...)
        @cuda threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(args...)
    else
        if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
            @inbounds begin
                ClimaComms.@threaded device for gcol in 1:ncol
                    Optics.build_cloud_mask!(as.cld_mask_lw, as.cld_frac, as.random_lw, gcol, as.cld_mask_type)
                end
            end
        end
        @inbounds begin
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
                    compute_optical_props!(op, as, src_lw, gcol, igpt, lkp_args...)
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
    lkp_args...,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
            Optics.build_cloud_mask!(as.cld_mask_lw, as.cld_frac, as.random_lw, gcol, as.cld_mask_type)
        end
        for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
            compute_optical_props!(op, as, src_lw, gcol, igpt, lkp_args...)
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
        lkp_args...,
    )

Solver for the shortwave radiation problem
"""
function solve_sw!(slv::Solver, max_threads::Int, lkp_args...)
    (; as, op, bcs_sw, src_sw, flux_sw, fluxb_sw) = slv
    (; nlay, ncol) = as
    context = RTE.context(slv)
    FT = RTE.float_type(slv)
    DA = RTE.array_type(slv)
    nargs = length(lkp_args)
    @assert nargs < 3
    if nargs > 0
        @assert lkp_args[1] isa LookUpSW
        if nargs == 2
            @assert lkp_args[2] isa LookUpCld
        end
        (; n_gpt, major_gpt2bnd, solar_src_scaled) = lkp_args[1]
        flux = fluxb_sw
    else
        n_gpt = 1
        major_gpt2bnd = DA(UInt8[1])
        solar_src_scaled = DA(FT[1])
        flux = flux_sw
    end
    # solving radiative transfer equation
    if slv.op isa OneScalar
        rte_sw_noscat_solve!(
            context,
            flux,
            flux_sw,
            op,
            bcs_sw,
            nlay,
            ncol,
            solar_src_scaled,
            max_threads,
            as,
            lkp_args...,
        ) # no-scattering solver
    else
        rte_sw_2stream_solve!(
            context,
            flux,
            flux_sw,
            op,
            bcs_sw,
            src_sw,
            nlay,
            ncol,
            major_gpt2bnd,
            solar_src_scaled,
            max_threads,
            as,
            lkp_args...,
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
        nlay,
        ncol,
        solar_src_scaled::AbstractArray{FT,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lkp_args...,
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
    nlay,
    ncol,
    solar_src_scaled::AbstractArray{FT, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lkp_args...,
) where {FT <: AbstractFloat}
    nlev = nlay + 1
    n_gpt = length(solar_src_scaled)
    device = ClimaComms.device(context)
    # setting references for flux_sw
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, flux_sw, op, bcs_sw, nlay, ncol, solar_src_scaled, as, lkp_args...)
        @cuda threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    else
        @inbounds begin
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                    compute_optical_props!(op, as, gcol, igpt, lkp_args...)
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
    lkp_args...,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(solar_src_scaled)
    # setting references for flux_sw
    if gcol ≤ ncol
        for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol, igpt, lkp_args...)
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
        nlay,
        ncol,
        major_gpt2bnd::AbstractArray{UInt8,1},
        solar_src_scaled::AbstractArray{FT,1},
        max_threads,
        as::AbstractAtmosphericState{FT},
        lkp_args...,
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
    nlay,
    ncol,
    major_gpt2bnd::AbstractArray{UInt8, 1},
    solar_src_scaled::AbstractArray{FT, 1},
    max_threads,
    as::AbstractAtmosphericState{FT},
    lkp_args...,
) where {FT <: AbstractFloat}
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, flux_sw, op, bcs_sw, src_sw, nlay, ncol, major_gpt2bnd, solar_src_scaled, as, lkp_args...)
        @cuda threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(args...)
    else
        @inbounds begin
            if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
                ClimaComms.@threaded device for gcol in 1:ncol
                    Optics.build_cloud_mask!(as.cld_mask_sw, as.cld_frac, as.random_sw, gcol, as.cld_mask_type)
                end
            end
            # setting references for flux_sw
            for igpt in 1:n_gpt
                ClimaComms.@threaded device for gcol in 1:ncol
                    igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
                    compute_optical_props!(op, as, gcol, igpt, lkp_args...)
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
    lkp_args...,
) where {FT <: AbstractFloat}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(major_gpt2bnd)
    if gcol ≤ ncol
        if as isa AtmosphericState && as.cld_mask_type isa AbstractCloudMask
            Optics.build_cloud_mask!(as.cld_mask_sw, as.cld_frac, as.random_sw, gcol, as.cld_mask_type)
        end
        for igpt in 1:n_gpt
            igpt == 1 && set_flux_to_zero!(flux_sw, gcol)
            compute_optical_props!(op, as, gcol, igpt, lkp_args...)
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
