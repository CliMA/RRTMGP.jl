module RTESolver

using StaticArrays
using CUDA
using Adapt
using UnPack

using ..Device: array_type, array_device, CUDADevice, CPU
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
# -------------------------------------------------------------------------------------------------
# Solver for the longwave radiation problem
# -------------------------------------------------------------------------------------------------
function solve_lw!(
    slv::Solver{FT,I},
    max_threads::I,
    lkp_args...,
) where {I<:Int,FT<:AbstractFloat}
    @unpack as, op, bcs_lw, src_lw, flux_lw, fluxb_lw = slv
    @unpack nlay, ncol = as

    nargs = length(lkp_args)
    @assert nargs < 3
    if nargs > 0
        @assert lkp_args[1] isa LookUpLW
        if nargs == 2
            @assert lkp_args[2] isa LookUpCld
        end
        n_gpt = lkp_args[1].n_gpt
        flux = fluxb_lw
    else
        n_gpt = 1
        flux = flux_lw
    end
    ibnd = 1
    set_flux_to_zero!(flux_lw)

    for igpt = 1:n_gpt
        if nargs > 0 && lkp_args[1] isa LookUpLW
            ibnd = I(lkp_args[1].major_gpt2bnd[igpt])
        end
        # computing optical properties
        compute_optical_props!(op, as, src_lw, igpt, lkp_args...)

        if op isa OneScalar
            rte_lw_noscat_solve!(
                flux,
                src_lw,
                bcs_lw,
                op,
                nlay,
                ncol,
                igpt,
                ibnd,
                max_threads,
            ) # no-scattering solver
        else
            rte_lw_2stream_solve!(
                flux,
                src_lw,
                bcs_lw,
                op,
                nlay,
                ncol,
                igpt,
                ibnd,
                max_threads,
            ) # 2-stream solver
        end

        nargs > 0 && add_to_flux!(flux_lw, flux)
    end
    flux_lw.flux_net .= flux_lw.flux_up .- flux_lw.flux_dn
    return nothing
end
# -------------------------------------------------------------------------------------------------
#
# LW fluxes, no scattering, mu (cosine of integration angle) specified by column
# Does radiation calculation at user-supplied angles; converts radiances to flux
# using user-supplied weights
#
# ---------------------------------------------------------------    
function rte_lw_noscat_solve!(
    flux::FluxLW{FT},
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
    nlay::I,
    ncol::I,
    igpt::I,
    ibnd::I,
    max_threads::I,
) where {FT<:AbstractFloat,I<:Int}
    nlev = nlay + 1
    device = array_device(op.τ)
    if device === CUDADevice() # launch CUDA kernel
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, src_lw, bcs_lw, op, nlay, ncol, igpt, ibnd)
        @cuda threads = (tx) blocks = (bx) rte_lw_noscat_solve_CUDA!(args...)
    else # use Julia native multithreading
        @inbounds Threads.@threads for gcol = 1:ncol
            rte_lw_noscat_source!(src_lw, op, gcol, nlay)
            rte_lw_noscat_transport!(
                src_lw,
                bcs_lw,
                op,
                gcol,
                flux,
                igpt,
                ibnd,
                nlay,
                nlev,
            )
        end
    end
    #------------------------------------------------
    return nothing
end
# ---------------------------------------------------------------    
function rte_lw_noscat_solve_CUDA!(
    flux::FluxLW{FT},
    src_lw::SourceLWNoScat{FT},
    bcs_lw::LwBCs{FT},
    op::OneScalar{FT},
    nlay::I,
    ncol::I,
    igpt::I,
    ibnd::I,
) where {FT<:AbstractFloat,I<:Int}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    if gcol ≤ ncol
        rte_lw_noscat_source!(src_lw, op, gcol, nlay)
        rte_lw_noscat_transport!(
            src_lw,
            bcs_lw,
            op,
            gcol,
            flux,
            igpt,
            ibnd,
            nlay,
            nlev,
        )
    end
    return nothing
end
#-------------------------------------------------------------------------------------------------
# Longwave two-stream calculation:
# RRTMGP provides source functions at each level using the spectral mapping
# of each adjacent layer. Combine these for two-stream calculations
# lw combine sources (rte_lw_2stream_combine_sources!)
# combine RRTMGP-specific sources at levels
# compute layer reflectance, transmittance
# compute total source function at levels using linear-in-tau (rte_lw_2stream_source!)
# transport (adding_lw!)
#-------------------------------------------------------------------------------------------------
function rte_lw_2stream_solve!(
    flux::FluxLW{FT},
    src_lw::SourceLW2Str{FT},
    bcs_lw::LwBCs{FT},
    op::TwoStream{FT},
    nlay::I,
    ncol::I,
    igpt::I,
    ibnd::I,
    max_threads::I,
) where {FT<:AbstractFloat,I<:Int}
    nlev = nlay + 1
    device = array_device(op.τ)
    if device === CUDADevice() # launch CUDA kernel
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, src_lw, bcs_lw, op, nlay, ncol, igpt, ibnd)
        @cuda threads = (tx) blocks = (bx) rte_lw_2stream_solve_CUDA!(args...)
    else # use Julia native multithreading
        Threads.@threads for gcol = 1:ncol
            rte_lw_2stream_combine_sources!(src_lw, gcol, nlev, ncol)
            rte_lw_2stream_source!(op, src_lw, gcol, nlay, ncol)
            adding_lw!(flux, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
        end
    end
    #------------------------------------------------
    return nothing
end
# -------------------------------------------------------------------------------------------------
function rte_lw_2stream_solve_CUDA!(
    flux::FluxLW{FT},
    src_lw::SourceLW2Str{FT},
    bcs_lw::LwBCs{FT},
    op::TwoStream{FT},
    nlay::I,
    ncol::I,
    igpt::I,
    ibnd::I,
) where {FT<:AbstractFloat,I<:Int}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    if gcol ≤ ncol
        rte_lw_2stream_combine_sources!(src_lw, gcol, nlev, ncol)
        rte_lw_2stream_source!(op, src_lw, gcol, nlay, ncol)
        adding_lw!(flux, src_lw, bcs_lw, gcol, igpt, ibnd, nlev, ncol)
    end
end
# -------------------------------------------------------------------------------------------------
# Solver for the shortwave radiation problem
# -------------------------------------------------------------------------------------------------
function solve_sw!(
    slv::Solver{FT,I},
    max_threads::I,
    lkp_args...,
) where {I<:Int,FT<:AbstractFloat}
    @unpack as, op, bcs_sw, src_sw, flux_sw, fluxb_sw = slv
    @unpack nlay, ncol = as

    nargs = length(lkp_args)
    @assert nargs < 3
    if nargs > 0
        @assert lkp_args[1] isa LookUpSW
        if nargs == 2
            @assert lkp_args[2] isa LookUpCld
        end
        n_gpt = lkp_args[1].n_gpt
        flux = fluxb_sw
    else
        n_gpt = 1
        flux = flux_sw
    end
    ibnd = 1
    solar_frac = FT(1)
    set_flux_to_zero!(slv.flux_sw)
    for igpt = 1:n_gpt
        if nargs > 0 && lkp_args[1] isa LookUpSW
            ibnd = I(lkp_args[1].major_gpt2bnd[igpt])
            solar_frac = lkp_args[1].solar_src_scaled[igpt]
        end
        # computing optical properties
        compute_optical_props!(op, as, igpt, lkp_args...)

        # solving radiative transfer equation
        if slv.op isa OneScalar
            rte_sw_noscat_solve!(
                flux,
                op,
                bcs_sw,
                nlay,
                ncol,
                igpt,
                ibnd,
                solar_frac,
                max_threads,
            ) # no-scattering solver
        else
            rte_sw_2stream_solve!(
                flux,
                op,
                bcs_sw,
                src_sw,
                nlay,
                ncol,
                igpt,
                ibnd,
                solar_frac,
                max_threads,
            ) # 2-stream solver
            flux.flux_dn .+= flux.flux_dn_dir
        end
        nargs > 0 && add_to_flux!(flux_sw, flux)
    end
    if slv.op isa TwoStream
        flux_sw.flux_net .= flux_sw.flux_up .- flux_sw.flux_dn
    end

    return nothing
end
#--------------------------------------------------------------------------------------------------
# Extinction-only i.e. solar direct beam
function rte_sw_noscat_solve!(
    flux::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    nlay::I,
    ncol::I,
    igpt::I,
    ibnd::I,
    solar_frac::FT,
    max_threads::I,
) where {FT<:AbstractFloat,I<:Int}
    nlev = nlay + 1
    device = array_device(op.τ)
    if device === CUDADevice() # launching CUDA kernel
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, op, bcs_sw, nlay, ncol, igpt, ibnd, solar_frac)
        @cuda threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    else # launch juila native multithreading
        @inbounds Threads.@threads for gcol = 1:ncol
            rte_sw_noscat_solve_kernel!(
                flux,
                op,
                bcs_sw,
                solar_frac,
                gcol,
                nlev,
            )
        end
    end
    #-----------------------------------------------------------------------------------
    return nothing
end
#--------------------------------------------------------------------------------------------------
function rte_sw_noscat_solve_CUDA!(
    flux::FluxSW{FT},
    op::OneScalar{FT},
    bcs_sw::SwBCs{FT},
    nlay::I,
    ncol::I,
    igpt::I,
    ibnd::I,
    solar_frac::FT,
) where {FT<:AbstractFloat,I<:Int}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    if gcol ≤ ncol
        rte_sw_noscat_solve_kernel!(flux, op, bcs_sw, solar_frac, gcol, nlev)
    end
    return nothing
end
#--------------------------------------------------------------------------------------------------
function rte_sw_2stream_solve!(
    flux::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    nlay::I,
    ncol::I,
    igpt::I,
    ibnd::I,
    solar_frac::FT,
    max_threads::I,
) where {FT<:AbstractFloat,I<:Int}
    nlev = nlay + 1
    device = array_device(op.τ)
    if device === CUDADevice()
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        args = (flux, op, bcs_sw, src_sw, nlay, ncol, igpt, ibnd, solar_frac)
        @cuda threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(args...)
    else # launch julia native multithreading
        @inbounds Threads.@threads for gcol = 1:ncol
            sw_two_stream!(op, src_sw, bcs_sw, gcol, nlay) # Cell properties: transmittance and reflectance for direct and diffuse radiation
            sw_source_2str!(src_sw, bcs_sw, gcol, flux, solar_frac, ibnd, nlay) # Direct-beam and source for diffuse radiation
            adding_sw!(src_sw, bcs_sw, gcol, flux, igpt, ibnd, nlev) # Transport
        end
    end
    #--------------------------------------
    return nothing
end
#--------------------------------------------------------------------------------------------------
function rte_sw_2stream_solve_CUDA!(
    flux::FluxSW{FT},
    op::TwoStream{FT},
    bcs_sw::SwBCs{FT},
    src_sw::SourceSW2Str{FT},
    nlay::I,
    ncol::I,
    igpt::I,
    ibnd::I,
    solar_frac::FT,
) where {FT<:AbstractFloat,I<:Int}
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    if gcol ≤ ncol
        sw_two_stream!(op, src_sw, bcs_sw, gcol, nlay) # Cell properties: transmittance and reflectance for direct and diffuse radiation
        sw_source_2str!(src_sw, bcs_sw, gcol, flux, solar_frac, ibnd, nlay) # Direct-beam and source for diffuse radiation
        adding_sw!(src_sw, bcs_sw, gcol, flux, igpt, ibnd, nlev) # Transport
    end
end
#--------------------------------------------------------------------------------------------------
end
