module GrayRTESolver

using StaticArrays
using KernelAbstractions
using CUDA
using Adapt

using ..Device: array_type, array_device
using ..GrayAngularDiscretizations
using ..GrayFluxes
using ..GrayBCs
using ..GraySources
using ..GrayAtmos
using ..GrayOptics

export gray_atmos_lw!,
    rte_lw_noscat_gray_solve!, rte_lw_2stream_gray_solve!, adding_gray!


include("GrayRTESolverKernels.jl")
# -------------------------------------------------------------------------------------------------
# Solver for the longwave gray radiation problem
# -------------------------------------------------------------------------------------------------

function gray_atmos_lw!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D};
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}
    # setting references
    flux = gray_rrtmgp.flux
    flux_up = flux.flux_up               # upward flux
    flux_dn = flux.flux_dn               # downward flux
    flux_net = flux.flux_net             # net flux

    as = gray_rrtmgp.as                  # gray atmospheric state
    optical_props = gray_rrtmgp.op       # optical properties
    source = gray_rrtmgp.src             # Planck sources
    nlay = as.nlay                       # number of layers per column
    # computing optical properties
    gas_optics_gray_atmos!(as, optical_props, source)
    # solving radiative transfer equation
    if typeof(gray_rrtmgp.op) == GrayOneScalar{FT,FTA2D}
        rte_lw_noscat_gray_solve!(gray_rrtmgp, max_threads = max_threads) # no-scattering solver
    else
        rte_lw_2stream_gray_solve!(gray_rrtmgp, max_threads = max_threads) # 2-stream solver
    end
    return nothing
end
# -------------------------------------------------------------------------------------------------
#
# LW fluxes, no scattering, mu (cosine of integration angle) specified by column
#   Does radiation calculation at user-supplied angles; converts radiances to flux
#   using user-supplied weights
#
# ---------------------------------------------------------------    
function rte_lw_noscat_gray_solve!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D};
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
}
    # Setting references 
    τ = gray_rrtmgp.op.τ                   # Optical thickness

    ncol = size(gray_rrtmgp.flux.flux_up, 2)         # number of columns
    nlev = size(gray_rrtmgp.flux.flux_up, 1)         # number of levels
    nlay = nlev - 1                 # number of layers
    ngpt = 1                        # number of g-points (only 1 for gray radiation model

    device = array_device(gray_rrtmgp.flux.flux_up)
    comp_stream = Event(device)
    #------Launching source computation kernel-----------------------
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlay)

    workgroup = (thr_x, thr_y)
    ndrange = (nlay, ncol)

    comp_stream = rte_lw_noscat_gray_source_kernel!(device, workgroup)(
        gray_rrtmgp,
        Val(nlay),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    #------Launching transport computation kernel-----------------------
    thr_x = min(max_threads, ncol)

    workgroup = (thr_x)
    ndrange = (ncol)

    comp_stream = rte_lw_noscat_gray_transport_kernel!(device, workgroup)(
        gray_rrtmgp,
        Val(nlay),
        Val(nlev),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    wait(comp_stream)
    #------------------------------------------------
    return nothing
end

function rte_lw_2stream_gray_solve!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D};
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
}
    # Setting references
    flux_up = gray_rrtmgp.flux.flux_up
    flux_dn = gray_rrtmgp.flux.flux_dn
    flux_net = gray_rrtmgp.flux.flux_net
    source = gray_rrtmgp.src
    sfc_emis = gray_rrtmgp.bcs.sfc_emis
    inc_flux = gray_rrtmgp.bcs.inc_flux
    ncol, nlev, ngpt = size(flux_up, 2), size(flux_up, 1), 1
    nlay = nlev - 1
    #----------------------------------------------------------------------------------
    #    if inc_flux ≠ nothing
    #        flux_dn[i_lev_top, 1] = inc_flux[1]
    #    end
    #-----------------------------------------------------------------------------------
    device = array_device(gray_rrtmgp.src.lev_source)
    comp_stream = Event(device)

    #------Launching source computation kernel-----------------------
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlev)

    workgroup = (thr_x, thr_y)
    ndrange = (nlev, ncol)

    comp_stream =
        rte_lw_2stream_gray_combine_sources_kernel!(device, workgroup)(
            gray_rrtmgp.src,
            Val(nlev),
            Val(ncol),
            ndrange = ndrange,
            dependencies = (comp_stream,),
        )
    # Cell properties: reflection, transmission for diffuse radiation
    # Coupling coefficients needed for source function
    # -------------------------------------------------------------------------------------------------
    #
    # Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
    #    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
    #
    # Equations are developed in Meador and Weaver, 1980,
    #    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
    #
    # -------------------------------------------------------------------------------------------------    
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlay)

    workgroup = (thr_x, thr_y)
    ndrange = (nlay, ncol)

    comp_stream = rte_lw_2stream_gray_source_kernel!(device, workgroup)(
        gray_rrtmgp,
        Val(nlay),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    # -------------------------------------------------------------------------------------------------    
    # Transport
    thr_x = min(max_threads, ncol)

    workgroup = (thr_x)
    ndrange = (ncol)

    comp_stream = adding_gray_kernel!(device, workgroup)(
        gray_rrtmgp,
        Val(nlev),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    wait(comp_stream)
    #------------------------------------------------
    return nothing
end

end
