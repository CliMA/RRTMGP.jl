module GrayRTESolver

using ..Device: array_type, array_device
using ..GrayAngularDiscretizations
using ..GraySources
using ..GrayAtmos
using ..GrayOptics

export gray_atmos_lw!,
    rte_lw_noscat_gray_solve!, rte_lw_2stream_gray_solve!, adding_gray!

using KernelAbstractions
using CUDA


include("GrayRTESolverKernels.jl")
# -------------------------------------------------------------------------------------------------
# Solver for the longwave gray radiation problem
# -------------------------------------------------------------------------------------------------

function gray_atmos_lw!(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D,FTA3D,B};
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    B<:Bool,
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
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D,FTA3D,B};
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    I<:Int,
    B<:Bool,
}
    # Setting references 
    flux_up = gray_rrtmgp.flux.flux_up     # upward flux
    flux_dn = gray_rrtmgp.flux.flux_dn     # downward flux
    flux_net = gray_rrtmgp.flux.flux_net   # net flux
    top_at_1 = gray_rrtmgp.as.top_at_1     # is top-of-atmos located at point # 1?
    source = gray_rrtmgp.src               # Planck sources
    sfc_emis = gray_rrtmgp.bcs.sfc_emis    # Surface emissivity
    inc_flux = gray_rrtmgp.bcs.inc_flux    # Surface flux
    τ = gray_rrtmgp.op.τ                   # Optical thickness

    ncol = size(flux_up, 2)         # number of columns
    nlev = size(flux_up, 1)         # number of levels
    nlay = nlev - 1                 # number of layers
    ngpt = 1                        # number of g-points (only 1 for gray radiation model

    i_lev_top = top_at_1 ? 1 : nlev # index for top level of column
    n_sfc = top_at_1 ? nlev : 1     # index for bottom level of column
    n_μ = gray_rrtmgp.angle_disc.n_gauss_angles
    Ds = gray_rrtmgp.angle_disc.gauss_Ds
    w_μ = gray_rrtmgp.angle_disc.gauss_wts

    lev_source_up = top_at_1 ? source.lev_source_dec : source.lev_source_inc # Mapping increasing/decreasing indices to up/down
    lev_source_dn = top_at_1 ? source.lev_source_inc : source.lev_source_dec
    lay_source = source.lay_source # Planck source at average layer temperature [W/m^2]
    sfc_source = source.sfc_source # Surface source function [W/m^2]
    source_up = source.source_up
    source_dn = source.source_dn
    trans = source.trans

    device = array_device(flux_up)
    comp_stream = Event(device)
    #------Launching source computation kernel-----------------------
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlay)

    workgroup = (thr_x, thr_y)
    ndrange = (nlay, ncol)

    comp_stream = rte_lw_noscat_gray_source_kernel!(device, workgroup)(
        lay_source,
        lev_source_up,
        lev_source_dn,
        τ,
        trans,
        source_up,
        source_dn,
        Ds,
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
        flux_up,
        flux_dn,
        flux_net,
        trans,
        source_up,
        source_dn,
        sfc_source,
        sfc_emis,
        inc_flux,
        τ,
        n_μ,
        Ds,
        w_μ,
        top_at_1,
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
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D,FTA3D,B};
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    I<:Int,
    B<:Bool,
}
    # Setting references
    flux_up = gray_rrtmgp.flux.flux_up
    flux_dn = gray_rrtmgp.flux.flux_dn
    flux_net = gray_rrtmgp.flux.flux_net
    top_at_1 = gray_rrtmgp.as.top_at_1
    source = gray_rrtmgp.src
    sfc_emis = gray_rrtmgp.bcs.sfc_emis
    inc_flux = gray_rrtmgp.bcs.inc_flux
    ncol, nlev, ngpt = size(flux_up, 2), size(flux_up, 1), 1
    nlay = nlev - 1
    lev_src_inc, lev_src_dec = source.lev_source_inc, source.lev_source_dec
    lev_source = source.lev_source
    lay_source, sfc_source = source.lay_source, source.sfc_source
    i_lev_top = top_at_1 ? 1 : nlev
    n_sfc = top_at_1 ? nlev : 1

    τ = gray_rrtmgp.op.τ
    g = gray_rrtmgp.op.g
    ssa = gray_rrtmgp.op.ssa

    Rdif, Tdif = source.Rdif, source.Tdif
    src_up, src_dn = source.src_up, source.src_dn # Source function for diffuse radiation
    albedo, src = source.albedo, source.src
    #----------------------------------------------------------------------------------
    if inc_flux ≠ nothing
        flux_dn[i_lev_top, 1] = inc_flux[1]
    end
    #-----------------------------------------------------------------------------------
    lw_diff_sec = FT(1.66)
    #-----------------------------------------------------------------------------------
    device = array_device(lev_source)
    comp_stream = Event(device)

    #------Launching source computation kernel-----------------------
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlev)

    workgroup = (thr_x, thr_y)
    ndrange = (nlev, ncol)

    comp_stream =
        rte_lw_2stream_gray_combine_sources_kernel!(device, workgroup)(
            lev_src_inc,
            lev_src_dec,
            lev_source,
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
        τ,
        ssa,
        g,
        lev_source,
        src_up,
        src_dn,
        Tdif,
        Rdif,
        lw_diff_sec,
        top_at_1,
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
        top_at_1,
        sfc_emis,
        Rdif,
        Tdif,
        src_up,
        src_dn,
        sfc_source,
        flux_up,
        flux_dn,
        flux_net,
        albedo,
        src,
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
