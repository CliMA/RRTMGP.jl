module RTESolver

using StaticArrays
using KernelAbstractions
using CUDA
using Adapt
using UnPack

using ..Device: array_type, array_device
using ..AngularDiscretizations
using ..Fluxes
using ..BCs
using ..Sources
using ..RTE
using ..Optics
using ..LookUpTables

export solve_lw!, solve_sw!


include("RTESolverKernels.jl")
# -------------------------------------------------------------------------------------------------
# Solver for the longwave gray radiation problem
# -------------------------------------------------------------------------------------------------
function solve_lw!(
    slv::Solver{FT,I,FTA1D,FTA2D};
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}
    # computing optical properties
    compute_optical_props!(slv.as, slv.op, islw = true, sf = slv.src_lw)
    # solving radiative transfer equation
    if typeof(slv.op) == OneScalar{FT,FTA2D}
        rte_lw_noscat_solve!(slv, isgray = true, max_threads = max_threads) # no-scattering solver
    else
        rte_lw_2stream_solve!(slv, isgray = true, max_threads = max_threads) # 2-stream solver
    end
    slv.flux_lw.flux_net .= slv.flux_lw.flux_up .- slv.flux_lw.flux_dn
    return nothing
end
# -------------------------------------------------------------------------------------------------
# Solver for the longwave radiation problem with full gas optics (RRTMGP)
# -------------------------------------------------------------------------------------------------
function solve_lw!(
    slv::Solver{FT,I,FTA1D,FTA2D},
    lkp::LookUpLW{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D};
    lkp_cld::Union{LookUpCld{I,B,FT,FTA1D,FTA2D,FTA3D,FTA4D},Nothing} = nothing,
    max_threads::I = Int(256),
) where {
    I<:Int,
    B<:Bool,
    FT<:AbstractFloat,
    UI8<:UInt8,
    UI8A1D<:AbstractArray{UI8,1},
    IA1D<:AbstractArray{I,1},
    IA2D<:AbstractArray{I,2},
    IA3D<:AbstractArray{I,3},
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    FTA4D<:AbstractArray{FT,4},
}
    n_gpt = lkp.n_gpt
    set_flux_to_zero!(slv.flux_lw)
    for igpt = 1:n_gpt
        ibnd = I(lkp.major_gpt2bnd[igpt])
        # computing optical properties
        compute_optical_props!(
            slv.as,
            lkp,
            slv.op,
            igpt,
            islw = true,
            sf = slv.src_lw,
            lkp_cld = lkp_cld,
            max_threads = max_threads,
        )
        # solving radiative transfer equation
        if typeof(slv.op) == OneScalar{FT,FTA2D}
            rte_lw_noscat_solve!(
                slv,
                isgray = false,
                igpt = igpt,
                ibnd = ibnd,
                max_threads = max_threads,
            ) # no-scattering solver
        else
            rte_lw_2stream_solve!(
                slv,
                isgray = false,
                igpt = igpt,
                ibnd = ibnd,
                max_threads = max_threads,
            ) # 2-stream solver
        end
        add_to_flux!(slv.flux_lw, slv.fluxb_lw)
    end
    slv.flux_lw.flux_net .= slv.flux_lw.flux_up .- slv.flux_lw.flux_dn
    return nothing
end
# -------------------------------------------------------------------------------------------------
#
# LW fluxes, no scattering, mu (cosine of integration angle) specified by column
#   Does radiation calculation at user-supplied angles; converts radiances to flux
#   using user-supplied weights
#
# ---------------------------------------------------------------    
function rte_lw_noscat_solve!(
    slv::Solver{FT,I,FTA1D,FTA2D};
    isgray::Bool = false,
    igpt::Int = 1,
    ibnd::Int = 1,
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
}
    # Setting references 
    τ = slv.op.τ         # Optical thickness
    ncol = slv.as.ncol   # number of columns
    nlay = slv.as.nlay
    nlev = nlay + 1      # number of levels

    if isgray
        flux = slv.flux_lw
    else
        flux = slv.fluxb_lw
    end

    device = array_device(slv.op.τ)
    comp_stream = Event(device)
    #------Launching source computation kernel-----------------------
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlay)

    workgroup = (thr_x, thr_y)
    ndrange = (nlay, ncol)

    comp_stream = rte_lw_noscat_source_kernel!(device, workgroup)(
        slv,
        Val(nlay),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    #------Launching transport computation kernel-----------------------
    thr_x = min(max_threads, ncol)

    workgroup = (thr_x)
    ndrange = (ncol)

    comp_stream = rte_lw_noscat_transport_kernel!(device, workgroup)(
        slv,
        flux,
        igpt,
        ibnd,
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

function rte_lw_2stream_solve!(
    slv::Solver{FT,I,FTA1D,FTA2D};
    isgray::Bool = false,
    igpt::Int = 1,
    ibnd::Int = 1,
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
}
    ncol, nlay = slv.as.ncol, slv.as.nlay
    nlev = nlay + 1
    islw = true # for adding kernel

    flux = isgray ? slv.flux_lw : slv.fluxb_lw
    #-----------------------------------------------------------------------------------
    device = array_device(slv.op.τ)
    comp_stream = Event(device)
    #------Launching source computation kernel-----------------------

    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlev)

    workgroup = (thr_x, thr_y)
    ndrange = (nlev, ncol)

    comp_stream = rte_lw_2stream_combine_sources_kernel!(device, workgroup)(
        slv.src_lw,
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

    comp_stream = rte_lw_2stream_source_kernel!(device, workgroup)(
        slv,
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

    comp_stream = adding_kernel!(device, workgroup)(
        slv,
        flux,
        islw,
        igpt,
        ibnd,
        Val(nlev),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    wait(comp_stream)
    #------------------------------------------------
    return nothing
end

# -------------------------------------------------------------------------------------------------
# Solver for the shortwave gray radiation problem
# -------------------------------------------------------------------------------------------------
function solve_sw!(
    slv::Solver{FT,I,FTA1D,FTA2D};
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}
    # computing optical properties
    compute_optical_props!(slv.as, slv.op, islw = false, sf = slv.src_sw)

    # solving radiative transfer equation
    if typeof(slv.op) == OneScalar{FT,FTA2D}
        rte_sw_noscat_solve!(slv, isgray = true, max_threads = max_threads) # no-scattering solver
    else
        rte_sw_2stream_solve!(slv, isgray = true, max_threads = max_threads) # 2-stream solver
    end
    if slv.op isa TwoStream
        slv.flux_sw.flux_net .= slv.flux_sw.flux_up .- slv.flux_sw.flux_dn
    end
    return nothing
end
# -------------------------------------------------------------------------------------------------
# Solver for the shortwave radiation problem with full gas optics (RRTMGP)
# -------------------------------------------------------------------------------------------------
function solve_sw!(
    slv::Solver{FT,I,FTA1D,FTA2D},
    lkp::LookUpSW{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D};
    lkp_cld::Union{LookUpCld{I,B,FT,FTA1D,FTA2D,FTA3D,FTA4D},Nothing} = nothing,
    max_threads::I = Int(256),
) where {
    I<:Int,
    B<:Bool,
    FT<:AbstractFloat,
    UI8<:UInt8,
    UI8A1D<:AbstractArray{UI8,1},
    IA1D<:AbstractArray{I,1},
    IA2D<:AbstractArray{I,2},
    IA3D<:AbstractArray{I,3},
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    FTA4D<:AbstractArray{FT,4},
}
    n_gpt = lkp.n_gpt
    set_flux_to_zero!(slv.flux_sw)
    for igpt = 1:n_gpt
        ibnd = I(lkp.major_gpt2bnd[igpt])
        # computing optical properties
        compute_optical_props!(
            slv.as,
            lkp,
            slv.op,
            igpt,
            islw = false,
            sf = slv.src_sw,
            lkp_cld = lkp_cld,
            max_threads = max_threads,
        )

        solar_frac = lkp.solar_src_scaled[igpt]
        # solving radiative transfer equation
        if typeof(slv.op) == OneScalar{FT,FTA2D}
            rte_sw_noscat_solve!(
                slv,
                isgray = false,
                solar_frac = solar_frac,
                max_threads = max_threads,
            ) # no-scattering solver
        else
            rte_sw_2stream_solve!(
                slv,
                isgray = false,
                igpt = igpt,
                ibnd = ibnd,
                solar_frac = solar_frac,
                max_threads = max_threads,
            ) # 2-stream solver
            slv.fluxb_sw.flux_dn .+= slv.fluxb_sw.flux_dn_dir
        end
        add_to_flux!(slv.flux_sw, slv.fluxb_sw)
    end
    if slv.op isa TwoStream
        slv.flux_sw.flux_net .= slv.flux_sw.flux_up .- slv.flux_sw.flux_dn
    end

    return nothing
end
#--------------------------------------------------------------------------------------------------
function rte_sw_noscat_solve!(
    slv::Solver{FT,I,FTA1D,FTA2D};
    isgray::Bool = false,
    solar_frac::FT = FT(1),
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
}
    nlev, ncol = slv.as.nlay + 1, slv.as.ncol
    flux = (isgray ? slv.flux_sw : slv.fluxb_sw)
    #-----------------------------------------------------------------------------------
    device = array_device(slv.op.τ)
    comp_stream = Event(device)
    thr_x = min(max_threads, ncol)
    workgroup = (thr_x)
    ndrange = (ncol)

    comp_stream = rte_sw_noscat_solve_kernel!(device, workgroup)(
        slv,
        flux,
        solar_frac,
        Val(nlev),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    wait(comp_stream)
    return nothing
end
#--------------------------------------------------------------------------------------------------
function rte_sw_2stream_solve!(
    slv::Solver{FT,I,FTA1D,FTA2D};
    isgray::Bool = false,
    igpt::Int = 1,
    ibnd::Int = 1,
    solar_frac::FT = FT(1),
    max_threads::I = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
}
    nlay, ncol = slv.as.nlay, slv.as.ncol
    nlev = nlay + 1
    flux = (isgray ? slv.flux_sw : slv.fluxb_sw)
    islw = false
    #--------------------------------------
    device = array_device(slv.op.τ)
    comp_stream = Event(device)
    #----launching sw_two_stream_kernel----------------------------------
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlay)

    workgroup = (thr_x, thr_y)
    ndrange = (nlay, ncol)

    comp_stream = sw_two_stream_kernel!(device, workgroup)(
        slv,
        Val(nlay),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    #--------------------------------------
    thr_x = min(max_threads, ncol)
    workgroup = (thr_x)
    ndrange = (ncol)

    comp_stream = sw_source_2str_kernel!(device, workgroup)(
        slv,
        flux,
        solar_frac,
        ibnd,
        Val(nlay),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    # -------------------------------------------------------------------------------------------------    
    # adding function only computes diffuse flux
    thr_x = min(max_threads, ncol)

    workgroup = (thr_x)
    ndrange = (ncol)
    comp_stream = adding_kernel!(device, workgroup)(
        slv,
        flux,
        islw,
        igpt,
        ibnd,
        Val(nlev),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    wait(comp_stream)
    #--------------------------------------
    return nothing
end
#--------------------------------------------------------------------------------------------------
end
