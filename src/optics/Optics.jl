module Optics

using DocStringExtensions
using KernelAbstractions
using CUDA
using ..Device: array_type, array_device
using Adapt
#---------------------------------------
using ..Vmrs
using ..LookUpTables
using ..AtmosphericStates
using ..Sources
#---------------------------------------
using CLIMAParameters
using CLIMAParameters.Planet: molmass_dryair, molmass_water, grav
#---------------------------------------

export AbstractOpticalProps,
    OneScalar,
    TwoStream,
    init_optical_props,
    compute_col_dry!,
    compute_optical_props!


abstract type AbstractOpticalProps{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} end

"""
    OneScalar{FT    <: AbstractFloat,
                  FTA2D <: AbstractArray{FT,2}} <: AbstractOpticalProps{FT,FTA2D

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OneScalar{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractOpticalProps{FT,FTA2D}
    "Optical Depth"
    τ::FTA2D
end
Adapt.@adapt_structure OneScalar

function OneScalar(
    ::Type{FT},
    ncol::Int,
    nlay::Int,
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return OneScalar{FT,DA{FT,2}}(DA{FT,2}(undef, nlay, ncol))
end

"""
    TwoStream{FT<:AbstractFloat,
                     FTA2D<:AbstractArray{FT,2}} <: AbstractOpticalProps{FT,FTA2D}

Two stream approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStream{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractOpticalProps{FT,FTA2D}
    "Optical depth"
    τ::FTA2D
    "Single-scattering albedo"
    ssa::FTA2D
    "Asymmetry parameter"
    g::FTA2D
end
Adapt.@adapt_structure TwoStream

function TwoStream(
    ::Type{FT},
    ncol::Int,
    nlay::Int,
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return TwoStream{FT,DA{FT,2}}(
        DA{FT,2}(zeros(nlay, ncol)),
        DA{FT,2}(zeros(nlay, ncol)),
        DA{FT,2}(zeros(nlay, ncol)),
    )
end

function init_optical_props(
    opc::Symbol,
    ::Type{FT},
    ::Type{DA},
    ncol::Int,
    nlay::Int,
) where {FT<:AbstractFloat,DA}
    if opc == :OneScalar
        return OneScalar(FT, ncol, nlay, DA)
    else
        return TwoStream(FT, ncol, nlay, DA)
    end

end

function compute_col_dry!(
    p_lev::FTA2D,
    t_lay::FTA2D,
    col_dry::FTA2D,
    param_set::AbstractEarthParameterSet,
    vmr_h2o::Union{FTA2D,Nothing} = nothing,
    lat::Union{FTA1D,Nothing} = nothing,
    max_threads::Int = Int(256),
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}
    nlay = size(t_lay, 1)
    ncol = size(t_lay, 2)

    mol_m_dry = FT(molmass_dryair(param_set))
    mol_m_h2o = FT(molmass_water(param_set))
    avogadro = FT(avogad())
    helmert1 = FT(grav(param_set))
    #------Launching computation kernel-----------------------
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlay)

    device = array_device(p_lev)
    comp_stream = Event(device)
    workgroup = (thr_x, thr_y)
    ndrange = (nlay, ncol)

    comp_stream = compute_col_dry_kernel!(device, workgroup)(
        Val(nlay),
        Val(ncol),
        p_lev,
        t_lay,
        col_dry,
        mol_m_dry,
        mol_m_h2o,
        avogadro,
        helmert1,
        vmr_h2o,
        lat,
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    wait(comp_stream)

    return nothing
end

function compute_optical_props!(
    as::AtmosphericState{FT,FTA1D,FTA1DN,FTA2D,VMR,I},
    lkp::AbstractLookUp{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D},
    op::AbstractOpticalProps{FT,FTA2D},
    igpt::Int;
    islw::Bool = true,
    sf::Union{Nothing,AbstractSource{FT}} = nothing,
    max_threads = Int(256),
) where {
    I<:Int,
    FT<:AbstractFloat,
    UI8<:UInt8,
    UI8A1D<:AbstractArray{UI8,1},
    IA1D<:AbstractArray{I,1},
    IA2D<:AbstractArray{I,2},
    IA3D<:AbstractArray{I,3},
    FTA1D<:AbstractArray{FT,1},
    FTA1DN<:Union{AbstractArray{FT,1},Nothing},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    FTA4D<:AbstractArray{FT,4},
    VMR<:AbstractVmr{FT},
}
    nlay = as.nlay
    ncol = as.ncol
    #------Launching computation kernel-----------------------
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlay)

    device = array_device(as.t_lay)
    comp_stream = Event(device)
    workgroup = (thr_x, thr_y)
    ndrange = (nlay, ncol)

    comp_stream = compute_optical_props_kernel!(device, workgroup)(
        Val(nlay),
        Val(ncol),
        as,
        lkp,
        op,
        igpt,
        islw,
        sf,
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    wait(comp_stream)

    return nothing
end

function compute_optical_props!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I},
    op::AbstractOpticalProps{FT,FTA2D};
    islw::Bool = true,
    sf::Union{AbstractSource{FT},Nothing} = nothing,
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}
    ncol = as.ncol
    nlay = as.nlay
    τ = op.τ
    #----Launcing KA Kernel---------------------------
    max_threads = 256
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlay)

    device = array_device(τ)
    comp_stream = Event(device)
    workgroup = (thr_x, thr_y)
    ndrange = (nlay, ncol)
    #----------------------------------
    comp_stream = compute_optical_props_kernel!(device, workgroup)(
        as,
        τ,
        islw,
        sf,
        Val(nlay),
        Val(ncol),
        ndrange = ndrange,
        dependencies = (comp_stream,),
    )
    wait(comp_stream)
    #----------------------------------
    if typeof(op) == TwoStream{FT,FTA2D}
        op.ssa .= FT(0)
        op.g .= FT(0)
    end
    return nothing
end

include("OpticsKernels.jl")

end
