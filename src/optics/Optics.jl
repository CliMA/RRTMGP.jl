module Optics

using DocStringExtensions
using CUDA
using ..Device: array_device, CPU, CUDADevice
using Adapt
#---------------------------------------
using ..Vmrs
using ..LookUpTables
using ..AtmosphericStates
using ..Sources
using ..AngularDiscretizations
import ..Parameters as RP
#---------------------------------------

export AbstractOpticalProps, OneScalar, TwoStream, compute_col_dry!, compute_optical_props!

"""
    AbstractOpticalProps{FT,FTA2D}

Optical properties for one scalar and two stream calculations.
"""
abstract type AbstractOpticalProps{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} end

"""
    OneScalar{FT,FTA1D,FTA2D,I,AD} <: AbstractOpticalProps{FT,FTA2D}

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OneScalar{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}, AD <: AngularDiscretization} <:
       AbstractOpticalProps{FT, FTA2D}
    "Optical Depth"
    τ::FTA2D
    "Angular discretization"
    angle_disc::AD
end
Adapt.@adapt_structure OneScalar

function OneScalar(::Type{FT}, ncol::Int, nlay::Int, ::Type{DA}) where {FT <: AbstractFloat, DA}
    I = Int
    τ = DA{FT, 2}(undef, nlay, ncol)
    ad = AngularDiscretization(FT, DA, I(1))

    return OneScalar{eltype(τ), typeof(τ), typeof(ad)}(τ, ad)
end

"""
    TwoStream{FT,FTA2D} <: AbstractOpticalProps{FT,FTA2D}

Two stream approximation for optical properties, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStream{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} <: AbstractOpticalProps{FT, FTA2D}
    "Optical depth"
    τ::FTA2D
    "Single-scattering albedo"
    ssa::FTA2D
    "Asymmetry parameter"
    g::FTA2D
end
Adapt.@adapt_structure TwoStream

function TwoStream(::Type{FT}, ncol::Int, nlay::Int, ::Type{DA}) where {FT <: AbstractFloat, DA}
    τ = DA{FT, 2}(zeros(nlay, ncol))
    ssa = DA{FT, 2}(zeros(nlay, ncol))
    g = DA{FT, 2}(zeros(nlay, ncol))
    return TwoStream{eltype(τ), typeof(τ)}(τ, ssa, g)
end

"""
    compute_col_dry!(
        p_lev,
        col_dry,
        param_set,
        vmr_h2o,
        lat,
        max_threads = Int(256),
    )

This function computes the column amounts of dry or moist air.

"""
function compute_col_dry!(
    p_lev::FTA2D,
    col_dry::FTA2D,
    param_set::RP.ARP,
    vmr_h2o::Union{AbstractArray{FT, 2}, Nothing} = nothing,
    lat::Union{AbstractArray{FT, 1}, Nothing} = nothing,
    max_threads::Int = Int(256),
) where {FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}}
    nlay, ncol = size(col_dry)
    mol_m_dry = FT(RP.molmass_dryair(param_set))
    mol_m_h2o = FT(RP.molmass_water(param_set))
    avogadro = FT(RP.avogad(param_set))
    helmert1 = FT(RP.grav(param_set))
    args = (p_lev, mol_m_dry, mol_m_h2o, avogadro, helmert1, vmr_h2o, lat)
    device = array_device(p_lev)
    if device === CUDADevice() # launching CUDA kernel
        tx = min(nlay * ncol, max_threads)
        bx = cld(nlay * ncol, tx)
        @cuda threads = (tx) blocks = (bx) compute_col_dry_CUDA!(col_dry, args...)
    else # launching Julia native multithreading kernel
        Threads.@threads for icnt in 1:(nlay * ncol)
            glaycol = ((icnt % nlay == 0) ? nlay : (icnt % nlay), cld(icnt, nlay))
            compute_col_dry_kernel!(col_dry, args..., glaycol)
        end
    end
    return nothing
end

function compute_col_dry_CUDA!(col_dry, args...)
    glx = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlay, ncol = size(col_dry)
    if glx ≤ nlay * ncol
        glaycol = ((glx % nlay == 0) ? nlay : (glx % nlay), cld(glx, nlay))
        compute_col_dry_kernel!(col_dry, args..., glaycol)
    end
    return nothing
end
#-----------------------------------------------------------------------------
function compute_optical_props_CUDA!(op, as, args...)
    glx = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlay, ncol = size(op.τ)
    if glx ≤ nlay * ncol
        glaycol = ((glx % nlay == 0) ? nlay : (glx % nlay), cld(glx, nlay))
        compute_optical_props_kernel!(op, as, glaycol, args...)
    end
    return nothing
end
#-----------------------------------------------------------------------------
"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::AtmosphericState{FT},
        sf::AbstractSourceLW{FT},
        gcol::I,
        igpt::I,
        lkp::LookUpLW{I,FT},
        lkp_cld::Union{LookUpCld,Nothing} = nothing,
    ) where {I<:Int,FT<:AbstractFloat}

Computes optical properties for the longwave problem.
"""
function compute_optical_props!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    sf::AbstractSourceLW{FT},
    gcol::I,
    igpt::I,
    lkp::LookUpLW{I, FT},
    lkp_cld::Union{LookUpCld, Nothing} = nothing,
) where {I <: Int, FT <: AbstractFloat}
    nlay = as.nlay
    lkp_args = (lkp_cld === nothing) ? (lkp,) : (lkp, lkp_cld)
    for ilay in 1:nlay
        glaycol = (ilay, gcol)
        compute_optical_props_kernel!(op, as, glaycol, sf, igpt, lkp_args...)
    end
    return nothing
end

"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::AtmosphericState{FT},
        gcol::I,
        igpt::I,
        lkp::LookUpSW{I,FT},
        lkp_cld::Union{LookUpCld,Nothing} = nothing,
    ) where {I<:Int,FT<:AbstractFloat}

Computes optical properties for the shortwave problem.
"""
function compute_optical_props!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    gcol::I,
    igpt::I,
    lkp::LookUpSW{I, FT},
    lkp_cld::Union{LookUpCld, Nothing} = nothing,
) where {I <: Int, FT <: AbstractFloat}
    nlay = as.nlay
    lkp_args = (lkp_cld === nothing) ? (lkp,) : (lkp, lkp_cld)
    for ilay in 1:nlay
        glaycol = (ilay, gcol)
        compute_optical_props_kernel!(op, as, glaycol, igpt, lkp_args...)
    end
    return nothing
end
#-----------------------------------------------------------------------------
"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::GrayAtmosphericState{FT},
        sf::AbstractSourceLW{FT},
        gcol::I,
        igpt::I = 1,
    ) where {FT<:AbstractFloat}

Computes optical properties for the longwave gray radiation problem.
"""
function compute_optical_props!(
    op::AbstractOpticalProps{FT},
    as::GrayAtmosphericState{FT},
    sf::AbstractSourceLW{FT},
    gcol::I,
    igpt::I = 1,
) where {FT <: AbstractFloat, I <: Int}
    nlay = as.nlay
    for ilay in 1:nlay
        glaycol = (ilay, gcol)
        compute_optical_props_kernel!(op, as, glaycol, sf)
    end
    return nothing
end

"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::GrayAtmosphericState{FT},
        gcol::I,
        igpt::I = 1,
    ) where {FT<:AbstractFloat}

Computes optical properties for the shortwave gray radiation problem.
"""
function compute_optical_props!(
    op::AbstractOpticalProps{FT},
    as::GrayAtmosphericState{FT},
    gcol::I,
    igpt::I = 1,
) where {FT <: AbstractFloat, I <: Int}
    nlay = as.nlay
    for ilay in 1:nlay
        glaycol = (ilay, gcol)
        compute_optical_props_kernel!(op, as, glaycol)
    end
    return nothing
end

include("OpticsKernels.jl")

end
