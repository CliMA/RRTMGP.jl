module GrayOptics

using DocStringExtensions
using ..Device: array_type
using Adapt

export AbstractGrayOpticalProps, GrayOneScalar, GrayTwoStream

abstract type AbstractGrayOpticalProps{
    FT<:AbstractFloat,
    FTA2D<:AbstractArray{FT,2},
} end

"""
    GrayOneScalar{FT    <: AbstractFloat,
                  FTA2D <: AbstractArray{FT,2}} <: AbstractGrayOpticalProps{FT,FTA2D

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayOneScalar{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractGrayOpticalProps{FT,FTA2D}
    "Optical Depth"
    τ::FTA2D
end

function Adapt.adapt_structure(to, x::GrayOneScalar)
    FT = eltype(x.τ)
    FTA2D = typeof(adapt(to, x.τ))
    GrayOneScalar{FT,FTA2D}(adapt(to, x.τ))
end

function GrayOneScalar(
    ::Type{FT},
    ncol::Int,
    nlay::Int,
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return GrayOneScalar{FT,DA{FT,2}}(DA{FT,2}(undef, nlay, ncol))
end

"""
    GrayTwoStream{FT<:AbstractFloat,
                     FTA2D<:AbstractArray{FT,2}} <: AbstractGrayOpticalProps{FT,FTA2D}

Two stream approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayTwoStream{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractGrayOpticalProps{FT,FTA2D}
    "Optical depth"
    τ::FTA2D
    "Single-scattering albedo"
    ssa::FTA2D
    "Asymmetry parameter"
    g::FTA2D
end

function Adapt.adapt_structure(to, x::GrayTwoStream)
    FT = eltype(x.τ)
    FTA2D = typeof(adapt(to, x.τ))
    GrayTwoStream{FT,FTA2D}(adapt(to, x.τ), adapt(to, x.ssa), adapt(to, x.g))
end

function GrayTwoStream(
    ::Type{FT},
    ncol::Int,
    nlay::Int,
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return GrayTwoStream{FT,DA{FT,2}}(
        DA{FT,2}(undef, nlay, ncol),
        DA{FT,2}(undef, nlay, ncol),
        DA{FT,2}(undef, nlay, ncol),
    )
end

end
