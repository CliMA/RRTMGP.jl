module GrayFluxes

using DocStringExtensions
using ..Device: array_type
using Adapt

export AbstractGrayFlux, GrayFluxLW

abstract type AbstractGrayFlux{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} end

"""
    GrayFluxLW{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}

Upward, downward and net fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayFluxLW{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractGrayFlux{FT,FTA2D}
    "upward flux `[W/m²]` `(nlev,ncol)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(nlev,ncol)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(nlev,ncol)`"
    flux_net::FTA2D
end

function Adapt.adapt_structure(to, x::GrayFluxLW)
    FT = eltype(x.flux_up)
    FTA2D = typeof(adapt(to, x.flux_up))
    GrayFluxLW{FT,FTA2D}(
        adapt(to, x.flux_up),
        adapt(to, x.flux_dn),
        adapt(to, x.flux_net),
    )
end

function GrayFluxLW(
    ncol::Int,
    nlay::Int,
    nlev::Int,
    ::Type{FT},
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return GrayFluxLW{FT,DA{FT,2}}(
        DA{FT}(undef, nlev, ncol),
        DA{FT}(undef, nlev, ncol),
        DA{FT}(undef, nlev, ncol),
    )
end

end
