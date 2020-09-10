module GrayFluxes

using DocStringExtensions
using ..Device: array_type
using ..GrayOptics

export GrayFluxLW, gray_flux_sw, AbstractGrayFlux

abstract type AbstractGrayFlux{FT<:AbstractFloat} end

"""
    GrayFlux{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}

Upward, downward and net fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayFluxLW{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractGrayFlux{FT}
    "upward flux `[W/m²]` `(nlev,ncol)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(nlev,ncol)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(nlev,ncol)`"
    flux_net::FTA2D

    function GrayFluxLW(
        ncol::Int,
        nlay::Int,
        nlev::Int,
        ::Type{FT},
        ::Type{DA},
    ) where {FT<:AbstractFloat,DA}
        return new{FT,DA{FT,2}}(
            DA{FT}(undef, nlev, ncol),
            DA{FT}(undef, nlev, ncol),
            DA{FT}(undef, nlev, ncol),
        )
    end
end

struct GrayFluxSW2Str{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractGrayFlux{FT}
    "upward flux `[W/m²]` `(nlev,ncol)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(nlev,ncol)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(nlev,ncol)`"
    flux_net::FTA2D
    "downward direct flux `[W/m²]` `(nlev,ncol)`"
    flux_dn_dir::FTA2D
end

function GrayFluxSW2Str(
    ncol::Int,
    nlay::Int,
    nlev::Int,
    ::Type{FT},
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return GrayFluxSW2Str{FT,DA{FT,2}}(
        DA{FT}(undef, nlev, ncol),
        DA{FT}(undef, nlev, ncol),
        DA{FT}(undef, nlev, ncol),
        DA{FT}(undef, nlev, ncol),
    )
end

struct GrayFluxSWNoScat{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractGrayFlux{FT}
    "downward direct flux `[W/m²]` `(nlev,ncol)`"
    flux_dn_dir::FTA2D
end

function GrayFluxSWNoScat(
    ncol::Int,
    nlay::Int,
    nlev::Int,
    ::Type{FT},
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return GrayFluxSWNoScat{FT,DA{FT,2}}(DA{FT}(undef, nlev, ncol),)
end

function gray_flux_sw(ncol, nlay, nlev, FT, DA, optical_props)
    if optical_props == GrayOneScalar
        return GrayFluxSWNoScat(ncol, nlay, nlev, FT, DA)
    else
        return GrayFluxSW2Str(ncol, nlay, nlev, FT, DA)
    end
end

end
