module GrayFluxes

using DocStringExtensions
using ..Device: array_type

export GrayFlux

"""
    GrayFlux{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}

Upward, downward and net fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayFlux{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    "upward flux `[W/m²]` `(nlev,ncol)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(nlev,ncol)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(nlev,ncol)`"
    flux_net::FTA2D

    function GrayFlux(
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

end
