module GrayFluxes

using ..Device: array_type

export GrayFlux

struct GrayFlux{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    flux_up::FTA2D
    flux_dn::FTA2D
    flux_net::FTA2D
    hr_lay::FTA2D

    function GrayFlux(
        ncol::Int,
        nlay::Int,
        nlev::Int,
        ::Type{FT},
        ::Type{DA},
    ) where {FT<:AbstractFloat,DA}
        return new{FT,DA{FT,2}}(
            DA{FT}(undef, ncol, nlev),
            DA{FT}(undef, ncol, nlev),
            DA{FT}(undef, ncol, nlev),
            DA{FT}(undef, ncol, nlay),
        )
    end
end

end
