module GrayBCs

using ..Device: array_type
using DocStringExtensions

export AbstractGrayBCs, GrayLwBCs, GraySwBCs

abstract type AbstractGrayBCs{FT<:AbstractFloat,FTA1D<:AbstractArray{FT,1}} end

"""
    GrayLwBCs{FT,FTA1D}

Longwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayLwBCs{FT,FTA1D} <: AbstractGrayBCs{FT,FTA1D}
    "Surface emissivity `[W/m²]` `(ncol)`"
    sfc_emis::FTA1D
    "incident flux at top of atmosphere `[W/m²]` `(ncol)`"
    inc_flux::Union{FTA1D,Nothing}

    function GrayLwBCs(::Type{DA}, sfc_emis, inc_flux = nothing) where {DA}
        FT = eltype(sfc_emis)

        @assert length(size(sfc_emis)) == 1
        sfc_emis = DA{FT}(sfc_emis)
        if inc_flux ≠ nothing
            @assert length(size(inc_flux)) == 1
            inc_flux = DA{FT}(inc_flux)
        end

        return new{FT,DA{FT,1}}(sfc_emis, inc_flux)
    end
end

"""
    GraySwBCs{FT,FTA1D}

Shortwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GraySwBCs{FT,FTA1D} <: AbstractGrayBCs{FT,FTA1D}
    "top of atmosphere flux"
    toa_flux::FTA1D
    "surface albedo for specular (direct) radiation"
    sfc_alb_direct::FTA1D
    "surface albedo for diffuse radiation"
    sfc_alb_diffuse::FTA1D
    "incident diffuse flux at top of domain `[W/m2]` `(ncol, ngpt)`"
    inc_flux_diffuse::Union{FTA1D,Nothing}
    "zenith angle `radians`"
    zenith::FTA1D

    function GraySwBCs(
        ::Type{DA},
        toa_flux::FTA1D,
        sfc_alb_direct::FTA1D,
        sfc_alb_diffuse::FTA1D,
        zenith::FTA1D,
        inc_flux_diffuse::Union{FTA1D,Nothing} = nothing,
    ) where {DA,FT<:AbstractFloat,FTA1D<:AbstractArray{FT,1}}

        @assert length(size(toa_flux)) == 1
        @assert length(size(sfc_alb_direct)) == 1
        @assert length(size(sfc_alb_diffuse)) == 1
        @assert length(size(zenith)) == 1

        toa_flux = DA{FT}(toa_flux)
        sfc_alb_direct = DA{FT}(sfc_alb_direct)
        sfc_alb_diffuse = DA{FT}(sfc_alb_diffuse)
        zenith = DA{FT}(zenith)

        if inc_flux_diffuse ≠ nothing
            @assert length(size(inc_flux_diffuse)) == 1
            inc_flux_diffuse = DA{FT}(inc_flux_diffuse)
        end

        return new{FT,DA{FT,1}}(
            toa_flux,
            sfc_alb_direct,
            sfc_alb_diffuse,
            inc_flux_diffuse,
            zenith,
        )
    end
end

end
