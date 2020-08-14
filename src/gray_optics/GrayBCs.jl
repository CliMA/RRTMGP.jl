module GrayBCs

using ..Device: array_type
using DocStringExtensions

export AbstractGrayBCs, GrayLwBCs

abstract type AbstractGrayBCs{FT<:AbstractFloat,FTA1D<:AbstractArray{FT,1}} end

"""
    GrayLwBCs{FT}

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

#=
"""
    GraySwBCs{FT}

Longwave boundary conditions

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
    inc_flux_diffuse::FTA1D
end

=#
end
