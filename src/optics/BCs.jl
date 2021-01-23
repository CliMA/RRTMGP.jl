module BCs

using ..Device: array_type
using DocStringExtensions
using Adapt

export AbstractBCs, LwBCs, SwBCs

abstract type AbstractBCs{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA1DN<:Union{Nothing,AbstractArray{FT,1}},
} end

"""
    LwBCs{FT,FTA1D}

Longwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LwBCs{FT,FTA1D,FTA1DN} <: AbstractBCs{FT,FTA1D,FTA1DN}
    "Surface emissivity `[W/m²]` `(ncol)`"
    sfc_emis::FTA1D
    "incident flux at top of atmosphere `[W/m²]` `(ncol)`"
    inc_flux::FTA1DN
end
LwBCs(sfc_emis, inc_flux) = LwBCs{eltype(sfc_emis), typeof(sfc_emis), typeof(inc_flux)}(sfc_emis, inc_flux)
Adapt.@adapt_structure LwBCs

"""
    SwBCs{FT,FTA1D}
Shortwave boundary conditions
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SwBCs{FT,FTA1D,FTA1DN} <: AbstractBCs{FT,FTA1D,FTA1DN}
    "top of atmosphere flux"
    toa_flux::FTA1D
    "surface albedo for specular (direct) radiation"
    sfc_alb_direct::FTA1D
    "surface albedo for diffuse radiation"
    sfc_alb_diffuse::FTA1D
    "incident diffuse flux at top of domain `[W/m2]` `(ncol, ngpt)`"
    inc_flux_diffuse::FTA1DN
    "zenith angle `radians`"
    zenith::FTA1D
end
Adapt.@adapt_structure SwBCs

function SwBCs(
    ::Type{DA},
    ::Type{FT},
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

    return SwBCs{FT,typeof(toa_flux),typeof(inc_flux_diffuse)}(
        toa_flux,
        sfc_alb_direct,
        sfc_alb_diffuse,
        inc_flux_diffuse,
        zenith,
    )
end

end
