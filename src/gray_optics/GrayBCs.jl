module GrayBCs

using ..Device: array_type
using DocStringExtensions
using Adapt

export AbstractGrayBCs, GrayLwBCs

abstract type AbstractGrayBCs{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA1DN<:Union{Nothing,AbstractArray{FT,1}},
} end

"""
    GrayLwBCs{FT,FTA1D}

Longwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayLwBCs{FT,FTA1D,FTA1DN} <: AbstractGrayBCs{FT,FTA1D,FTA1DN}
    "Surface emissivity `[W/m²]` `(ncol)`"
    sfc_emis::FTA1D
    "incident flux at top of atmosphere `[W/m²]` `(ncol)`"
    inc_flux::FTA1DN
end

function Adapt.adapt_structure(to, x::GrayLwBCs)
    FT = eltype(x.sfc_emis)
    FTA1D = typeof(adapt(to, x.sfc_emis))
    FTA1DN = typeof(adapt(to, x.inc_flux))
    GrayLwBCs{FT,FTA1D,FTA1DN}(adapt(to, x.sfc_emis), adapt(to, x.inc_flux))
end

function GrayLwBCs(::Type{DA}, sfc_emis, inc_flux = nothing) where {DA}
    FT = eltype(sfc_emis)

    @assert length(size(sfc_emis)) == 1
    sfc_emis = DA{FT}(sfc_emis)
    if inc_flux ≠ nothing
        @assert length(size(inc_flux)) == 1
        inc_flux = DA{FT}(inc_flux)
    end

    return GrayLwBCs{FT,DA{FT,1},typeof(inc_flux)}(sfc_emis, inc_flux)
end

"""
    GraySwBCs{FT,FTA1D}
Shortwave boundary conditions
# Fields
$(DocStringExtensions.FIELDS)
"""
struct GraySwBCs{FT,FTA1D,FTA1DN} <: AbstractGrayBCs{FT,FTA1D,FTA1DN}
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

function Adapt.adapt_structure(to, x::GraySwBCs)
    FT = eltype(x.toa_flux)
    FTA1D = typeof(adapt(to, x.toa_flux))
    FTA1DN = typeof(adapt(to, x.inc_flux_diffuse))
    GraySwBCs{FT,FTA1D,FTA1DN}(
        adapt(to, x.toa_flux),
        adapt(to, x.sfc_alb_direct),
        adapt(to, x.sfc_alb_diffuse),
        adapt(to, x.inc_flux_diffuse),
        adapt(to, x.zenith),
    )
end

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

    return new{FT,typeof(toa_flux),typeof(inc_flux_diffuse)}(
        toa_flux,
        sfc_alb_direct,
        sfc_alb_diffuse,
        inc_flux_diffuse,
        zenith,
    )
end

end
