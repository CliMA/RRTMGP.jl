module BCs

using ..Device: array_type
using DocStringExtensions
using Adapt

export AbstractBCs, LwBCs, SwBCs

abstract type AbstractBCs{FT<:AbstractFloat} end

"""
    LwBCs{FT,FTA1D,FTA2DN}

Longwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LwBCs{FT,FTA2D,FTA2DN} <: AbstractBCs{FT}
    "Surface emissivity `[W/m²]` `(nbnd, ncol)`"
    sfc_emis::FTA2D
    "incident flux at top of atmosphere `[W/m²]` `(ncol, ngpt)`"
    inc_flux::FTA2DN
end
LwBCs(sfc_emis, inc_flux) =
    LwBCs{eltype(sfc_emis),typeof(sfc_emis),typeof(inc_flux)}(
        sfc_emis,
        inc_flux,
    )
Adapt.@adapt_structure LwBCs

"""
    SwBCs{FT,FTA1D}
Shortwave boundary conditions
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SwBCs{FT,FTA1D,FTA1DN,FTA2D} <: AbstractBCs{FT}
    "zenith angle"
    zenith::FTA1D
    "top of atmosphere flux"
    toa_flux::FTA1D
    "surface albedo for specular (direct) radiation `(nbnd, ncol)`"
    sfc_alb_direct::FTA2D
    "incident diffuse flux at top of domain `[W/m2]` `(ncol, ngpt)`"
    inc_flux_diffuse::FTA1DN
    "surface albedo for diffuse radiation `(nbnd, ncol)`"
    sfc_alb_diffuse::FTA2D
end
SwBCs(zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse) =
    SwBCs{
        eltype(zenith),
        typeof(zenith),
        typeof(inc_flux_diffuse),
        typeof(sfc_alb_direct),
    }(
        zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
Adapt.@adapt_structure SwBCs

function SwBCs(
    ::Type{DA},
    ::Type{FT},
    zenith::FTA1D,
    toa_flux::FTA1D,
    sfc_alb_direct::FTA2D,
    inc_flux_diffuse::Union{FTA1D,Nothing} = nothing,
    sfc_alb_diffuse::Union{FTA2D,Nothing} = nothing,
) where {
    DA,
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}

    zenith = DA{FT,1}(zenith)
    toa_flux = DA{FT,1}(toa_flux)
    sfc_alb_direct = DA{FT,2}(sfc_alb_direct)

    if sfc_alb_diffuse isa AbstractArray
        sfc_alb_diffuse = DA{FT,2}(sfc_alb_diffuse)
    else
        sfc_alb_diffuse = deepcopy(sfc_alb_direct)
    end

    if inc_flux_diffuse isa AbstractArray
        inc_flux_diffuse = DA{FT,1}(inc_flux_diffuse)
    end

    return SwBCs{
        FT,
        typeof(toa_flux),
        typeof(inc_flux_diffuse),
        typeof(sfc_alb_direct),
    }(
        zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
end

end
