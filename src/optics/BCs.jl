module BCs

using DocStringExtensions
using Adapt

export LwBCs, SwBCs

"""
    LwBCs{D, DN}

Longwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LwBCs{D, DN}
    "Surface emissivity `[W/m²]` `(nbnd, ncol)`"
    sfc_emis::D
    "incident flux at top of atmosphere `[W/m²]` `(ncol, ngpt)`"
    inc_flux::DN
end
Adapt.@adapt_structure LwBCs

"""
    SwBCs{FTA1D, FTA1DN, FTA2D}

Shortwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SwBCs{FTA1D, FTA1DN, FTA2D}
    "cosine of zenith angle `(ncol)`"
    cos_zenith::FTA1D
    "top of atmosphere flux `(ncol)`"
    toa_flux::FTA1D
    "surface albedo for specular (direct) radiation `(nbnd, ncol)`"
    sfc_alb_direct::FTA2D
    "incident diffuse flux at top of domain `[W/m2]` `(ncol, ngpt)`"
    inc_flux_diffuse::FTA1DN
    "surface albedo for diffuse radiation `(nbnd, ncol)`"
    sfc_alb_diffuse::FTA2D
end
Adapt.@adapt_structure SwBCs

end
