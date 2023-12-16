module BCs

using DocStringExtensions
using Adapt

export LwBCs, SwBCs, SwBCsc

"""
    LwBCs{FT,FTA1D,FTA2DN}

Longwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LwBCs{FT, FTA2D, FTA2DN}
    "Surface emissivity `[W/m²]` `(nbnd, ncol)`"
    sfc_emis::FTA2D
    "incident flux at top of atmosphere `[W/m²]` `(ncol, ngpt)`"
    inc_flux::FTA2DN
end
LwBCs(sfc_emis, inc_flux) = LwBCs{eltype(sfc_emis), typeof(sfc_emis), typeof(inc_flux)}(sfc_emis, inc_flux)
Adapt.@adapt_structure LwBCs

"""
    SwBCs{FT,FTA1D}
Shortwave boundary conditions
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SwBCs{FT, FTA1D, FTA1DN, FTA2D}
    "cosine of zenith angle `(ncol)`"
    cos_zenith::FTA1D
    "top of atmosphere flux `(ncol)`"
    toa_flux::FTA1D
    "surface albedo for specular (direct) radiation `(nbnd, ncol)`"
    sfc_alb_direct::FTA2D
    "incident diffuse flux at top of domain `[W/m2]` `(ngpt, ncol)`"
    inc_flux_diffuse::FTA1DN
    "surface albedo for diffuse radiation `(nbnd, ncol)`"
    sfc_alb_diffuse::FTA2D
end
SwBCs(cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse) =
    SwBCs{eltype(cos_zenith), typeof(cos_zenith), typeof(inc_flux_diffuse), typeof(sfc_alb_direct)}(
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
Adapt.@adapt_structure SwBCs

struct SwBCsc{FTVA1D, FTA1D, FTA1DN}
    cos_zenith::FTVA1D
    toa_flux::FTVA1D
    sfc_alb_direct::FTA1D
    inc_flux_diffuse::FTA1DN
    sfc_alb_diffuse::FTA1D
end
Adapt.@adapt_structure SwBCsc

function SwBCsc(swbcs::SwBCs, gcol::Int)
    inc_flux_diffuse = swbcs.inc_flux_diffuse isa Nothing ? nothing : view(swbcs.inc_flux_diffuse, :, gcol)
    return SwBCsc(
        view(swbcs.cos_zenith, gcol:gcol),
        view(swbcs.toa_flux, gcol:gcol),
        view(swbcs.sfc_alb_direct, :, gcol),
        inc_flux_diffuse,
        view(swbcs.sfc_alb_diffuse, :, gcol),
    )
end

end
