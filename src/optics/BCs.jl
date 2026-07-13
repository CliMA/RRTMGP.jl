module BCs

using Adapt
import ..RRTMGPGridParams

export LwBCs, SwBCs

"""
    LwBCs{FT, FTA2D, FTA2DN}

Longwave boundary conditions.

# Fields
- `sfc_emis`: Surface emissivity `(nbnd, ncol)`.
- `inc_flux`: Incident flux at the top of the atmosphere [W/m²] `(ncol, ngpt)`.
"""
struct LwBCs{FT, FTA2D, FTA2DN}
    sfc_emis::FTA2D
    inc_flux::FTA2DN
end
LwBCs(sfc_emis, inc_flux) =
    LwBCs{eltype(sfc_emis), typeof(sfc_emis), typeof(inc_flux)}(
        sfc_emis,
        inc_flux,
    )
Adapt.@adapt_structure LwBCs

"""
    SwBCs{FT, FTA1D, FTA1DN, FTA2D}

Shortwave boundary conditions.

# Fields
- `cos_zenith`: Cosine of the solar zenith angle `(ncol)`.
- `toa_flux`: Top-of-atmosphere flux `(ncol)`.
- `sfc_alb_direct`: Surface albedo for specular (direct) radiation `(nbnd, ncol)`.
- `inc_flux_diffuse`: Incident diffuse flux at the top of the domain [W/m²] `(ncol, ngpt)`.
- `sfc_alb_diffuse`: Surface albedo for diffuse radiation `(nbnd, ncol)`.
"""
struct SwBCs{FT, FTA1D, FTA1DN, FTA2D}
    cos_zenith::FTA1D
    toa_flux::FTA1D
    sfc_alb_direct::FTA2D
    inc_flux_diffuse::FTA1DN
    sfc_alb_diffuse::FTA2D
end
SwBCs(cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse) =
    SwBCs{
        eltype(cos_zenith),
        typeof(cos_zenith),
        typeof(inc_flux_diffuse),
        typeof(sfc_alb_direct),
    }(
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
Adapt.@adapt_structure SwBCs

end
