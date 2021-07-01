module Fluxes

using Adapt
using DocStringExtensions

export AbstractFlux,
    FluxLW, FluxSW, init_flux_sw, set_flux_to_zero!, add_to_flux!

abstract type AbstractFlux{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} end


"""
    FluxLW{FT,FTA2D}

Upward, downward and net longwave fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxLW{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractFlux{FT,FTA2D}
    "upward flux `[W/m²]` `(nlev,ncol)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(nlev,ncol)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(nlev,ncol)`"
    flux_net::FTA2D
end
FluxLW(flux_up, flux_dn, flux_net) =
    FluxLW{eltype(flux_up),typeof(flux_up)}(flux_up, flux_dn, flux_net)
Adapt.@adapt_structure FluxLW

function FluxLW(
    ncol::Int,
    nlay::Int,
    ::Type{FT},
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return FluxLW{FT,DA{FT,2}}(
        DA{FT}(undef, nlay + 1, ncol),
        DA{FT}(undef, nlay + 1, ncol),
        DA{FT}(undef, nlay + 1, ncol),
    )
end

"""
    FluxSW{FT,FTA2D}

Upward, downward and net shortwave fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxSW{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractFlux{FT,FTA2D}
    "upward flux `[W/m²]` `(nlev,ncol)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(nlev,ncol)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(nlev,ncol)`"
    flux_net::FTA2D
    "direct downward flux `[W/m²]` `(nlev,ncol)`"
    flux_dn_dir::FTA2D
end
FluxSW(flux_up, flux_dn, flux_net, flux_dn_dir) =
    FluxSW{eltype(flux_up),typeof(flux_up)}(
        flux_up,
        flux_dn,
        flux_net,
        flux_dn_dir,
    )
Adapt.@adapt_structure FluxSW

function FluxSW(
    ncol::Int,
    nlay::Int,
    ::Type{FT},
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return FluxSW{FT,DA{FT,2}}(
        DA{FT}(undef, nlay + 1, ncol),
        DA{FT}(undef, nlay + 1, ncol),
        DA{FT}(undef, nlay + 1, ncol),
        DA{FT}(undef, nlay + 1, ncol),
    )
end

# sets components of flux struct to zero
function set_flux_to_zero!(flux::FluxLW{FT}) where {FT<:AbstractFloat}
    flux.flux_up .= FT(0)
    flux.flux_dn .= FT(0)
    flux.flux_net .= FT(0)
    return nothing
end

function set_flux_to_zero!(flux::FluxSW{FT}) where {FT<:AbstractFloat}
    flux.flux_up .= FT(0)
    flux.flux_dn .= FT(0)
    flux.flux_net .= FT(0)
    flux.flux_dn_dir .= FT(0)
    return nothing
end

# flux1 .+= flux2
function add_to_flux!(flux1::FluxLW, flux2::FluxLW)
    flux1.flux_up .+= flux2.flux_up
    flux1.flux_dn .+= flux2.flux_dn
    flux1.flux_net .+= flux2.flux_net
    return nothing
end

function add_to_flux!(flux1::FluxSW, flux2::FluxSW)
    flux1.flux_up .+= flux2.flux_up
    flux1.flux_dn .+= flux2.flux_dn
    flux1.flux_net .+= flux2.flux_net
    flux1.flux_dn_dir .+= flux2.flux_dn_dir
    return nothing
end
#----------------------------------------------
end
