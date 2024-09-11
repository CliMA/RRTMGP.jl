module Fluxes

using Adapt
using DocStringExtensions

export AbstractFlux, FluxLW, FluxSW, set_flux_to_zero!, add_to_flux!

abstract type AbstractFlux{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} end


"""
    FluxLW{FT,FTA2D}

Upward, downward and net longwave fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxLW{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} <: AbstractFlux{FT, FTA2D}
    "upward flux `[W/m²]` `(ncol, nlev)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(ncol, nlev)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(ncol, nlev)`"
    flux_net::FTA2D
end
FluxLW(flux_up, flux_dn, flux_net) = FluxLW{eltype(flux_up), typeof(flux_up)}(flux_up, flux_dn, flux_net)
Adapt.@adapt_structure FluxLW

function FluxLW(ncol::Int, nlay::Int, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}
    flux_up = DA{FT}(undef, ncol, nlay + 1)
    flux_dn = DA{FT}(undef, ncol, nlay + 1)
    flux_net = DA{FT}(undef, ncol, nlay + 1)
    return FluxLW{FT, typeof(flux_net)}(flux_up, flux_dn, flux_net)
end

"""
    FluxSW{FT,FTA2D}

Upward, downward and net shortwave fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxSW{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} <: AbstractFlux{FT, FTA2D}
    "upward flux `[W/m²]` `(ncol, nlev)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(ncol, nlev)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(ncol, nlev)`"
    flux_net::FTA2D
    "direct downward flux `[W/m²]` `(ncol, nlev)`"
    flux_dn_dir::FTA2D
end
FluxSW(flux_up, flux_dn, flux_net, flux_dn_dir) =
    FluxSW{eltype(flux_up), typeof(flux_up)}(flux_up, flux_dn, flux_net, flux_dn_dir)
Adapt.@adapt_structure FluxSW

function FluxSW(ncol::Int, nlay::Int, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}
    flux_up = DA{FT}(undef, ncol, nlay + 1)
    flux_dn = DA{FT}(undef, ncol, nlay + 1)
    flux_net = DA{FT}(undef, ncol, nlay + 1)
    flux_dn_dir = DA{FT}(undef, ncol, nlay + 1)
    return FluxSW{FT, typeof(flux_net)}(flux_up, flux_dn, flux_net, flux_dn_dir)
end

"""
    set_flux_to_zero!(flux::FluxLW{FT}) where {FT<:AbstractFloat}

Set longwave flux to zero

"""
function set_flux_to_zero!(flux::FluxLW{FT}) where {FT <: AbstractFloat}
    flux.flux_up .= FT(0)
    flux.flux_dn .= FT(0)
    flux.flux_net .= FT(0)
    return nothing
end

"""
    set_flux_to_zero!(flux::FluxLW{FT}, gcol::Int) where {FT<:AbstractFloat}

Set longwave flux for column `gcol` to zero
"""
function set_flux_to_zero!(flux::FluxLW{FT}, gcol::Int) where {FT <: AbstractFloat}
    (; flux_up, flux_dn, flux_net) = flux
    nlev = size(flux_up, 2)
    @inbounds for ilev in 1:nlev
        flux_up[gcol, ilev] = FT(0)
        flux_dn[gcol, ilev] = FT(0)
        flux_net[gcol, ilev] = FT(0)
    end
    return nothing
end

"""
    set_flux_to_zero!(flux::FluxSW{FT}) where {FT<:AbstractFloat}

Set shortwave flux to zero

"""
function set_flux_to_zero!(flux::FluxSW{FT}) where {FT <: AbstractFloat}
    flux.flux_up .= FT(0)
    flux.flux_dn .= FT(0)
    flux.flux_net .= FT(0)
    flux.flux_dn_dir .= FT(0)
    return nothing
end
"""
    set_flux_to_zero!(flux::FluxSW{FT}, gcol::Int) where {FT<:AbstractFloat}

Set shortwave flux for column `gcol` to zero

"""
function set_flux_to_zero!(flux::FluxSW{FT}, gcol::Int) where {FT <: AbstractFloat}
    (; flux_up, flux_dn, flux_net, flux_dn_dir) = flux
    nlev = size(flux_up, 2)
    @inbounds for ilev in 1:nlev
        flux_up[gcol, ilev] = FT(0)
        flux_dn[gcol, ilev] = FT(0)
        flux_net[gcol, ilev] = FT(0)
        flux_dn_dir[gcol, ilev] = FT(0)
    end
    return nothing
end

"""
    add_to_flux!(flux1::FluxLW, flux2::FluxLW)

add longwave flux2 to longwave flux1 
flux1 .+= flux2

"""
function add_to_flux!(flux1::FluxLW, flux2::FluxLW)
    flux1.flux_up .+= flux2.flux_up
    flux1.flux_dn .+= flux2.flux_dn
    flux1.flux_net .+= flux2.flux_net
    return nothing
end
"""
    add_to_flux!(flux1::FluxLW, flux2::FluxLW, gcol::Int)

add longwave flux2 to longwave flux1 for column `gcol`
flux1 .+= flux2

"""
function add_to_flux!(flux1::FluxLW, flux2::FluxLW, gcol::Int)
    flux_up1, flux_dn1, flux_net1 = flux1.flux_up, flux1.flux_dn, flux1.flux_net
    flux_up2, flux_dn2, flux_net2 = flux2.flux_up, flux2.flux_dn, flux2.flux_net
    nlev = size(flux_up1, 2)
    @inbounds for ilev in 1:nlev
        flux_up1[gcol, ilev] += flux_up2[gcol, ilev]
        flux_dn1[gcol, ilev] += flux_dn2[gcol, ilev]
        flux_net1[gcol, ilev] += flux_net2[gcol, ilev]
    end
    return nothing
end

"""
    add_to_flux!(flux1::FluxSW, flux2::FluxSW)

add shortwave flux2 to shortwave flux1 
flux1 .+= flux2

"""
function add_to_flux!(flux1::FluxSW, flux2::FluxSW)
    flux1.flux_up .+= flux2.flux_up
    flux1.flux_dn .+= flux2.flux_dn
    flux1.flux_net .+= flux2.flux_net
    flux1.flux_dn_dir .+= flux2.flux_dn_dir
    return nothing
end
"""
    add_to_flux!(flux1::FluxSW, flux2::FluxSW, gcol::Int)

add shortwave flux2 to shortwave flux1 for column `gcol`
flux1 .+= flux2

"""
function add_to_flux!(flux1::FluxSW, flux2::FluxSW, gcol)
    flux_up1, flux_dn1, flux_net1, flux_dn_dir1 = flux1.flux_up, flux1.flux_dn, flux1.flux_net, flux1.flux_dn_dir
    flux_up2, flux_dn2, flux_net2, flux_dn_dir2 = flux2.flux_up, flux2.flux_dn, flux2.flux_net, flux2.flux_dn_dir
    nlev = size(flux_up1, 2)
    @inbounds for ilev in 1:nlev
        flux_up1[gcol, ilev] += flux_up2[gcol, ilev]
        flux_dn1[gcol, ilev] += flux_dn2[gcol, ilev]
        flux_net1[gcol, ilev] += flux_net2[gcol, ilev]
        flux_dn_dir1[gcol, ilev] += flux_dn_dir2[gcol, ilev]
    end
    return nothing
end

end
