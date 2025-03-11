module Fluxes

using Adapt
using DocStringExtensions

export AbstractFlux, FluxLW, FluxSW, set_flux!, set_net_flux!, add_to_flux!, compute_net_flux!

abstract type AbstractFlux end


"""
    FluxLW{D} <: AbstractFlux

Upward, downward and net longwave fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxLW{D} <: AbstractFlux
    "upward flux `[W/m²]` `(ncol, nlev)`"
    flux_up::D
    "downward flux `[W/m²]` `(ncol, nlev)`"
    flux_dn::D
    "net flux `[W/m²]` `(ncol, nlev)`"
    flux_net::D
end
Adapt.@adapt_structure FluxLW

function FluxLW(ncol::Int, nlay::Int, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}
    flux_up = DA{FT}(undef, ncol, nlay + 1)
    flux_dn = DA{FT}(undef, ncol, nlay + 1)
    flux_net = DA{FT}(undef, ncol, nlay + 1)
    return FluxLW(flux_up, flux_dn, flux_net)
end

"""
    FluxSW{D} <: AbstractFlux

Upward, downward, net and direct downward shortwave fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxSW{D} <: AbstractFlux
    "upward flux `[W/m²]` `(ncol, nlev)`"
    flux_up::D
    "downward flux `[W/m²]` `(ncol, nlev)`"
    flux_dn::D
    "net flux `[W/m²]` `(ncol, nlev)`"
    flux_net::D
    "direct downward flux `[W/m²]` `(ncol, nlev)`"
    flux_dn_dir::D
end
Adapt.@adapt_structure FluxSW

function FluxSW(ncol::Int, nlay::Int, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}
    flux_up = DA{FT}(undef, ncol, nlay + 1)
    flux_dn = DA{FT}(undef, ncol, nlay + 1)
    flux_net = DA{FT}(undef, ncol, nlay + 1)
    flux_dn_dir = DA{FT}(undef, ncol, nlay + 1)
    return FluxSW(flux_up, flux_dn, flux_net, flux_dn_dir)
end

# return number of levels and columns
@inline get_dims(flux::AbstractFlux) = (size(flux.flux_up, 2), size(flux.flux_up, 1))

"""
    set_flux!(flux::AbstractFlux{FT}, val::FT, gcol::Int) where {FT <: AbstractFloat}

Set `flux.flux_up` and `flux.flux_dn` for column `gcol` to `val` 
"""
@inline function set_flux!(flux::AbstractFlux, val::FT, gcol::Int) where {FT}
    (; flux_up, flux_dn) = flux
    nlev, _ = get_dims(flux)
    @inbounds for ilev in 1:nlev
        flux_up[gcol, ilev] = val
        flux_dn[gcol, ilev] = val
    end
    return nothing
end

"""
    set_flux!(flux1::F, flux2::F, gcol::Int) where {F <: AbstractFlux}

Copies `flux2` into `flux1`
"""
@inline function set_flux!(flux1::F, flux2::F, gcol::Int) where {F <: AbstractFlux}
    flux_up1, flux_dn1 = flux1.flux_up, flux1.flux_dn
    flux_up2, flux_dn2 = flux2.flux_up, flux2.flux_dn
    nlev, _ = get_dims(flux1)
    @inbounds begin
        for ilev in 1:nlev
            flux_up1[gcol, ilev] = flux_up2[gcol, ilev]
            flux_dn1[gcol, ilev] = flux_dn2[gcol, ilev]
        end
    end
    return nothing
end

"""
    set_net_flux!(flux::AbstractFlux, val::FT, gcol) where {FT}

Sets net flux (`flux.flux_net`) to `val`
"""
@inline function set_net_flux!(flux::AbstractFlux, val::FT, gcol) where {FT}
    (; flux_net) = flux
    nlev, _ = get_dims(flux)
    @inbounds begin
        for ilev in 1:nlev
            flux_net[gcol, ilev] = val
        end
    end
    return nothing
end

"""
    compute_net_flux!(flux::AbstractFlux, gcol)

Computes net flux for column `gcol`

`flux.flux_net` = `flux.flux_up` - `flux.flux_dn`
"""
@inline function compute_net_flux!(flux::AbstractFlux, gcol)
    (; flux_up, flux_dn, flux_net) = flux
    nlev, _ = get_dims(flux)
    @inbounds begin
        for ilev in 1:nlev
            flux_net[gcol, ilev] = flux_up[gcol, ilev] - flux_dn[gcol, ilev]
        end
    end
    return nothing
end

"""
    add_to_flux!(flux1::F, flux2::F, gcol::Int) where {F <: AbstractFlux}

Add flux2 to flux1 for column `gcol`
"""
@inline function add_to_flux!(flux1::F, flux2::F, gcol::Int) where {F <: AbstractFlux}
    flux_up1, flux_dn1 = flux1.flux_up, flux1.flux_dn
    flux_up2, flux_dn2 = flux2.flux_up, flux2.flux_dn
    nlev, _ = get_dims(flux1)
    @inbounds begin
        for ilev in 1:nlev
            flux_up1[gcol, ilev] += flux_up2[gcol, ilev]
            flux_dn1[gcol, ilev] += flux_dn2[gcol, ilev]
        end
    end
    return nothing
end

end
