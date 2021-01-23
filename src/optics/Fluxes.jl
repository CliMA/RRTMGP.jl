module Fluxes

using StaticArrays
using KernelAbstractions
using CUDA
using Adapt
using DocStringExtensions
using ..Device: array_type, array_device

export AbstractFlux,
    FluxLW,
    FluxSW2Str,
    FluxSWNoScat,
    FluxSW2Str,
    init_flux_sw,
    set_flux_to_zero!,
    add_to_flux!

abstract type AbstractFlux{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} end


"""
    FluxLW{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}

Upward, downward and net fluxes at each level.

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

struct FluxSWNoScat{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractFlux{FT,FTA2D}
    flux_dn_dir::FTA2D
end
FluxSWNoScat(flux_dn_dir) =
    FluxSWNoScat{eltype(flux_dn_dir),typeof(flux_dn_dir)}(flux_dn_dir)
Adapt.@adapt_structure FluxSWNoScat

function FluxSWNoScat(
    ncol::Int,
    nlay::Int,
    ::Type{FT},
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return FluxSWNoScat{FT,DA{FT,2}}(DA{FT}(undef, nlay + 1, ncol),)
end

struct FluxSW2Str{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractFlux{FT,FTA2D}
    flux_up::FTA2D
    flux_dn::FTA2D
    flux_net::FTA2D
    flux_dn_dir::FTA2D
end
FluxSW2Str(flux_up, flux_dn, flux_net, flux_dn_dir) =
    FluxSW2Str{eltype(flux_up),typeof(flux_up)}(
        flux_up,
        flux_dn,
        flux_net,
        flux_dn_dir,
    )
Adapt.@adapt_structure FluxSW2Str

function FluxSW2Str(
    ncol::Int,
    nlay::Int,
    ::Type{FT},
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    return FluxSW2Str{FT,DA{FT,2}}(
        DA{FT}(undef, nlay + 1, ncol),
        DA{FT}(undef, nlay + 1, ncol),
        DA{FT}(undef, nlay + 1, ncol),
        DA{FT}(undef, nlay + 1, ncol),
    )
end

function init_flux_sw(
    ncol::Int,
    nlay::Int,
    ::Type{FT},
    ::Type{DA},
    opc::Symbol,
) where {FT<:AbstractFloat,DA}
    if opc == :OneScalar
        return FluxSWNoScat(ncol, nlay, FT, DA)
    else
        return FluxSW2Str(ncol, nlay, FT, DA)
    end
end

# sets components of flux struct to zero
function set_flux_to_zero!(
    flux::FluxLW{FT,FTA2D},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    flux.flux_up .= FT(0)
    flux.flux_dn .= FT(0)
    flux.flux_net .= FT(0)
    return nothing
end

function set_flux_to_zero!(
    flux::FluxSWNoScat{FT,FTA2D},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    flux.flux_dn_dir .= FT(0)
    return nothing
end

function set_flux_to_zero!(
    flux::FluxSW2Str{FT,FTA2D},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    flux_up .= FT(0)
    flux_dn .= FT(0)
    flux_net .= FT(0)
    flux.flux_dn_dir .= FT(0)
    return nothing
end

# flux1 += flux2
function add_to_flux!(
    flux1::FluxLW{FT,FTA2D},
    flux2::FluxLW{FT,FTA2D},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    flux1.flux_up .+= flux2.flux_up
    flux1.flux_dn .+= flux2.flux_dn
    flux1.flux_net .+= flux2.flux_net
    return nothing
end

function add_to_flux!(
    flux1::FluxSWNoScat{FT,FTA2D},
    flux2::FluxSWNoScat{FT,FTA2D},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    flux1.flux_dn_dir .+= flux2.flux_dn_dir
    return nothing
end

function add_to_flux!(
    flux1::FluxSW2Str{FT,FTA2D},
    flux2::FluxSW2Str{FT,FTA2D},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    flux1.flux_up .+= flux2.flux_up
    flux1.flux_dn .+= flux2.flux_dn
    flux1.flux_net .+= flux2.flux_net
    flux1.flux_dn_dir .+= flux2.flux_dn_dir
    return nothing
end
#----------------------------------------------
end
