module Fluxes

using Adapt
using DocStringExtensions
using ClimaComms
import ..RRTMGPGridParams

export AbstractFlux,
    FluxLW,
    FluxSW,
    FluxBand,
    set_flux_to_zero!,
    set_band_flux_to_zero!,
    accumulate_band_flux!,
    add_to_flux!,
    apply_metric_scaling!,
    compute_net_flux!

abstract type AbstractFlux{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} end


"""
    FluxLW{FT,FTA2D}

Upward, downward and net longwave fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxLW{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} <:
       AbstractFlux{FT, FTA2D}
    "upward flux `[W/m²]` `(nlev,ncol)`"
    flux_up::FTA2D
    "downward flux `[W/m²]` `(nlev,ncol)`"
    flux_dn::FTA2D
    "net flux `[W/m²]` `(nlev,ncol)`"
    flux_net::FTA2D
end
FluxLW(flux_up, flux_dn, flux_net) =
    FluxLW{eltype(flux_up), typeof(flux_up)}(flux_up, flux_dn, flux_net)
Adapt.@adapt_structure FluxLW

function FluxLW(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    flux_up = DA{FT}(undef, nlay + 1, ncol)
    flux_dn = DA{FT}(undef, nlay + 1, ncol)
    flux_net = DA{FT}(undef, nlay + 1, ncol)
    return FluxLW{FT, typeof(flux_net)}(flux_up, flux_dn, flux_net)
end

"""
    FluxSW{FT,FTA2D}

Upward, downward and net shortwave fluxes at each level.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxSW{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} <:
       AbstractFlux{FT, FTA2D}
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
    FluxSW{eltype(flux_up), typeof(flux_up)}(
        flux_up,
        flux_dn,
        flux_net,
        flux_dn_dir,
    )
Adapt.@adapt_structure FluxSW

function FluxSW(grid_params::RRTMGPGridParams)
    (; nlay, ncol) = grid_params
    FT = eltype(grid_params)
    DA = ClimaComms.array_type(grid_params)
    flux_up = DA{FT}(undef, nlay + 1, ncol)
    flux_dn = DA{FT}(undef, nlay + 1, ncol)
    flux_net = DA{FT}(undef, nlay + 1, ncol)
    flux_dn_dir = DA{FT}(undef, nlay + 1, ncol)
    return FluxSW{FT, typeof(flux_net)}(flux_up, flux_dn, flux_net, flux_dn_dir)
end

"""
    FluxBand{FT, FTA3D}

Optional per-band upward and downward radiative fluxes at each level,
`(nlev, ncol, n_bnd)`. Only allocated when spectrally-resolved fluxes are requested.
Band `b`'s slice `[:, :, b]` has the same `(nlev, ncol)` layout as the broadband
fluxes, and summing over the band dimension recovers the broadband up/down fluxes.

# Fields
- `flux_up`: upward flux per band [W/m²], `(nlev, ncol, n_bnd)`.
- `flux_dn`: downward flux per band [W/m²], `(nlev, ncol, n_bnd)`.
"""
struct FluxBand{FT <: AbstractFloat, FTA3D <: AbstractArray{FT, 3}}
    flux_up::FTA3D
    flux_dn::FTA3D
end
Adapt.@adapt_structure FluxBand

function FluxBand(grid_params::RRTMGPGridParams, n_bnd::Int)
    (; nlay, ncol) = grid_params
    FT = eltype(grid_params)
    DA = ClimaComms.array_type(grid_params)
    flux_up = DA{FT}(undef, nlay + 1, ncol, n_bnd)
    flux_dn = DA{FT}(undef, nlay + 1, ncol, n_bnd)
    return FluxBand{FT, typeof(flux_up)}(flux_up, flux_dn)
end

# Zero a per-band flux buffer before a solve (no-op when spectral fluxes are off).
set_band_flux_to_zero!(::Nothing) = nothing
function set_band_flux_to_zero!(band::FluxBand{FT}) where {FT}
    band.flux_up .= FT(0)
    band.flux_dn .= FT(0)
    return nothing
end

# Add one g-point's up/down flux (column `gcol`, per-g-point scratch `flux_up`/`flux_dn`
# of shape `(nlev, ncol)`) into its band `ibnd`. No-op when spectral fluxes are off, so
# the broadband path pays nothing (the branch is specialized away on `::Nothing`).
@inline accumulate_band_flux!(::Nothing, flux_up, flux_dn, gcol, ibnd, nlev) =
    nothing
@inline function accumulate_band_flux!(
    band::FluxBand,
    flux_up,
    flux_dn,
    gcol,
    ibnd,
    nlev,
)
    bu, bd = band.flux_up, band.flux_dn
    @inbounds for ilev in 1:nlev
        bu[ilev, gcol, ibnd] += flux_up[ilev, gcol]
        bd[ilev, gcol, ibnd] += flux_dn[ilev, gcol]
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
    nlev = size(flux_up, 1)
    @inbounds begin
        for ilev in 1:nlev
            flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
        end
    end
    return nothing
end

"""
    set_flux_to_zero!(flux::FluxLW{FT}, gcol::Int) where {FT<:AbstractFloat}

Set longwave flux for column `gcol` to zero
"""
@inline function set_flux_to_zero!(
    flux::FluxLW{FT},
    gcol::Int,
) where {FT <: AbstractFloat}
    (; flux_up, flux_dn, flux_net) = flux
    nlev = size(flux_up, 1)
    @inbounds for ilev in 1:nlev
        flux_up[ilev, gcol] = FT(0)
        flux_dn[ilev, gcol] = FT(0)
        flux_net[ilev, gcol] = FT(0)
    end
    return nothing
end

"""
    set_flux_to_zero!(flux::FluxSW{FT}, gcol::Int) where {FT<:AbstractFloat}

Set shortwave flux for column `gcol` to zero

"""
@inline function set_flux_to_zero!(
    flux::FluxSW{FT},
    gcol::Int,
) where {FT <: AbstractFloat}
    (; flux_up, flux_dn, flux_net, flux_dn_dir) = flux
    nlev = size(flux_up, 1)
    @inbounds for ilev in 1:nlev
        flux_up[ilev, gcol] = FT(0)
        flux_dn[ilev, gcol] = FT(0)
        flux_net[ilev, gcol] = FT(0)
        flux_dn_dir[ilev, gcol] = FT(0)
    end
    return nothing
end

"""
    add_to_flux!(flux1::FluxLW, flux2::FluxLW, gcol::Int)

add longwave flux2 to longwave flux1 for column `gcol`
flux1 .+= flux2

"""
@inline function add_to_flux!(flux1::FluxLW, flux2::FluxLW, gcol::Int)
    flux_up1, flux_dn1, flux_net1 = flux1.flux_up, flux1.flux_dn, flux1.flux_net
    flux_up2, flux_dn2, flux_net2 = flux2.flux_up, flux2.flux_dn, flux2.flux_net
    nlev = size(flux_up1, 1)
    @inbounds for ilev in 1:nlev
        flux_up1[ilev, gcol] += flux_up2[ilev, gcol]
        flux_dn1[ilev, gcol] += flux_dn2[ilev, gcol]
        flux_net1[ilev, gcol] += flux_net2[ilev, gcol]
    end
    return nothing
end

"""
    apply_metric_scaling!(flux::Union{FluxLW, FluxSW}, metric_scaling)

Apply metric scaling factor for longwave radiative fluxes. This accounts for geometric expansion of 
grid columns with increasing altitude when deep-atmosphere metric terms are used. Requires
`flux_up`, `flux_dn` and `flux_net` to be available in `flux`. The metric scaling factor is 
an optional argument that has a value `nothing` when shallow atmosphere approximations
are used (in which case this function returns nothing).
"""
function apply_metric_scaling!(
    flux::AbstractFlux,
    metric_scaling::AbstractArray,
)
    flux.flux_up .= flux.flux_up .* metric_scaling
    flux.flux_dn .= flux.flux_dn .* metric_scaling
    flux.flux_net .= flux.flux_net .* metric_scaling
    flux isa FluxSW && (flux.flux_dn_dir .= flux.flux_dn_dir .* metric_scaling)
    return nothing
end
function apply_metric_scaling!(
    flux::Union{FluxLW, FluxSW},
    metric_scaling::Nothing,
)
    nothing
end

# Deep-atmosphere scaling of the optional per-band fluxes. `metric_scaling` `(nlev, ncol)`
# broadcasts across the band dimension, mirroring the broadband scaling so that summing
# the bands still recovers the (scaled) broadband fluxes. No-op when either is absent.
apply_metric_scaling!(::Nothing, metric_scaling) = nothing
apply_metric_scaling!(band::FluxBand, metric_scaling::Nothing) = nothing
function apply_metric_scaling!(band::FluxBand, metric_scaling::AbstractArray)
    band.flux_up .= band.flux_up .* metric_scaling
    band.flux_dn .= band.flux_dn .* metric_scaling
    return nothing
end

end
