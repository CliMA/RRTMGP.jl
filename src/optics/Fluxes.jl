module Fluxes

using Adapt
using ClimaComms
import ..RRTMGPGridParams

export AbstractFlux,
    FluxLW,
    FluxSW,
    FluxBand,
    set_flux_to_zero!,
    set_band_flux_to_zero!,
    accumulate_band_flux!,
    apply_metric_scaling!,
    compute_net_flux!

"""
    AbstractFlux

Abstract type for broadband radiative-flux containers: [`FluxLW`](@ref) and
[`FluxSW`](@ref).
"""
abstract type AbstractFlux{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} end


"""
    FluxLW{FT, FTA2D}

Upward, downward and net longwave fluxes at each level.

# Fields
- `flux_up`: Upward flux [W/m²] `(nlev, ncol)`.
- `flux_dn`: Downward flux [W/m²] `(nlev, ncol)`.
- `flux_net`: Net flux [W/m²] `(nlev, ncol)`.
"""
struct FluxLW{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} <:
       AbstractFlux{FT, FTA2D}
    flux_up::FTA2D
    flux_dn::FTA2D
    flux_net::FTA2D
end
FluxLW(flux_up, flux_dn, flux_net) =
    FluxLW{eltype(flux_up), typeof(flux_up)}(flux_up, flux_dn, flux_net)
Adapt.@adapt_structure FluxLW

function FluxLW(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    flux_up = DA{FT}(undef, ncol, nlay + 1)
    flux_dn = DA{FT}(undef, ncol, nlay + 1)
    flux_net = DA{FT}(undef, ncol, nlay + 1)
    return FluxLW{FT, typeof(flux_net)}(flux_up, flux_dn, flux_net)
end

"""
    FluxSW{FT, FTA2D}

Upward, downward and net shortwave fluxes at each level.

# Fields
- `flux_up`: Upward flux [W/m²] `(ncol, nlev)`.
- `flux_dn`: Downward flux [W/m²] `(ncol, nlev)`.
- `flux_net`: Net flux [W/m²] `(ncol, nlev)`.
- `flux_dn_dir`: Direct downward flux [W/m²] `(ncol, nlev)`.
"""
struct FluxSW{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} <:
       AbstractFlux{FT, FTA2D}
    flux_up::FTA2D
    flux_dn::FTA2D
    flux_net::FTA2D
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
    flux_up = DA{FT}(undef, ncol, nlay + 1)
    flux_dn = DA{FT}(undef, ncol, nlay + 1)
    flux_net = DA{FT}(undef, ncol, nlay + 1)
    flux_dn_dir = DA{FT}(undef, ncol, nlay + 1)
    return FluxSW{FT, typeof(flux_net)}(flux_up, flux_dn, flux_net, flux_dn_dir)
end

"""
    FluxBand{FT, FTA3D}

Optional per-band upward, downward, and net radiative fluxes at each level,
`(ncol, nlev, n_bnd)`. Only allocated when spectrally-resolved fluxes are requested.
Band `b`'s slice `[:, :, b]` has the same `(ncol, nlev)` layout as the broadband
fluxes, and summing over the band dimension recovers the broadband fluxes.

# Fields
- `flux_up`: upward flux per band [W/m²], `(ncol, nlev, n_bnd)`.
- `flux_dn`: downward flux per band [W/m²], `(ncol, nlev, n_bnd)`.
- `flux_net`: net flux per band (`flux_up - flux_dn`) [W/m²], `(ncol, nlev, n_bnd)`.
"""
struct FluxBand{FT <: AbstractFloat, FTA3D <: AbstractArray{FT, 3}}
    flux_up::FTA3D
    flux_dn::FTA3D
    flux_net::FTA3D
end
Adapt.@adapt_structure FluxBand

function FluxBand(grid_params::RRTMGPGridParams, n_bnd::Int)
    (; nlay, ncol) = grid_params
    FT = eltype(grid_params)
    DA = ClimaComms.array_type(grid_params)
    flux_up = DA{FT}(undef, ncol, nlay + 1, n_bnd)
    flux_dn = DA{FT}(undef, ncol, nlay + 1, n_bnd)
    flux_net = DA{FT}(undef, ncol, nlay + 1, n_bnd)
    return FluxBand{FT, typeof(flux_up)}(flux_up, flux_dn, flux_net)
end

# Zero a per-band flux buffer before a solve (no-op when spectral fluxes are off).
set_band_flux_to_zero!(::Nothing) = nothing
function set_band_flux_to_zero!(band::FluxBand{FT}) where {FT}
    band.flux_up .= FT(0)
    band.flux_dn .= FT(0)
    band.flux_net .= FT(0)
    return nothing
end

# Add one g-point's up/down flux (column `gcol`, per-g-point scratch `flux_up`/`flux_dn`
# of shape `(ncol, nlev)`) into its band `ibnd`. No-op when spectral fluxes are off, so
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
        bu[gcol, ilev, ibnd] += flux_up[gcol, ilev]
        bd[gcol, ilev, ibnd] += flux_dn[gcol, ilev]
    end
    return nothing
end

"""
    compute_net_flux!(flux::AbstractFlux, gcol, nlev)
    compute_net_flux!(flux::AbstractFlux, gcol)

Compute the net flux for column `gcol` across `nlev` levels:

`flux.flux_net` = `flux.flux_up` - `flux.flux_dn`
"""
@inline function compute_net_flux!(flux::AbstractFlux, gcol, nlev)
    (; flux_up, flux_dn, flux_net) = flux
    @inbounds begin
        for ilev in 1:nlev
            flux_net[gcol, ilev] = flux_up[gcol, ilev] - flux_dn[gcol, ilev]
        end
    end
    return nothing
end

@inline compute_net_flux!(flux::AbstractFlux, gcol) =
    compute_net_flux!(flux, gcol, size(flux.flux_up, 2))

"""
    set_flux_to_zero!(flux::FluxLW{FT}, gcol::Int, nlev::Int) where {FT<:AbstractFloat}
    set_flux_to_zero!(flux::FluxLW{FT}, gcol::Int) where {FT<:AbstractFloat}

Set longwave flux for column `gcol` to zero across `nlev` levels.
"""
@inline function set_flux_to_zero!(
    flux::FluxLW{FT},
    gcol::Int,
    nlev::Int,
) where {FT <: AbstractFloat}
    (; flux_up, flux_dn, flux_net) = flux
    @inbounds for ilev in 1:nlev
        flux_up[gcol, ilev] = FT(0)
        flux_dn[gcol, ilev] = FT(0)
        flux_net[gcol, ilev] = FT(0)
    end
    return nothing
end

@inline set_flux_to_zero!(flux::FluxLW, gcol::Int) =
    set_flux_to_zero!(flux, gcol, size(flux.flux_up, 2))

"""
    set_flux_to_zero!(flux::FluxSW{FT}, gcol::Int, nlev::Int) where {FT<:AbstractFloat}
    set_flux_to_zero!(flux::FluxSW{FT}, gcol::Int) where {FT<:AbstractFloat}

Set shortwave flux for column `gcol` to zero across `nlev` levels.
"""
@inline function set_flux_to_zero!(
    flux::FluxSW{FT},
    gcol::Int,
    nlev::Int,
) where {FT <: AbstractFloat}
    (; flux_up, flux_dn, flux_net, flux_dn_dir) = flux
    @inbounds for ilev in 1:nlev
        flux_up[gcol, ilev] = FT(0)
        flux_dn[gcol, ilev] = FT(0)
        flux_net[gcol, ilev] = FT(0)
        flux_dn_dir[gcol, ilev] = FT(0)
    end
    return nothing
end

@inline set_flux_to_zero!(flux::FluxSW, gcol::Int) =
    set_flux_to_zero!(flux, gcol, size(flux.flux_up, 2))

"""
    apply_metric_scaling!(flux::Union{FluxLW, FluxSW}, metric_scaling)

Apply the metric scaling factor to radiative fluxes (longwave or shortwave; for a `FluxSW` the
direct beam `flux_dn_dir` is scaled too). This accounts for geometric expansion of 
grid columns with increasing altitude when deep-atmosphere metric terms are used. Requires
`flux_up`, `flux_dn` and `flux_net` to be available in `flux`. The metric scaling factor is 
an optional argument that has a value `nothing` when shallow atmosphere approximations
are used (in which case this function returns nothing).
"""
function apply_metric_scaling!(
    flux::AbstractFlux,
    metric_scaling::AbstractArray,
)
    sc = PermutedDimsArray(metric_scaling, (2, 1))
    flux.flux_up .= flux.flux_up .* sc
    flux.flux_dn .= flux.flux_dn .* sc
    flux.flux_net .= flux.flux_net .* sc
    flux isa FluxSW && (flux.flux_dn_dir .= flux.flux_dn_dir .* sc)
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
    sc = PermutedDimsArray(metric_scaling, (2, 1))
    band.flux_up .= band.flux_up .* sc
    band.flux_dn .= band.flux_dn .* sc
    return nothing
end

end
