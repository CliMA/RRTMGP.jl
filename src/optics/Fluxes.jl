module Fluxes

using Adapt
using ClimaComms
import ..RRTMGPGridParams

export AbstractFlux,
    FluxLW,
    FluxSW,
    FluxBand,
    FluxPresentation,
    update_presentation!,
    transpose_into!,
    transpose_sum_into!,
    set_flux_to_zero!,
    set_band_flux_to_zero!,
    accumulate_band_flux!,
    apply_metric_scaling!,
    compute_net_flux!

"""
    lazy_transpose(x)

Lazily present the `(ncol, nlev[, n_bnd])` internal buffer `x` as
`(nlev, ncol[, n_bnd])` without copying. Equivalent to
`PermutedDimsArray(x, perm)`, but with the permutation in the type parameters:
before Julia 1.12, the `PermutedDimsArray(x, perm)` constructor does not
constant-fold `perm`, which would make every flux getter type-unstable (and
allocating) on Julia 1.9–1.11.
"""
lazy_transpose(x::AbstractArray{T, 2}) where {T} =
    PermutedDimsArray{T, 2, (2, 1), (2, 1), typeof(x)}(x)
lazy_transpose(x::AbstractArray{T, 3}) where {T} =
    PermutedDimsArray{T, 3, (2, 1, 3), (2, 1, 3), typeof(x)}(x)

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
- `flux_up`: Upward flux [W/m²] `(ncol, nlev)`.
- `flux_dn`: Downward flux [W/m²] `(ncol, nlev)`.
- `flux_net`: Net flux [W/m²] `(ncol, nlev)`.
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
`(nlev, ncol, n_bnd)`. Only allocated when spectrally-resolved fluxes are requested.
Summing over the band dimension recovers the broadband fluxes.

Unlike the broadband compute buffers, the band buffers keep the host-facing
vertical-first layout: they are an opt-in diagnostic that the `spectral_*`
getters expose as plain views, and keeping them `(nlev, ncol, n_bnd)` avoids
doubling their (large) memory with separate presentation copies. The
accumulation writes are uncoalesced on the GPU, a cost paid only when
per-band fluxes are requested.

# Fields
- `flux_up`: upward flux per band [W/m²], `(nlev, ncol, n_bnd)`.
- `flux_dn`: downward flux per band [W/m²], `(nlev, ncol, n_bnd)`.
- `flux_net`: net flux per band (`flux_up - flux_dn`) [W/m²], `(nlev, ncol, n_bnd)`.
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
    flux_up = DA{FT}(undef, nlay + 1, ncol, n_bnd)
    flux_dn = DA{FT}(undef, nlay + 1, ncol, n_bnd)
    flux_net = DA{FT}(undef, nlay + 1, ncol, n_bnd)
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
        bu[ilev, gcol, ibnd] += flux_up[gcol, ilev]
        bd[ilev, gcol, ibnd] += flux_dn[gcol, ilev]
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
    _scale_by_transpose!(flux.flux_up, metric_scaling)
    _scale_by_transpose!(flux.flux_dn, metric_scaling)
    _scale_by_transpose!(flux.flux_net, metric_scaling)
    flux isa FluxSW && _scale_by_transpose!(flux.flux_dn_dir, metric_scaling)
    return nothing
end

# Multiply the column-first `(ncol, nlev[, n_bnd])` array `a` by the
# vertical-first `(nlev, ncol)` scaling `sc`. Device arrays broadcast through
# the lazy transpose (a single kernel); host arrays use explicit loops, since
# broadcasts with a permuted-wrapper operand allocate a few bytes on
# Julia ≤ 1.11, tripping the zero-allocation tests.
_scale_by_transpose!(a::AbstractArray, sc) =
    (a .= a .* lazy_transpose(sc); nothing)
function _scale_by_transpose!(a::Array{FT, 2}, sc::Array{FT, 2}) where {FT}
    ncol, nlev = size(a)
    @inbounds for ilev in 1:nlev, gcol in 1:ncol
        a[gcol, ilev] *= sc[ilev, gcol]
    end
    return nothing
end
"""
    FluxPresentation{FTA2D, FD}

Host-facing `(nlev, ncol)` copies of a broadband flux set, filled by
[`update_presentation!`](@ref) from the column-first `(ncol, nlev)` compute
buffers at the end of every `update_lw_fluxes!`/`update_sw_fluxes!` (and hence
of [`update_fluxes!`](@ref RRTMGP.update_fluxes!)). The Layer-2 flux getters
return plain domain-masked views of these arrays, so the getter contract
(materializable with `Array`, broadcastable, reducible — also on the GPU, where
lazily transposed views of the compute buffers would fall outside the
wrapper types CUDA.jl dispatches on) holds without per-getter laziness.

# Fields
- `flux_up`: upward flux [W/m²] `(nlev, ncol)`.
- `flux_dn`: downward flux [W/m²] `(nlev, ncol)`.
- `flux_net`: net flux [W/m²] `(nlev, ncol)`.
- `flux_dn_dir`: direct downward flux [W/m²] `(nlev, ncol)`, or `nothing`
  (longwave).
"""
struct FluxPresentation{FTA2D <: AbstractArray, FD}
    flux_up::FTA2D
    flux_dn::FTA2D
    flux_net::FTA2D
    flux_dn_dir::FD
end
Adapt.@adapt_structure FluxPresentation

function FluxPresentation(grid_params::RRTMGPGridParams; direct::Bool)
    (; nlay, ncol) = grid_params
    FT = eltype(grid_params)
    DA = ClimaComms.array_type(grid_params)
    alloc() = DA{FT}(undef, nlay + 1, ncol)
    return FluxPresentation(
        alloc(),
        alloc(),
        alloc(),
        direct ? alloc() : nothing,
    )
end

"""
    update_presentation!(pres::FluxPresentation, flux::AbstractFlux)

Fill the `(nlev, ncol)` presentation arrays from the `(ncol, nlev)` compute
buffers of `flux` (transposing copies; the direct beam only when present).
"""
function update_presentation!(pres::FluxPresentation, flux::AbstractFlux)
    transpose_into!(pres.flux_up, flux.flux_up)
    transpose_into!(pres.flux_dn, flux.flux_dn)
    transpose_into!(pres.flux_net, flux.flux_net)
    isnothing(pres.flux_dn_dir) ||
        transpose_into!(pres.flux_dn_dir, flux.flux_dn_dir)
    return nothing
end

"""
    transpose_into!(dest, src)

Copy the 2D array `src` into `dest` with the two dimensions swapped
(`dest[i, j] = src[j, i]`). Device arrays use a single broadcast kernel; host
arrays use explicit loops, since broadcasts with a permuted-wrapper operand
allocate a few bytes on Julia ≤ 1.11, tripping the zero-allocation tests.
"""
transpose_into!(dest::AbstractArray, src) =
    (dest .= lazy_transpose(src); nothing)
function transpose_into!(dest::Array{FT, 2}, src::Array{FT, 2}) where {FT}
    n1, n2 = size(dest)
    @inbounds for j in 1:n2, i in 1:n1
        dest[i, j] = src[j, i]
    end
    return nothing
end

"""
    transpose_sum_into!(dest, a, b)

Set `dest[i, j] = a[j, i] + b[j, i]` (transposing sum of two 2D arrays); the
device/host split follows [`transpose_into!`](@ref).
"""
transpose_sum_into!(dest::AbstractArray, a, b) =
    (dest .= lazy_transpose(a) .+ lazy_transpose(b); nothing)
function transpose_sum_into!(
    dest::Array{FT, 2},
    a::Array{FT, 2},
    b::Array{FT, 2},
) where {FT}
    n1, n2 = size(dest)
    @inbounds for j in 1:n2, i in 1:n1
        dest[i, j] = a[j, i] + b[j, i]
    end
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
    # The band buffers share `metric_scaling`'s vertical-first layout, so the
    # `(nlev, ncol)` scaling broadcasts across the band dimension directly.
    band.flux_up .= band.flux_up .* metric_scaling
    band.flux_dn .= band.flux_dn .* metric_scaling
    return nothing
end

end
