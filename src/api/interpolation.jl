#####
##### Level <-> layer interpolation and bottom extrapolation
#####
#
# Methods for obtaining cell-face (level) pressures/temperatures from cell-center
# (layer) pressures/temperatures, used when a caller provides only layer values.

# TODO: Remove comment after ClimaAtmos update
#
# Lifted from the ClimaAtmos RRTMGP interface so RRTMGP can do this itself,
# without ClimaCore. The arithmetic is kept identical to preserve bit-for-bit
# reproducibility of downstream fluxes.

# Each scheme assumes a constant lapse rate ∂T/∂z = L between the two known points
# (z₁, p₁, T₁) and (z₂, p₂, T₂), giving T(z) = T₁ + (T₂ - T₁) / (z₂ - z₁) * (z - z₁).
# Pressure follows hydrostatic balance when T₁ == T₂ and an isentropic p ∝ T^B law
# otherwise:
#   p(z) = T₁ == T₂ ? p₁ * (p₂ / p₁)^((z - z₁) / (z₂ - z₁))
#                   : p₁ * (p₂ / p₁)^(log(T(z) / T₁) / log(T₂ / T₁)).
# UniformZ, UniformP, and GeometricMean are the midpoint-in-height,
# midpoint-in-pressure, and z = z₁ + (z₂ - z₁) / (√(T₂/T₁) + 1) special cases that
# avoid needing altitudes; ArithmeticMean takes the plain average of p and T. The
# full derivation is a candidate for a docs Explanation page alongside the
# "Functional core" page.

"""
    AbstractInterpolation

Strategy for obtaining cell-face (level) pressures and temperatures from
cell-center (layer) values, or vice versa, used when only one of the two is
provided. The scheme interpolates on interior faces and extrapolates on boundary
faces.

Subtypes:
- `NoInterpolation`: levels are supplied directly; do not interpolate.
- `ArithmeticMean`: arithmetic mean of the two adjacent layers.
- `GeometricMean`: geometric mean of the two adjacent layers.
- `UniformZ`: assume the face lies midway in height between the layers.
- `UniformP`: assume the face lies midway in pressure between the layers.
- `BestFit`: constant-lapse-rate fit using layer altitudes; requires `center_z` and `face_z`.

The derivations are summarized in the comment block above.
"""
abstract type AbstractInterpolation end
struct NoInterpolation <: AbstractInterpolation end
struct ArithmeticMean <: AbstractInterpolation end
struct GeometricMean <: AbstractInterpolation end
struct UniformZ <: AbstractInterpolation end
struct UniformP <: AbstractInterpolation end
struct BestFit <: AbstractInterpolation end

# UseSurfaceTempAtBottom and HydrostaticBottom assume a dry ideal gas in an
# isentropic process, p(z) = p⁺ * (T(z) / T⁺)^(cₚ / R), where (p⁺, T⁺) is the first
# layer above the bottom face. UseSurfaceTempAtBottom sets T(z) = Tₛ (the surface
# temperature); HydrostaticBottom uses the dry-adiabatic lapse rate, giving
# T(z) = T⁺ + (g / cₚ) * (z⁺ - z). The full derivation is a candidate for a docs
# Explanation page alongside the "Functional core" page.

"""
    AbstractBottomExtrapolation

Strategy for obtaining the bottom cell-face (level) pressure and temperature from
the layer values above it.

Subtypes:
- `SameAsInterpolation`: extrapolate using the interpolation scheme.
- `UseSurfaceTempAtBottom`: set the bottom-face air temperature to the surface temperature.
- `HydrostaticBottom`: assume a dry-adiabatic lapse rate; requires `center_z` and `face_z`.
"""
abstract type AbstractBottomExtrapolation end
struct SameAsInterpolation <: AbstractBottomExtrapolation end
struct UseSurfaceTempAtBottom <: AbstractBottomExtrapolation end
struct HydrostaticBottom <: AbstractBottomExtrapolation end

requires_z(::Any) = false
requires_z(::Union{BestFit, HydrostaticBottom}) = true

uniform_z_p(T, p₁, T₁, p₂, T₂) =
    T₁ == T₂ ? sqrt(p₁ * p₂) : p₁ * (p₂ / p₁)^(log(T / T₁) / log(T₂ / T₁))
best_fit_p(T, z, p₁, T₁, z₁, p₂, T₂, z₂) =
    T₁ == T₂ ? p₁ * (p₂ / p₁)^((z - z₁) / (z₂ - z₁)) :
    p₁ * (p₂ / p₁)^(log(T / T₁) / log(T₂ / T₁))

function interp!(::ArithmeticMean, p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
    @. T = (Tꜜ + Tꜛ) / 2
    @. p = (pꜜ + pꜛ) / 2
end
function interp!(::GeometricMean, p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
    @. T = sqrt(Tꜜ * Tꜛ)
    @. p = sqrt(pꜜ * pꜛ)
end
function interp!(::UniformZ, p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
    @. T = (Tꜜ + Tꜛ) / 2
    @. p = uniform_z_p(T, pꜜ, Tꜜ, pꜛ, Tꜛ)
end
function interp!(::UniformP, p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
    @. p = (pꜜ + pꜛ) / 2
    @. T = Tꜜ * (Tꜛ / Tꜜ)^(log(p / pꜜ) / log(pꜛ / pꜜ)) # assume that pꜜ != pꜛ
end
function interp!(::BestFit, p, T, z, pꜜ, Tꜜ, zꜜ, pꜛ, Tꜛ, zꜛ)
    @. T = Tꜜ + (Tꜛ - Tꜜ) * (z - zꜜ) / (zꜛ - zꜜ)
    @. p = best_fit_p(T, z, pꜜ, Tꜜ, zꜜ, pꜛ, Tꜛ, zꜛ)
end

function extrap!(::ArithmeticMean, p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
    @. T = (3 * T⁺ - T⁺⁺) / 2
    @. p = (3 * p⁺ - p⁺⁺) / 2
end
function extrap!(::GeometricMean, p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
    @. T = sqrt(T⁺^3 / T⁺⁺)
    @. p = sqrt(p⁺^3 / p⁺⁺)
end
function extrap!(::UniformZ, p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
    @. T = (3 * T⁺ - T⁺⁺) / 2
    @. p = uniform_z_p(T, p⁺, T⁺, p⁺⁺, T⁺⁺)
end
function extrap!(::UniformP, p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
    @. p = (3 * p⁺ - p⁺⁺) / 2
    @. T = T⁺ * (T⁺⁺ / T⁺)^(log(p / p⁺) / log(p⁺⁺ / p⁺)) # assume that p⁺ != p⁺⁺
end
function extrap!(::BestFit, p, T, z, p⁺, T⁺, z⁺, p⁺⁺, T⁺⁺, z⁺⁺, Tₛ, params)
    @. T = T⁺ + (T⁺⁺ - T⁺) * (z - z⁺) / (z⁺⁺ - z⁺)
    @. p = best_fit_p(T, z, p⁺, T⁺, z⁺, p⁺⁺, T⁺⁺, z⁺⁺)
end
function extrap!(::UseSurfaceTempAtBottom, p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
    cₚ = RP.cp_d(params)
    R = RP.R_d(params)
    @. T = Tₛ
    @. p = p⁺ * (T / T⁺)^(cₚ / R)
end
function extrap!(
    ::HydrostaticBottom,
    p,
    T,
    z,
    p⁺,
    T⁺,
    z⁺,
    p⁺⁺,
    T⁺⁺,
    z⁺⁺,
    Tₛ,
    params,
)
    FT = eltype(p)
    g = FT(RP.grav(params))
    cₚ = FT(RP.cp_d(params))
    R = FT(RP.R_d(params))
    @. T = T⁺ + g / cₚ * (z⁺ - z)
    @. p = p⁺ * (T / T⁺)^(cₚ / R)
end
