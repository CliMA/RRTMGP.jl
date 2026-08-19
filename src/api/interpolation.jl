#####
##### Level <-> layer interpolation and bottom extrapolation
#####
#
# Methods for obtaining cell-face (level) pressures/temperatures from cell-center
# (layer) pressures/temperatures, used when a caller provides only layer values.

# Each scheme assumes a constant lapse rate ∂T/∂z = L between the two known points
# (z₁, p₁, T₁) and (z₂, p₂, T₂), giving T(z) = T₁ + (T₂ - T₁) / (z₂ - z₁) * (z - z₁).
# Pressure follows hydrostatic balance when T₁ == T₂ and the constant-lapse-rate
# hydrostatic power law p ∝ T^B otherwise:
#   p(z) = T₁ == T₂ ? p₁ * (p₂ / p₁)^((z - z₁) / (z₂ - z₁))
#                   : p₁ * (p₂ / p₁)^(log(T(z) / T₁) / log(T₂ / T₁)).
# UniformZ, UniformP, and GeometricMean are the midpoint-in-height,
# midpoint-in-pressure, and z = z₁ + (z₂ - z₁) / (√(T₂/T₁) + 1) special cases that
# avoid needing altitudes; ArithmeticMean takes the plain average of p and T. The
# full derivation is on the "Level interpolation" docs page
# (docs/src/interpolation.md).

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

The derivations, and guidance on choosing a scheme, are on the
"Level interpolation" docs page.
"""
abstract type AbstractInterpolation end

"""
    NoInterpolation()

Take the level (cell-face) pressures and temperatures as the caller supplied
them; interpolate nothing.
"""
struct NoInterpolation <: AbstractInterpolation end

"""
    ArithmeticMean()

Set each interior face to the arithmetic mean of its two adjacent layers, in
both pressure and temperature, and extrapolate linearly at boundary faces.
"""
struct ArithmeticMean <: AbstractInterpolation end

"""
    GeometricMean()

Set each interior face to the geometric mean of its two adjacent layers, in both
pressure and temperature (the constant-lapse-rate state at the log-pressure
midpoint), and extrapolate logarithmically at boundary faces.
"""
struct GeometricMean <: AbstractInterpolation end

"""
    UniformZ()

Place each interior face midway in height between its two adjacent layers:
average the temperature, then take the pressure from the constant-lapse-rate
power law.
"""
struct UniformZ <: AbstractInterpolation end

"""
    UniformP()

Place each interior face midway in pressure between its two adjacent layers:
average the pressure, then invert the constant-lapse-rate power law for the
temperature.
"""
struct UniformP <: AbstractInterpolation end

"""
    BestFit()

Place each face at its true altitude on the constant-lapse-rate profile through
the two adjacent layers. Requires `center_z` and `face_z`.
"""
struct BestFit <: AbstractInterpolation end

# UseSurfaceTempAtBottom and HydrostaticBottom assume a dry ideal gas in an
# isentropic process, p(z) = p⁺ * (T(z) / T⁺)^(cₚ / R), where (p⁺, T⁺) is the first
# layer above the bottom face. UseSurfaceTempAtBottom sets T(z) = Tₛ (the surface
# temperature); HydrostaticBottom uses the dry-adiabatic lapse rate, giving
# T(z) = T⁺ + (g / cₚ) * (z⁺ - z). The full derivation is on the
# "Level interpolation" docs page (docs/src/interpolation.md).

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

"""
    SameAsInterpolation()

Extrapolate the bottom face with the same scheme used for the other boundary
faces, so its air temperature is set by the atmosphere alone and may differ from
the ground temperature.
"""
struct SameAsInterpolation <: AbstractBottomExtrapolation end

"""
    UseSurfaceTempAtBottom()

Set the bottom-face air temperature to the surface temperature (the limit of
strong turbulent heat exchange at the surface) and take the pressure from the
dry isentrope through the first layer.
"""
struct UseSurfaceTempAtBottom <: AbstractBottomExtrapolation end

"""
    HydrostaticBottom()

Extend the first layer down to the bottom face at the dry-adiabatic lapse rate,
with the isentropic pressure. Requires `center_z` and `face_z`.
"""
struct HydrostaticBottom <: AbstractBottomExtrapolation end

"""
    requires_z(scheme)

Return whether `scheme` needs layer and face altitudes (`center_z`, `face_z`) at
construction: `true` for [`BestFit`](@ref) and [`HydrostaticBottom`](@ref),
`false` for every other scheme.
"""
requires_z(::Any) = false
requires_z(::Union{BestFit, HydrostaticBottom}) = true

"""
    uniform_z_p(T, p₁, T₁, p₂, T₂)

Return the pressure at temperature `T` on the constant-lapse-rate hydrostatic
power law through `(p₁, T₁)` and `(p₂, T₂)`, degenerating to the geometric mean
of the pressures in the isothermal limit `T₁ == T₂`.
"""
uniform_z_p(T, p₁, T₁, p₂, T₂) =
    T₁ == T₂ ? sqrt(p₁ * p₂) : p₁ * (p₂ / p₁)^(log(T / T₁) / log(T₂ / T₁))
"""
    best_fit_p(T, z, p₁, T₁, z₁, p₂, T₂, z₂)

Return the pressure at temperature `T` and altitude `z` on the
constant-lapse-rate hydrostatic power law through `(p₁, T₁, z₁)` and
`(p₂, T₂, z₂)`, using the altitudes in the isothermal limit `T₁ == T₂`.
"""
best_fit_p(T, z, p₁, T₁, z₁, p₂, T₂, z₂) =
    T₁ == T₂ ? p₁ * (p₂ / p₁)^((z - z₁) / (z₂ - z₁)) :
    p₁ * (p₂ / p₁)^(log(T / T₁) / log(T₂ / T₁))

"""
    interp!(scheme, p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
    interp!(::BestFit, p, T, z, pꜜ, Tꜜ, zꜜ, pꜛ, Tꜛ, zꜛ)

Fill the interior-face pressures `p` and temperatures `T` from the layers below
(`pꜜ`, `Tꜜ`) and above (`pꜛ`, `Tꜛ`) according to `scheme`. Add a method to
support a new [`AbstractInterpolation`](@ref).
"""
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

"""
    extrap!(scheme, p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
    extrap!(::BestFit, p, T, z, p⁺, T⁺, z⁺, p⁺⁺, T⁺⁺, z⁺⁺, Tₛ, params)

Fill a boundary-face pressure `p` and temperature `T` from the nearest layer
(`p⁺`, `T⁺`), the next one in (`p⁺⁺`, `T⁺⁺`), the surface temperature `Tₛ`, and
`params`, according to `scheme`. Add a method to support a new
[`AbstractInterpolation`](@ref) or [`AbstractBottomExtrapolation`](@ref).
"""
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
