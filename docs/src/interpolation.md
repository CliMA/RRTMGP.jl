# Level interpolation and extrapolation

The RTE solvers need pressures and temperatures both at layer centers, where
the optical properties are computed, and at the level faces, where the fluxes
are defined and the longwave source varies. A host that carries only
layer-center values asks RRTMGP to reconstruct the faces: the `interpolation`
and `bottom_extrapolation` options of [`RRTMGPSolver`](@ref
RRTMGP.RRTMGPSolver) and [`solve`](@ref RRTMGP.solve) select a scheme, and
[`interpolate_levels!`](@ref RRTMGP.interpolate_levels!) applies it inside
every [`update_fluxes!`](@ref RRTMGP.update_fluxes!) call, so a driver that
marches the layer temperatures keeps consistent faces automatically (the
[radiative-convective equilibrium tutorial](tutorials/manabe_rce.md) works
this way). This page derives the schemes and gives guidance on choosing one.

## The underlying model

Every scheme treats the atmosphere between two known points ``(p_1, T_1)`` and
``(p_2, T_2)`` (two layer centers, or a layer center and the surface) as a
hydrostatic ideal-gas column with a constant lapse rate, so the temperature is
linear in height:

```math
T(z) = T_1 + \frac{T_2 - T_1}{z_2 - z_1}\,(z - z_1).
```

Hydrostatic balance, ``dp/p = -g\,dz/(R_d T)``, integrates under this linear
``T(z)`` to a power law between pressure and temperature,

```math
p(T) = p_1 \left(\frac{T}{T_1}\right)^{B},
\qquad
B = \frac{\ln(p_2/p_1)}{\ln(T_2/T_1)},
```

where requiring the curve to pass through both endpoints fixes the exponent
``B`` without knowing ``g``, ``R_d``, or the lapse rate individually. In the
isothermal limit ``T_1 = T_2``, the power law degenerates and the pressure
instead decays exponentially with height,
``p(z) = p_1\, e^{-(z - z_1)/H}`` with scale height ``H = R_d T_1 / g``.
Fitting the two endpoints eliminates ``H`` and gives the form used in the
code,

```math
p(z) = p_1 \left(\frac{p_2}{p_1}\right)^{(z - z_1)/(z_2 - z_1)},
```

which is the same exponential: the base ``p_2/p_1`` is constant and the
exponent is linear in ``z``.

The schemes differ only in what they assume about where the face lies between
the two layer centers, and in whether they use the true face altitude.

## Interior faces

For a face between the adjacent layer centers ``(p_1, T_1)`` (below) and
``(p_2, T_2)`` (above):

| Scheme | Assumption | Face values |
|:-------|:-----------|:------------|
| `NoInterpolation` | levels are supplied by the host | none computed (default) |
| `ArithmeticMean` | none | ``T = (T_1 + T_2)/2``, ``p = (p_1 + p_2)/2`` |
| `GeometricMean` | face midway in ``\ln p`` | ``T = \sqrt{T_1 T_2}``, ``p = \sqrt{p_1 p_2}`` |
| `UniformZ` | face midway in height | ``T = (T_1 + T_2)/2``, ``p = p(T)`` from the power law |
| `UniformP` | face midway in pressure | ``p = (p_1 + p_2)/2``, ``T = p^{-1}(p)`` by inverting the power law |
| `BestFit` | face at its true altitude ``z`` | ``T(z)`` from the linear profile, ``p = p(T)`` from the power law |

`ArithmeticMean` averages pressure and temperature independently; the pair is
not exactly hydrostatic, but the scheme is the simplest and needs no
assumptions. The other schemes place the face on the constant-lapse-rate
curve, so their ``(p, T)`` pairs are hydrostatically consistent: under the
power law, ``\ln T`` is linear in ``\ln p``, so the geometric means of
`GeometricMean` are exactly the power-law state at the log-pressure midpoint.
`BestFit` is the only interior scheme that uses the actual face position
rather than assuming one; it requires the altitudes `center_z` and `face_z`
at construction.

## Boundary faces

At the top and bottom of the column, the face lies outside the layer centers
and the schemes extrapolate instead. With ``(p^+, T^+)`` the nearest layer and
``(p^{++}, T^{++})`` the next one in, `ArithmeticMean` and `UniformZ`
extrapolate the temperature linearly, ``T = (3T^+ - T^{++})/2``, and
`GeometricMean` and `UniformP` do the same in logarithms; each scheme then
completes the ``(p, T)`` pair the same way as in its interior rule.

The bottom face can instead be tied to the ground through
`bottom_extrapolation`:

- `SameAsInterpolation` (default): extrapolate like any boundary face. The
  bottom-face air temperature is then set by the atmosphere alone and may
  differ from the ground temperature, which enters separately through the
  surface-emission boundary condition. This difference is physical: pure
  radiative equilibrium produces air at the surface colder than the ground
  beneath it (see the
  [radiative-convective equilibrium tutorial](tutorials/manabe_rce.md)).
- `UseSurfaceTempAtBottom`: set the bottom-face air temperature to the ground
  temperature ``T_s``, the limit of strong turbulent heat exchange at the
  surface. The pressure follows the dry isentrope through the first layer,
  ``p = p^+ (T_s/T^+)^{c_p/R_d}``.
- `HydrostaticBottom`: extend the first layer downward at the dry-adiabatic
  lapse rate, ``T = T^+ + (g/c_p)(z^+ - z)``, with the same isentropic
  pressure; requires altitudes.

## Choosing a scheme

- **The host has level values**: keep the default `NoInterpolation` and write
  them through the getters.
- **Layer centers only, no altitudes**: `ArithmeticMean` is the simplest
  choice and adequate for smooth profiles. Use `GeometricMean`, `UniformZ`,
  or `UniformP` when hydrostatically consistent face pairs matter; they
  differ only in the assumed face position, and the differences shrink with
  layer thickness.
- **Layer centers and altitudes**: `BestFit` places each face where it
  actually is and is the most accurate of the schemes.
- **Bottom face**: keep the default when the air-ground temperature
  difference matters (radiative equilibrium, weak surface coupling); use
  `UseSurfaceTempAtBottom` to remove it (strong surface coupling);
  `HydrostaticBottom` extends the column dry-adiabatically when altitudes are
  available.

The interpolation runs on the device inside `update_fluxes!` on every call;
[`interpolate_levels!`](@ref RRTMGP.interpolate_levels!) can also be called on
its own. See [`AbstractInterpolation`](@ref RRTMGP.AbstractInterpolation) and
[`AbstractBottomExtrapolation`](@ref RRTMGP.AbstractBottomExtrapolation) for
the type reference.
