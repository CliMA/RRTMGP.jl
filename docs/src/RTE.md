# Radiative transfer

The RTE solvers integrate the radiative transfer equation through a
plane-parallel, horizontally homogeneous column and return the upward and
downward radiative fluxes and the heating rate. This page states the equations
they solve and how they solve them; the optical properties that enter as
coefficients are described under [Optics](@ref).

## The radiative transfer equation

Radiation is described by the spectral radiance ``I_\lambda``, the power per
unit area, solid angle, and wavelength flowing in a given direction ``\Omega``.
Along a ray, the medium changes the radiance by three processes: extinction
(absorption plus scattering out of the beam), thermal emission, and scattering
into the beam from all other directions. In terms of the extinction optical
depth ``\tau_\lambda`` measured along the ray, the radiative transfer equation
reads
```math
\frac{dI_\lambda(\Omega)}{d\tau_\lambda} =
(1 - \omega_{0,\lambda})\, B_\lambda(T) - I_\lambda(\Omega)
+ \omega_{0,\lambda} \int_{4\pi} P_\lambda(\Omega', \Omega)\,
I_\lambda(\Omega')\, d\Omega',
```
where ``B_\lambda(T)`` is the Planck function at the local temperature,
``\omega_{0,\lambda}`` is the single-scattering albedo (the fraction of the
extinction due to scattering), and ``P_\lambda(\Omega', \Omega)`` is the phase
function, the probability of scattering from direction ``\Omega'`` into
``\Omega``. The emission term ``B_\lambda`` follows from Kirchhoff's law under
local thermodynamic equilibrium (absorptivity = emissivity at every wavelength),
which is assumed to hold throughout the atmosphere RRTMGP models. The scattering
integral couples all directions, which is what makes the general equation
computationally demanding; in the atmosphere, scattering matters for clouds and
aerosols in both the longwave and shortwave and for Rayleigh scattering in the
shortwave. The two-stream approximation below reduces the angular coupling to
two hemispheric fluxes.

For longwave transfer in clear skies, scattering is negligible
(``\omega_{0,\lambda} \to 0``), and the equation reduces to the *Schwarzschild
equation*,
```math
\frac{dI_\lambda}{d\tau_\lambda} = B_\lambda(T) - I_\lambda,
```
with ``\tau_\lambda`` now the absorption optical depth.

## Two-stream approximation

Resolving the full angular dependence of ``I_\lambda`` is impractical in a
climate model, so the flux calculation collapses the radiation field into an
upward flux ``F_\lambda^\uparrow`` and a downward flux ``F_\lambda^\downarrow``.
In the plane-parallel column, the vertical optical depth
``\widehat{\tau}_\lambda`` increases downward from zero at the top of the
atmosphere. Neglecting scattering, the two streams obey
```math
\frac{dF_\lambda^\uparrow}{d\widehat{\tau}_\lambda} = D\,(F_\lambda^\uparrow - \pi B_\lambda),
\qquad
\frac{dF_\lambda^\downarrow}{d\widehat{\tau}_\lambda} = D\,(\pi B_\lambda - F_\lambda^\downarrow),
```
where the *diffusivity factor* ``D`` accounts for the slantwise paths of a
diffuse field: replacing the range of zenith angles by a single effective secant
``D`` closes the hemispheric integral that relates flux to radiance. With
scattering, the streams couple through a layer reflectance, and RRTMGP forms the
two-stream transmittance and reflectance from the optical thickness, single
scattering albedo, and asymmetry parameter following [meador1980](@citet).

## Angular discretization

The longwave no-scattering solver improves on the single diffusivity angle: it
integrates the Schwarzschild equation along a small set of discrete zenith
angles and sums the results with Gauss quadrature weights, so the hemispheric
flux is ``F_\lambda = \sum_i w_i\, I_\lambda(\mu_i)``. RRTMGP uses the
Gauss-Jacobi-5 nodes of [hogan2024](@citet) with one to four angles
([`AngularDiscretization`](@ref
RRTMGP.AngularDiscretizations.AngularDiscretization)). The default single angle
has secant ``D \approx 1.64``, close to Elsasser's classic diffusivity factor of
``1.66``, which the two-stream longwave solver adopts following
[fu1997](@citet); adding angles improves the accuracy of the angular integral.

## Radiative heating rate

The fluxes feed back on temperature through the divergence of the net upward
flux ``F^{\mathrm{net}} = F^\uparrow - F^\downarrow``. A layer warms where the
net flux converges and cools where it diverges:
```math
\rho\, c_p \frac{\partial T}{\partial t} = -\frac{dF^{\mathrm{net}}}{dz}.
```
Because the temperature tendency scales as ``1/\rho``, a given flux divergence
warms thin, high-altitude air more than dense air near the surface; this is one
reason ozone's shortwave absorption heats the stratosphere effectively.
[`heating_rate`](@ref RRTMGP.heating_rate) returns ``\partial T/\partial t`` in
K/s.

## How the equations are solved

The continuous equations are solved by discretizing the column into layers of
uniform optical properties, within which they admit exact solutions. The
resulting sweeps run independently for every spectral quadrature point and
column, which is what makes the solvers parallelize well.

The no-scattering longwave solver integrates the Schwarzschild equation along
each quadrature angle exactly: crossing a layer of vertical optical depth
``\widehat{\tau}`` multiplies the radiance by the slant-path transmittance
``t = e^{-D\widehat{\tau}}`` and adds the layer's own emission,
```math
I_{\mathrm{out}} = t\, I_{\mathrm{in}} + S,
```
where the source ``S`` follows from a Planck function that varies linearly in
optical depth across the layer ([clough1992](@citet)). One downward and one
upward sweep of this update per quadrature angle yields the fluxes.

In the two-stream solver, the coupled equations for ``F^\uparrow`` and
``F^\downarrow`` within a layer have exponential solutions ``e^{\pm k\tau}``,
which reduce each layer ``i`` to a diffuse reflectance ``R_i``, a transmittance
``T_i``, and upward and downward source terms ``S_i^{\pm}``
([meador1980](@citet); [toon1989](@citet) for the thermal source). The layers
are then coupled by the adding method of [shonk2008](@citet). With ``\alpha_i``
the albedo of the atmosphere–surface system below level ``i`` and ``G_i`` its
upward emission, a first pass climbs from the surface,
```math
\begin{aligned}
\alpha_{i+1} &= R_i + T_i^2 \beta_i \alpha_i,\\
G_{i+1} &= S_i^{+} + T_i \beta_i \left(G_i + \alpha_i S_i^{-}\right),
\end{aligned}
```
where ``\beta_i = (1 - \alpha_i R_i)^{-1}`` sums the infinite series of
reflections between layer ``i`` and the medium below it. A second pass descends
from the top-of-atmosphere boundary condition, computing the downwelling flux
```math
F_i^{\downarrow} = \beta_i \left(T_i F_{i+1}^{\downarrow} + R_i G_i +
S_i^{-}\right) + F_{i,\mathrm{dir}}^{\downarrow},
```
where the last term is the direct solar beam, attenuated by Beer's law, and from
it the upwelling flux ``F_i^{\uparrow} = \alpha_i F_i^{\downarrow} + G_i``. The
result is the discrete counterpart of the integral solution of the transfer
equation: the flux at each level accumulates every layer's emission, transmitted
and reflected through the layers between the emitting layer and that level.

## No-scattering and two-stream solvers

RRTMGP provides two solver families, distinguished by whether they represent
scattering:

- The **no-scattering solver** carries the optical thickness ``\tau`` alone
  ([`OneScalar`](@ref RRTMGP.Optics.OneScalar) optics). It integrates absorption
  and emission over the Gauss angles above and is used for clear-sky longwave
  transfer, where scattering is negligible; in the shortwave, it reduces to
  Beer's-law extinction of the direct solar beam.
- The **two-stream solver** carries the optical thickness, single scattering
  albedo, and asymmetry parameter ([`TwoStream`](@ref RRTMGP.Optics.TwoStream)
  optics) and represents multiple scattering. It is required wherever scattering
  is not negligible: clouds and aerosols in both bands, and Rayleigh scattering
  in the shortwave.

Both families exist for the longwave and the shortwave ([`NoScatLWRTE`](@ref
RRTMGP.RTE.NoScatLWRTE) / [`TwoStreamLWRTE`](@ref RRTMGP.RTE.TwoStreamLWRTE) and
[`NoScatSWRTE`](@ref RRTMGP.RTE.NoScatSWRTE) / [`TwoStreamSWRTE`](@ref
RRTMGP.RTE.TwoStreamSWRTE)), and [`solve_lw!`](@ref RRTMGP.RTESolver.solve_lw!)
/ [`solve_sw!`](@ref RRTMGP.RTESolver.solve_sw!) dispatch on the workspace type.

## Boundary conditions

The sweep needs conditions at both ends of the column. At the surface, the
longwave upwelling radiance is the surface Planck emission scaled by the surface
emissivity, and the shortwave reflection is set by the direct and diffuse
albedos; the incident solar flux and the cosine of the solar zenith angle enter
at the top. Any prescribed incident diffuse flux at the top defaults to zero, as
does the downwelling longwave flux at the top.

RRTMGP can also insert an **isothermal boundary layer**: an extra layer between
the host model's top and the minimum pressure of the gas-optics tables, held at
the model-top temperature, that represents the radiative effect of the
atmosphere above the model lid. It is enabled through [`RRTMGPGridParams`](@ref
RRTMGP.RRTMGPGridParams), and every getter masks it off so the host reads and
writes arrays sized to its own grid (see [The getter contract](@ref)).
