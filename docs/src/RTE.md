# Radiatve Transfer Equations

## Radiation heating rate

Radiative fluxes enter the thermal energy equation as follows:

```math
\begin{align}
\partial_t (\rho e) = - \nabla \cdot F_{net}
\end{align}
```

where ``F_{net} = \sum_j (F^{+}_j - F^{-}_j)``, and ``F^{+}_j`` and ``F^{-}_j`` denote the sum of shortwave and longwave upward and downward fluxes for spectral interval j, respectively.

## RTE solver

`solve_lw!` and `solve_sw!` solve the radiative transfor equation for longwave and shortwave radiation.

!!! note
    We have both no scattering and two stream solvers in the code. Only two stream solver is documented here.

## Computing fluxes through a vertically layered atmosphere

`rte_lw_2stream!` and `rte_sw_2stream!` compute longwave and shortwave radiative fluxes through a vertically layered atmosphere. The equations for upward and downward fluxes are calculated following [shonk2008](@cite):

```math
\begin{align}
F^{+}_i = \alpha_i F^{-}_i + G_i \\
F^{-}_i = \beta_i (T_i F^{-}_{i+1} + R_i G_i + S^{-}_i) + F^{-}_{i, dir} \\
\alpha_{i+1} = R_i + T_i^2 \beta_i \alpha_i \\
G_{i+1} = S^{+}_i + T_i \beta_i (G_i + \alpha_i S^{-}_i) \\
\end{align}
```

where ``\beta_k = 1 / (1 - \alpha_k R_k)``. ``T`` and ``R`` are transmittance and reflectance, and ``S^{+}`` and ``S^{-}`` represent upward and downward source functions. ``F^-_{dir}`` is the downward direct shortwave flux. The subscript ``i`` represents the level.

## Computing transmittance, reflectance, and total source functions

### Longwave

`lw_2stream_coeffs` calculates transmittance (``T``), reflectance (``R``), and source functions (``S^+`` and ``S^-``) from optical properties for longwave radiation. Transmittance and reflectance are calculated following [meador1980](@cite):

```math
\begin{align}
T = \frac{2 k e^{-k \tau}}{k + \gamma 1 + (k - \gamma_1) e^{-2 k \tau}} \\
R = \frac{\gamma_2 (1 - e^{-2 k \tau})}{k + \gamma 1 + (k - \gamma_1) e^{-2 k \tau}}
\end{align}
```

where ``\tau`` is the optical depth, and ``k = \sqrt{\gamma_1^2 - \gamma_2^2}``. ``\gamma_1`` and ``\gamma_2`` are coupling coefficients in the two-stream approximation that describe how isotropic the phase function is, and are determined by the optical properties (single scattering albedo and asymmetry parameter).

``S^+`` and ``S^-`` are calculated from optical properties and Planck source functions following [toon1989](@cite).

### Shortwave
`sw_2stream_coeffs` calculates diffuse transmittance (``T``), diffuse reflectance (``R``), direct transmittance (``T_{dir}``), and direct reflectance (``R_{dir}``) from optical properties for shortwave radiation. ``T`` and ``R`` are calculated following the same equations as for the longwave radiation. ``T_{dir}`` and ``R_{dir}`` is calculated following [meador1980](@cite):

```math
\begin{align}
T_{dir} = - A ((1 + k \mu) (\alpha_1 + k \gamma_4) e^{-\tau/\mu} - (1 - k \mu) (\alpha_1 - k \gamma_4) e^{-2k\tau} e^{-\tau/\mu} - 2k (\gamma_4 + \alpha_1 \mu) e^{-\tau/\mu}) \\
R_{dir} = A ((1 - k \mu) (\alpha_2 + k \gamma_3) - (1 + k \mu) (\alpha_2 - k \gamma_3) e^{-2 k \tau} - 2 (k \gamma_3 - \alpha_2 k \mu) e^{-k \tau}  e^{-\tau/\mu})
\end{align}
```

where ``A = \omega_0 / (1 - k^2 \mu^2) / (k (1 + e^{-2 k \tau}) + \gamma_1 (1 - e^{-2 k \tau}))``. ``\omega_0`` is the single scattering albedo. ``\mu`` is the cosine of solar zenith angle. ``\gamma_3``, and ``\gamma_4`` are coefficients in the two-stream approximation and determined by the optical properties. They are constrained to ``\gamma_3 + \gamma_4 = 1`` by energy conservation.

The direct downward flux (``F^{-}_{dir}``) and source functions (``S^+`` and ``S^-``) for each level are calculated from ``T_{dir}`` and ``R_{dir}``:

```math
\begin{align}
F^-_{dir, i} = e^{-\tau/\mu} F^-_{dir, i+1} \\
S^+_i = R_{dir, i} F^-_{dir, i+1} \\
S^-_i = T_{dir, i} F^-_{dir, i+1}
\end{align}
```
