# RRTMGP Optics

`compute_optical_props!` computes the optical properties and source functions given an atmospheric state. For the two-stream approximation, the optical properties include the optical thickness of the atmosphere layer (``\tau``), single-scattering albedo (``\omega_0``), and asymmetry parameter (``g``).

## Gray atmosphere optics

The gray atmosphere (more precisely, semi-gray atmosphere) approximates the radiative properties of the atmosphere by taking the optical properties to be independent of wavelength, separately in the longwave and shortwave bands. 

### Longwave

The optical depth (``d``) for longwave radiation follows [schneider2004](@cite):

```math
\begin{align}
d(\phi, p) = d_0(\phi) \left(\frac{p}{p_0}\right)^\alpha
\end{align}
```

Here, ``\alpha``  is ratio of the pressure scale height to the partial-pressure scale height of the infrared absorber (``\alpha=3.5`` for water vapor); ``\phi`` is the latitude; ``p`` is the pressure, and ``p_0 = 1000~\mathrm{hPa}`` is the reference pressure. The optical depth at ``\phi`` and ``p_0``, ``d_0``, is calculated from the latitude-dependent radiative equilibrium temperature ``T_s(\phi)``:

```math
\begin{align}
T_s(\phi) = T_e + \Delta T \left(\frac{1}{3} - \sin^2\phi\right) \\
d_0(\phi) = \left(\frac{T_s(\phi)}{T_t}\right)^4 - 1.
\end{align}
```

Here, ``T_e`` is the global-mean surface temperature in radiative equilibrium, ``T_t`` is the temperature at the top of the atmosphere, and ``\Delta T`` is the equator-to-pole temperature difference in radiative equilibrium. The default values are ``T_e = 300~\mathrm{K}``, ``T_t = 200~\mathrm{K}``, and ``\Delta T = 60~\mathrm{K}``.

The optical thickness of an atmosphere layer (the differential optical depth) of pressure thickness ``\Delta p`` is
```math
\begin{align}
\tau(\phi, p) = \alpha d_0(\phi) \left(\frac{p}{p_0}\right)^\alpha \frac{\Delta p}{p}.
\end{align}
```

The source function for longwave radiation is calculated as ``S = \sigma T^4 / \pi``, where ``T`` is the air temperature and ``\sigma`` is the Stefan–Boltzmann constant.

### Shortwave

The optical depth (``d``) for shortwave radiation follows [frierson2007](@cite) and [ogorman2008](@cite):

```math
\begin{align}
d(p) = \tau_0 \left(\frac{p}{p_0}\right)^2
\end{align}
```

where ``p`` is the pressure of the atmosphere layer and ``p_0 = 1000~\mathrm{hPa}`` is the reference pressure. The default value for ``\tau_0`` is 0.22. 

The differential optical thickness of an atmosphere layer of pressure thickness ``\Delta p`` is

```math
\begin{align}
\tau(p) = 2 \tau_0 \frac{p}{p_0} \frac{\Delta p}{p_0}.
\end{align}
```


The single scattering albedo and asymmetry parameter are zero.