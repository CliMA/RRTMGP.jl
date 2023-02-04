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

## Gas optics

RRTMGP calculates gas optics using a correlated-k method, in which the integration over frequency ``\nu`` of a complex line spectrum is replaced by an integration over a much smoother function of a new variable ``g``. The integral is approximated by the sum over ``G`` quadrature points (g-points). The mapping from ``\nu`` to ``g`` is computed for a set of longwave and shortwave bands. Within each band, the absorption is dominated by no more two gases (major species). The bands are non-overlapping, contiguous, and span the frequencies of radiation emitted by the Earth or Sun. For each band, the total optical thickness (``\tau``) includes contributions from major species (``\tau_\mathrm{major}``), minor species (``\tau_\mathrm{minor}``), and Rayleigh scattering (``\tau_\mathrm{rayleigh}``) for shortwave bands.

`compute_τ_ssa_lw_src!` calculates the optical properties for each g-point in the longwave and shortwave bands, as well as source functions for the longwave bands. The input atmospheric conditions include pressure (``p``), temperature (``T``), and volume mixing ratios of gases (``\chi``). The relationship between the volume mixing ratio and the mass fraction (``q``) of a gas is ``\chi = q / q_d * M_d / M``, where ``q_d = 1-q_t`` is the dry air mass fraction (total specific humidity ``q_t``), ``M_d`` is the molar mass of the dry air, and `M` is the molar mass of the gas. 

For the major species, the absorption coefficient is computed by linearly interpolating the tabulated values in ``\log(p)``, ``T``, and relative abundances of the two major species (``\eta``) calculated from the volume mixing ratios. The contributions from minor species and Rayleigh scattering are treated with less detail. For each minor species and Rayleigh scattering, a representative pressure is chosen, and the absorption coefficient is computed by linearly interpolating the tabulated values in ``T`` and ``\eta``. For the shortwave bands, the single scattering albedo is calculated as ``\omega_0 = \tau_\mathrm{rayleigh} / \tau``, and the asymmetry parameter is zero.

For longwave bands, the Planck fraction (defined as the fraction of the band-integrated Planck energy associated with each g-point) is computed by linearly interpolating the tabulated values in ``\log(p)``, ``T``, and ``\eta``. The band-integrated Planck function is uniquely determined by temperature, and is computed by linearly interpolating the tabulated values in ``T``. The shortwave source function for each g-point is assumed to be constant.

### Lookup tables
The lookup tables for the spectral map between ``\nu`` and ``g`` are stored in netcdf files [here](https://caltech.box.com/shared/static/wbtrwp44dyn08g7mozjf4fcyrexwbe6a.gz). Longwave and shortwave lookup tables are stored in `clearsky_lw.nc` and `clearsky_sw.nc`, respectively. `LookUpLW` and `LookUpSW` read the lookup tables. The tabulated information includes the spectral distretization (bands), the major and minor gases considered in each band, and tabulated data of absorption coefficients, Planck fraction, and band-integrated Planck function. The absorption coefficients and Planck fraction are computed at pressures ``1 \mathrm{Pa} ≤ p ≤ 109663 \mathrm{Pa}`` in increments of ``\log(p/1 \mathrm{Pa}) = 0.2``, temperatures ``160 \mathrm{K} ≤ T ≤ 355 \mathrm{K}`` in ``15 \mathrm{K}`` increments, and ``0 ≤ \eta ≤ 1`` in ``1/8`` increments. The band-integrated Plack function is computed at temperatures ``160 \mathrm{K} ≤ T ≤ 355 \mathrm{K}`` in ``1 \mathrm{K}`` increments. A more detailed description of the spectral structure of the lookup tables can be found in `Appendix A` in [pincus2019](@cite).