# RRTMGP Optics

`compute_optical_props!` computes the optical properties and source functions given an atmospheric state. For the two-stream approximation, the optical properties include the optical thickness of the atmosphere layer (``\tau``), single-scattering albedo (``\omega_0``), and asymmetry parameter (``g``).

## Gray atmosphere optics

The gray atmosphere (more precisely, semi-gray atmosphere) approximates the radiative properties of the atmosphere by taking the optical properties to be independent of wavelength, separately in the longwave and shortwave bands. 

### Longwave
Two options are currently supported for computing the optical depth for longwave radiation. The optical depth (``d``) for longwave radiation, with the "GrayOpticalThicknessSchneider2004" option, follows [schneider2004](@cite):

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


The optical depth (``d``) for longwave radiation, with the "GrayOpticalThicknessOGorman2008" option, follows [ogorman2008](@cite):

```math
\begin{align}
d(\phi, p) = \alpha \left[f_l \sigma + (1 - f_l) \sigma^4 \right] \left[ \tau_e + (\tau_p - \tau_e) \sin^2\phi \right]
\end{align}
```

where ``f_l = 0.2``, ``\sigma = p / p_0`` is pressure p normalized by
surface pressure ``p_0``, ``\phi`` is latitude, and the longwave optical thicknesses at the equator and at the pole are  $\tau_e = 7.2$ and $\tau_p = 1.8$, respectively. $\alpha$ is a scaling factor. 

The optical thickness of an atmosphere layer (the differential optical depth) of pressure thickness ``\Delta p`` is
```math
\begin{align}
\tau(\phi, p) = (\alpha * \frac{\Delta p}{p}) \left[f_l \sigma + 4 (1 - f_l) \sigma^4 \right] \left[ \tau_e + (\tau_p - \tau_e) \sin^2\phi \right]
\end{align}
```

The source function for longwave radiation is calculated as ``S = \sigma T^4 / \pi``, where ``T`` is the air temperature and ``\sigma`` is the Stefan–Boltzmann constant.

### Shortwave

Two options are currently supported for computing the optical depth for shortwave radiation. The optical depth (``d``) for shortwave radiation, with the "GrayOpticalThicknessSchneider2004" option, is 0.

The optical depth (``d``) for shortwave radiation, with the "GrayOpticalThicknessOGorman2008" option, follows [frierson2007](@cite) and [ogorman2008](@cite):

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
The lookup tables for the spectral map between ``\nu`` and ``g`` are stored in netcdf files [here](https://caltech.box.com/shared/static/wbtrwp44dyn08g7mozjf4fcyrexwbe6a.gz). Longwave and shortwave lookup tables are stored in `clearsky_lw.nc` and `clearsky_sw.nc`, respectively. `LookUpLW` and `LookUpSW` read the lookup tables. The tabulated information includes the spectral discretization (bands), the major and minor gases considered in each band, and tabulated data of absorption coefficients, Planck fraction, and band-integrated Planck function. The absorption coefficients and Planck fraction are computed at pressures ``1 \mathrm{Pa} ≤ p ≤ 109663 \mathrm{Pa}`` in increments of ``\log(p/1 \mathrm{Pa}) = 0.2``, temperatures ``160 \mathrm{K} ≤ T ≤ 355 \mathrm{K}`` in ``15 \mathrm{K}`` increments, and ``0 ≤ \eta ≤ 1`` in ``1/8`` increments. The band-integrated Plack function is computed at temperatures ``160 \mathrm{K} ≤ T ≤ 355 \mathrm{K}`` in ``1 \mathrm{K}`` increments. A more detailed description of the spectral structure of the lookup tables can be found in `Appendix A` in [pincus2019](@cite).

!!! note
    Absorption by both major and minor species is treated separately in the upper and lower atmosphere, defined as pressures above and below the pressure defined in `press_ref_trop` in the lookup tables (~9948 Pa). This corresponds to element 13 in array `press_ref`. There is a single large table of absorption coefficients. Table elements 1:13 in the pressure dimension correspond to pressures >= `press_ref_trop`, while table elements 14:60 correspond to pressures <= `press_ref_trop`. Elements 13 and 14 refer to the same pressure but for the lower and upper atmosphere, respectively. This is why the the absorption coefficients lookup table has a pressure dimension of 60, but there are only 59 reference pressure levels in `press_ref`.

## Cloud optics

`compute_cld_props` computes the optical properties of clouds, which is a combination of contributions from liquid and ice particles. The optical thickness (``\tau``), single scattering albedo (``\omega_0``), and asymmetry parameter (``g``) for clouds are:

```math
\begin{align}
\tau_c &= \tau_l + \tau_i \\
\omega_{0c} &= \frac{\tau_l \omega_{0l} + \tau_i \omega_{0i}}{\tau_l + \tau_i} \\
g_c &= \frac{\tau_l \omega_{0l} g_l + \tau_i \omega_{0i} g_i}{\tau_l \omega_{0l} + \tau_i \omega_{0i}} \\
\end{align}
```

where subscripts c, l, and i indicate cloud, liquid, and ice, respectively. The liquid cloud and ice cloud optical thicknesses are calculated by multiplying the cloud liquid and ice extinction coefficients (``\mu``) by the cloud liquid (lwp) and ice water paths (iwp), respectively. That is, ``\tau_l = \mu_l \mathrm{lwp}`` and ``\tau_i = \mu_i \mathrm{iwp}``. The extinction coefficient, single scattering albedo, and asymmetry parameter are computed by linearly interpolating the tabulated values in the liquid and ice particle effective radius. For ice partices, the optical properties also depend on the surface roughness, following [yang2013](@cite). The optical properties for shortwave bands are delta-scaled to reduce the biases in calculating radiative fluxes for highly asymmetric phase functions ([joseph1976](@cite)).

`add_cloud_optics_2stream` adds cloud optics to the gas optics for each g-point that sees cloud (the cloud mask is described below).

### Lookup tables
The lookup tables for cloud optics are stored in netcdf files [here](https://caltech.box.com/shared/static/wbtrwp44dyn08g7mozjf4fcyrexwbe6a.gz). Longwave and shortwave lookup tables are stored in `cloudysky_lw.nc` and `cloudysky_sw.nc`, respectively. `LookUpCld` reads the lookup tables. The tabulated information includes tabulated data of extinction coefficients, single scattering albedo, and asymmetry parameter for liquid and ice particles. The optical properties for liquid particles are computed at radius ``2.5 \mathrm{\mu m} ≤ r_l ≤ 21.5 \mathrm{\mu m}`` in increments of ``1 \mathrm{\mu m}``. The optical properties for ice particles are computed at radius ``10 \mathrm{\mu m} ≤ r_i ≤ 180 \mathrm{\mu m}`` in increments of ``10 \mathrm{\mu m}``.

### Cloud overlap method
Climate models predict grid-mean cloud fraction and cloud condensates, and the vertical structure of clouds needs to be prescribed using some overlap assumptions for calculating radiative fluxes. The Monte Carlo independent column approximation (McICA) allows us to have fractional cloud areas, by randomly assigning some wavelengths to see cloud and other wavelengths to see no cloud. `build_cloud_mask!` builds McICA-sampled cloud masks from cloud fraction. We use the maximum-random overlap method, which maximizes the cloud overlap between adjacent layers but randomly distributes clouds at different altitudes separated by clear sky. The random numbers are generated for each layer and each g-point. For a certain layer (layer 1), the g-points that have been assigned a random number larger than one minus the cloud fraction see cloud, and the other g-points see clear-sky. The random numbers are recalculated for the layer below (layer 2). If the g-point sees cloud in layer 1, its random number in layer 2 is changed to the one in layer 1. If the g-point sees clear-sky in layer 1, its random number in layer 2 is multiplied by one minus the cloud fraction of layer 1. This ensures that the random numbers in layer 2 are randomly distributed in the range [0, 1].