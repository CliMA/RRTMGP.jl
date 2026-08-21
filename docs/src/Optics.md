# RRTMGP optics

`compute_optical_props!` computes the optical properties and source functions
given an atmospheric state. For the two-stream approximation, the optical
properties include the optical thickness of the atmosphere layer (``\tau``),
single-scattering albedo (``\omega_0``), and asymmetry parameter (``g``).

## Gray atmosphere optics

The gray atmosphere (more precisely, semi-gray atmosphere) approximates the
radiative properties of the atmosphere by taking the optical properties to be
independent of wavelength, separately in the longwave and shortwave bands.

### Longwave

Two options are currently supported for computing the optical depth for longwave
radiation. The optical depth (``d``) for longwave radiation, with the
"GrayOpticalThicknessSchneider2004" option, follows [schneider2004](@cite):

```math
\begin{align}
d(\phi, p) = d_0(\phi) \left(\frac{p}{p_0}\right)^\alpha .
\end{align}
```

Here, ``\alpha`` is ratio of the pressure scale height to the partial-pressure
scale height of the infrared absorber (``\alpha=3.5`` for water vapor); ``\phi``
is the latitude; ``p`` is the pressure, and ``p_0`` is the surface pressure. The
optical depth at ``\phi`` and ``p_0``, ``d_0``, is calculated from the
latitude-dependent radiative equilibrium temperature ``T_s(\phi)``:

```math
\begin{align}
T_s(\phi) = T_e + \Delta T \left(\frac{1}{3} - \sin^2\phi\right), \\
d_0(\phi) = \left(\frac{T_s(\phi)}{T_t}\right)^4 - 1.
\end{align}
```

Here, ``T_e`` is the global-mean surface temperature in radiative equilibrium,
``T_t`` is the temperature at the top of the atmosphere, and ``\Delta T`` is the
equator-to-pole temperature difference in radiative equilibrium. The default
values are ``T_e = 300~\mathrm{K}``, ``T_t = 200~\mathrm{K}``, and
``\Delta T = 60~\mathrm{K}``.

The optical depth (``d``) for longwave radiation, with the
"GrayOpticalThicknessOGorman2008" option, follows [frierson2006](@cite) and
[ogorman2008](@cite):

```math
\begin{align}
d(\phi, p) = \alpha \left[f_l \sigma + (1 - f_l) \sigma^4 \right] \left[ \tau_e + (\tau_p - \tau_e) \sin^2\phi \right],
\end{align}
```

where ``f_l = 0.2``, ``\sigma = p / p_0`` is pressure ``p`` normalized by
surface pressure ``p_0``, ``\phi`` is latitude, and the longwave optical
thicknesses at the equator and at the pole are ``\tau_e = 7.2`` and
``\tau_p = 1.8``, respectively. The prefactor ``\alpha`` (default 1) rescales
the total longwave optical depth; [ogorman2008](@citet) vary it to simulate a
range of climates. This ``\alpha`` is unrelated to the exponent ``\alpha`` of
the Schneider option above; each matches the ``α`` field of its own parameter
struct.

The optical thickness of an atmosphere layer (the differential optical depth) of
pressure thickness ``\Delta p`` is

```math
\begin{align}
\tau(\phi, p) = \alpha \frac{\Delta p}{p} \left[f_l \sigma + 4 (1 - f_l) \sigma^4 \right] \left[ \tau_e + (\tau_p - \tau_e) \sin^2\phi \right].
\end{align}
```

The source function for longwave radiation is calculated as
``S = \sigma T^4 / \pi``, where ``T`` is the air temperature and ``\sigma`` is
the Stefan–Boltzmann constant.

### Shortwave

Two options are currently supported for computing the optical depth for
shortwave radiation. The optical depth (``d``) for shortwave radiation, with the
"GrayOpticalThicknessSchneider2004" option, is 0.

The optical depth (``d``) for shortwave radiation, with the
"GrayOpticalThicknessOGorman2008" option, follows [ogorman2008](@cite):

```math
\begin{align}
d(p) = \tau_0 \left(\frac{p}{p_0}\right)^2 ,
\end{align}
```

where ``p`` is the pressure of the atmosphere layer and ``p_0`` is the surface
pressure. The default value for ``\tau_0`` is 0.22.

The single scattering albedo and asymmetry parameter are zero.

## Gas optics

At any single frequency, gaseous absorption obeys Beer's law, and the transfer
equations of the [RTE page](RTE.md) apply directly. The difficulty is spectral:
across an absorption band, the monochromatic absorption coefficient ``k_\nu``
sweeps over orders of magnitude through thousands of lines, so the band-averaged
transmission over an absorber path ``u``,

```math
\bar{t}(u) = \frac{1}{\Delta\nu} \int_{\Delta\nu} e^{-k_\nu u}\, d\nu,
```

cannot be represented by any single mean absorption coefficient: the exponential
weights weak and strong lines differently at every path length. The
``k``-distribution method [lacis1991](@cite) replaces the integral over
frequency by an integral over the cumulative distribution of absorption
strength: with ``g(k)`` the fraction of the band where the absorption
coefficient is below ``k``, reordering the spectrum by strength turns the ragged
``k_\nu`` into a smooth, monotonic ``k(g)`` on ``g \in [0, 1]``, and

```math
\bar{t}(u) = \int_0^1 e^{-k(g)\,u}\, dg \approx \sum_{j=1}^{G} \omega_j\, e^{-k(g_j)\, u},
```

a quadrature over a handful of *g-points*. The *correlated*-``k`` assumption is
that the frequency-to-``g`` ordering is the same at all pressures and
temperatures along the path, so a single ``g`` coordinate serves the whole
column, and each g-point behaves like one monochromatic radiative-transfer
problem. RRTMGP's tables carry 16 g-points in each of 16 longwave and 14
shortwave bands (256 and 224 in total); the bands are non-overlapping,
contiguous, and span the thermal and solar spectra [pincus2019](@cite). The RTE
is solved once per g-point and column, and fluxes add over g-points; the
quadrature weights ``\omega_j`` enter through the per-g-point source terms
described below.

`compute_optical_props!` computes each g-point's optical properties from the
pressure ``p``, temperature ``T``, and gas volume mixing ratios ``\chi`` of the
atmospheric state. (The volume mixing ratio relates to the mass fraction ``q``
of a gas by ``\chi = q / q_d \cdot M_d / M``, where ``q_d = 1 - q_t`` is the
dry-air mass fraction, ``M_d`` the molar mass of dry air, and ``M`` that of the
gas.) The optical thickness of a layer splits into

```math
\tau_g = \tau_\mathrm{major} + \tau_\mathrm{minor} + \tau_\mathrm{rayleigh},
```

with the Rayleigh term present in the shortwave only.

**Major species.** Within each band, absorption is dominated by at most two
gases. Their relative abundance enters through the binary mixing parameter

```math
\eta = \frac{\chi_1}{\chi_1 + r\,\chi_2},
```

where ``r`` is the ratio of the two abundances in the reference atmosphere used
to build the tables, so ``\eta`` runs from 0 (only species 2) to 1 (only species
1). The absorption coefficient is tabulated on a ``(\eta, \log p, T)`` grid and
interpolated trilinearly (`compute_interp_frac_temp`/`_press`/`_η`);
``\tau_\mathrm{major}`` is the interpolated coefficient times the column amount
of the binary mixture.

**Minor species and Rayleigh scattering.** These contributions are linear in the
respective column amounts, with coefficients tabulated against ``(\eta, T)`` at
a representative pressure. Minor gases carry additional density and scaling-gas
corrections where the spectroscopy requires them (continua); Rayleigh scattering
is proportional to the column amount of air, water vapor included. In the
shortwave, the single scattering albedo of the gas optics is
``\omega_0 = \tau_\mathrm{rayleigh} / \tau_g``, and the asymmetry parameter is
zero (Rayleigh scattering is nearly symmetric fore–aft).

**Source terms.** In the longwave, each g-point emits the fraction
``f_g(\eta, p, T)`` (the *Planck fraction*, interpolated like the absorption
coefficient) of the band-integrated Planck emission ``B_\mathrm{band}(T)``,
which is tabulated against temperature at 1 K resolution; the fractions of a
band sum to one, so summing g-points recovers the band's blackbody emission. In
the shortwave, each g-point carries a fixed fraction of the total solar
irradiance, scaled to the incident flux the host prescribes at the top of the
atmosphere. These per-g-point sources are where the quadrature weights of the
``k``-distribution live, which is why fluxes add over g-points without further
weighting.

### Lookup tables

The lookup tables for the spectral map between ``\nu`` and ``g`` are NetCDF
files from the
[rrtmgp-data](https://github.com/earth-system-radiation/rrtmgp-data) repository
(v1.9), downloaded automatically as Julia artifacts. The longwave and shortwave
lookup tables are `rrtmgp-gas-lw-g256.nc` and `rrtmgp-gas-sw-g224.nc`,
respectively (see `RRTMGP.ArtifactPaths`). `LookUpLW` and `LookUpSW` read the
lookup tables. The tabulated information includes the spectral discretization
(bands), the major and minor gases considered in each band, and tabulated data
of absorption coefficients, Planck fraction, and band-integrated Planck
function. The absorption coefficients and Planck fraction are computed at
pressures ``1\,\mathrm{Pa} \le p \le 109663\,\mathrm{Pa}`` in increments of
``\log(p/1\,\mathrm{Pa}) = 0.2``, temperatures
``160\,\mathrm{K} \le T \le 355\,\mathrm{K}`` in ``15\,\mathrm{K}`` increments,
and ``0 \le \eta \le 1`` in ``1/8`` increments. The band-integrated Planck
function is computed at temperatures
``160\,\mathrm{K} \le T \le 355\,\mathrm{K}`` in ``1\,\mathrm{K}`` increments. A
more detailed description of the spectral structure of the lookup tables can be
found in `Appendix A` in [pincus2019](@cite).

!!! note
    Absorption by both major and minor species is treated separately in the
    upper and lower atmosphere, defined as pressures above and below the
    pressure defined in `press_ref_trop` in the lookup tables (~9948 Pa). This
    corresponds to element 13 in array `press_ref`. There is a single large
    table of absorption coefficients. Table elements 1:13 in the pressure
    dimension correspond to pressures >= `press_ref_trop`, while table elements
    14:60 correspond to pressures <= `press_ref_trop`. Elements 13 and 14 refer
    to the same pressure but for the lower and upper atmosphere, respectively.
    This is why the absorption coefficients lookup table has a pressure
    dimension of 60, but there are only 59 reference pressure levels in
    `press_ref`.

## Cloud optics

`compute_lookup_cld_liq_props` and `compute_lookup_cld_ice_props` compute the
optical properties of liquid and ice cloud particles, whose combination gives
the cloud optical properties. The optical thickness (``\tau``), single
scattering albedo (``\omega_0``), and asymmetry parameter (``g``) for clouds
are:

```math
\begin{align}
\tau_c &= \tau_l + \tau_i, \\
\omega_{0c} &= \frac{\tau_l \omega_{0l} + \tau_i \omega_{0i}}{\tau_l + \tau_i}, \\
g_c &= \frac{\tau_l \omega_{0l} g_l + \tau_i \omega_{0i} g_i}{\tau_l \omega_{0l} + \tau_i \omega_{0i}},
\end{align}
```

where subscripts c, l, and i indicate cloud, liquid, and ice, respectively. The
liquid cloud and ice cloud optical thicknesses are calculated by multiplying the
cloud liquid and ice extinction coefficients (``\mu``) by the cloud liquid (lwp)
and ice water paths (iwp), respectively. That is,
``\tau_l = \mu_l \mathrm{lwp}`` and ``\tau_i = \mu_i \mathrm{iwp}``. The
extinction coefficient, single scattering albedo, and asymmetry parameter are
computed by linearly interpolating the tabulated values in the liquid and ice
particle effective radius. For ice particles, the optical properties also depend
on the surface roughness, following [yang2013](@cite). The optical properties
for shortwave bands are delta-scaled to reduce the biases in calculating
radiative fluxes for highly asymmetric phase functions ([joseph1976](@cite)).

`add_cloud_optics_2stream!` adds cloud optics to the gas optics for each g-point
that sees cloud (the cloud mask is described below).

### Lookup tables

The lookup tables for cloud optics are NetCDF files from the same
[rrtmgp-data](https://github.com/earth-system-radiation/rrtmgp-data) artifacts.
The longwave and shortwave lookup tables are `rrtmgp-clouds-lw-bnd.nc` and
`rrtmgp-clouds-sw-bnd.nc`, respectively. `LookUpCld` reads the lookup tables.
The tabulated information includes tabulated data of extinction coefficients,
single scattering albedo, and asymmetry parameter for liquid and ice particles.
The optical properties for liquid particles are computed at radius
``2.5\,\mathrm{\mu m} \le r_l \le 21.5\,\mathrm{\mu m}`` in increments of
``1\,\mathrm{\mu m}``. The optical properties for ice particles are tabulated
against diameter, ``10\,\mathrm{\mu m} \le d_i \le 180\,\mathrm{\mu m}`` in
increments of ``10\,\mathrm{\mu m}``; `LookUpCld` halves these bounds on load,
so the interpolation and the clamp applied to `cloud_ice_effective_radius` work
in radius, ``5\,\mathrm{\mu m} \le r_i \le 90\,\mathrm{\mu m}``.

### Cloud overlap method

Climate models predict grid-mean cloud fraction and cloud condensate, so the
vertical overlap of clouds must be prescribed to compute radiative fluxes. The
Monte Carlo independent column approximation (McICA) represents fractional
cloudiness by assigning each g-point a random subcolumn that, in each layer,
either sees cloud or sees clear sky; `build_cloud_mask!` builds these masks from
the cloud fraction. We use the maximum-random overlap method: clouds in adjacent
layers overlap maximally, while clouds at different altitudes separated by clear
sky are randomly overlapped.

Proceeding downward for each g-point, layer ``k`` sees cloud where its random
number ``R_k`` exceeds one minus the layer's cloud fraction,
``R_k > 1 - \mathrm{CF}_k``, with

```math
R_k =
\begin{cases}
R_{k-1} & \text{if the g-point sees cloud in layer } k-1,\\
u_k \, (1 - \mathrm{CF}_{k-1}) & \text{otherwise},
\end{cases}
```

where ``u_k \sim U(0, 1)`` is an independent uniform draw; the mixture of the
two cases leaves ``R_k`` uniform on [0, 1]. Three properties follow. Each layer
is cloudy with probability ``\mathrm{CF}_k``, so the sampled mask reproduces the
prescribed cloud fractions in expectation. Vertically contiguous cloudy layers
overlap maximally: a g-point that is cloudy in one layer stays cloudy in the
next wherever the fractions allow, so a contiguous cloudy block has total cloud
cover ``\max_k \mathrm{CF}_k``. Cloudy blocks separated by clear air are
uncorrelated. Because every g-point carries an independent subcolumn, the
broadband flux averages over the samples: the McICA estimate is unbiased, and
its sampling noise largely cancels in the spectral integration. The all-sky
methods' `reset_rng_seed` option reseeds the generator on each
[`update_fluxes!`](@ref RRTMGP.update_fluxes!) call to make the sampling
reproducible (see the reproducibility caveat in [How to run on
GPUs](howto/gpu.md)).

## Aerosol optics

`add_aerosol_optics_1scalar!` and `add_aerosol_optics_2stream!` add aerosol
extinction to the gas (and cloud) optical properties, using the MERRA aerosol
tables shipped with the rrtmgp-data artifacts
(`rrtmgp-aerosols-merra-lw.nc`/`-sw.nc`, read by `LookUpAerosolMerra`). For each
species and band, the tables give the optical properties per unit column mass of
aerosol; the host supplies the per-species column mass densities and particle
sizes (the [`aerosol_column_mass_density`](@ref
RRTMGP.aerosol_column_mass_density) and [`aerosol_radius`](@ref
RRTMGP.aerosol_radius) getters), and layers without aerosol mass are skipped
through a precomputed mask.

The MERRA set comprises 15 species: dust and sea salt in five size bins each,
sulfate, and black and organic carbon each in a hydrophilic and a hydrophobic
variant. Dust and sea salt select their size bin from the supplied particle
radius. The hygroscopic species (sea salt, sulfate, and the hydrophilic carbon
variants) are additionally interpolated in relative humidity, which enters from
the atmospheric state; keeping `layer_relative_humidity` current is a host
responsibility (see [Driving RRTMGP from a host model](howto/host_model.md)).
Dust and the hydrophobic carbon variants are humidity-independent. Internally,
the species occupy fixed positions in the aerosol state arrays: hosts should map
their species through [`aerosol_index_map`](@ref RRTMGP.aerosol_index_map) (or
[`aerosol_names`](@ref RRTMGP.aerosol_names)) rather than assume an ordering,
and the named getters above handle this lookup.

What each solver family receives mirrors the treatment of gas optics. The
non-scattering longwave path adds only the aerosol absorption optical depth
``\tau (1 - \omega_0)``, while the two-stream path folds the full
``(\tau, \omega_0, g)`` of each species into the running optical properties by
the same ``\tau``-weighted mixing rules as for clouds, delta-scaled in the
shortwave. In the shortwave band containing 550 nm, the summed extinction and
scattering aerosol optical depths are retained as diagnostics
(`aod_sw_extinction`, `aod_sw_scattering`), the standard aerosol optical depth.
