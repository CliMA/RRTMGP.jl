# The getter contract

```@meta
CurrentModule = RRTMGP
```

A [`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver) owns its device buffers (state,
fluxes, scratch). Hosts read and write them through **named getters**. This
getter set forms the data-exchange interface, operating on device arrays in a
documented layout.

## The uniform contract

Every getter that carries a vertical dimension returns a **writable view into
the solver's own buffer**, following one invariant:

- **Layout: physical dimension first, `ncol` last.** Layer-centered quantities
  are `(nlay, ncol)`, level/face quantities are `(nlev, ncol)` (with
  `nlev = nlay + 1`), per-band fluxes are `(nlev, ncol, nbnd)`, and
  surface/column quantities are `(ncol,)` (except the spectral surface
  properties like emissivity and albedos, which are `(nbnd, ncol)`).
- **Domain-masked.** When the solver carries an internal isothermal boundary
  layer (`isothermal_boundary_layer = true`), the getters exclude that extra top
  layer/level, so the returned arrays are sized to the physical domain
  (`domain_nlay` / `domain_nlay + 1`) and line up with the host's grid. The
  getter returns a concrete `view` type, keeping [`update_fluxes!`](@ref
  RRTMGP.update_fluxes!) type-stable.
- **A view, not a copy.** Writing into an input getter mutates the solver's
  buffer in place; reading an output getter reads the buffer. `ncol` is the
  flattened column count (`Nij·Nij·Nh` for a GCM, `1` standalone). Two getters
  have different behavior: `volume_mixing_ratio` for a well-mixed gas (with the
  global-mean `VmrGM` storage) returns a host scalar copied off the device, and
  [`heating_rate`](@ref RRTMGP.heating_rate) computes its result on demand into
  a fresh array. Because the well-mixed scalar is a copy, write it with
  [`set_volume_mixing_ratio!`](@ref RRTMGP.set_volume_mixing_ratio!) rather than
  broadcasting into the getter.

Internal code that needs the full, boundary-extended buffer (rather than the
domain view) uses the raw struct fields (e.g., `solver.net_flux_buffer`) or
`parent(getter(solver))`.

## Writing inputs and reading outputs

Inputs are written *before* the solve, outputs read *after* it:

```julia
# write inputs in place (no copy), then solve
ᶜT = RRTMGP.layer_temperature(solver)          # (nlay, ncol) writable view
ᶜT .= my_temperatures
RRTMGP.update_fluxes!(solver, seed)             # runs the radiation update

F = RRTMGP.net_flux(solver)            # (nlev, ncol) view of the result
```

With ClimaCore, this forms the bridge: `array2field(getter(solver), space)`
wraps a getter's view as a `Field` (for an input, writes flow into the solver's
buffer; for an output, the `Field` reads the buffer in place), and
`field2array(dst) .= getter(solver)` copies an output into an existing
destination field. Both directions are copy-free: the getters are single-level
views of plain `(nlev, ncol)` buffers, so `array2field`'s lazy `reshape` wraps
them without materializing anything. Do **not** `copy(getter(solver))` before
wrapping: the view is directly wrappable, and copying allocates a full
`(nlev, ncol)` field on every call (reading many fluxes per step this way adds
up to large per-step allocations). Likewise avoid `parent(getter(solver))` as a
shortcut: `parent` returns the full boundary-extended buffer, undoing the domain
masking. On the GPU, the returned views are plain `SubArray`s of device arrays;
to bring one to the CPU, use `Array(getter(solver))` (or the element-wise `.=`
copy above). One subtlety: the solver's compute buffers are indexed
column-first, `(ncol, nlev)` (with the physical layout chosen per device for
performance), and `update_fluxes!` ends by copying them into the `(nlev, ncol)`
presentation arrays the flux getters expose; read fluxes after `update_fluxes!`,
not between manual Layer-1 solves.

## Responsibilities: what the host provides vs what RRTMGP derives

The host **writes** the inputs through the getters before calling
`update_fluxes!`:

- layer/level pressure and temperature, gas volume mixing ratios, cloud
  properties (effective radii, water paths, fraction), and aerosol properties;
- surface and solar boundary conditions (`surface_temperature`,
  `surface_emissivity`, the albedos, `cos_zenith`, `toa_flux`);
- for the deep-atmosphere GCM path, `deep_atmosphere_inverse_scaling` (computed
  from the host's geometry).

RRTMGP **derives** the rest inside `update_fluxes!`: the optional preparation
cascade (level values from centers, the isothermal boundary layer, input
clipping, and the dry-air column amount via `compute_col_gas!`), then the optics
and RTE solve, then the net-flux combine.

Two things are **deliberately the host's job**:

- **Relative humidity**: call `compute_relative_humidity!` if your gas optics
  need an up-to-date `layer_relative_humidity`.
- **Idealized profiles** (idealized water vapor or clouds): a host modeling
  choice, applied before the solve.

## Getter reference

Every getter has a docstring for REPL help; the tables below give the
overview, and the contract above supplies the shared rules. All vertical/layer
getters are domain-masked views into solver-owned buffers.

### Layer-center state: inputs, `(nlay, ncol)`

| Getter | Quantity |
|---|---|
| `layer_pressure` | layer-center pressure [Pa] |
| `layer_temperature` | layer-center temperature [K] |
| `layer_relative_humidity` | layer-center relative humidity |

### Level (face) state: inputs, `(nlev, ncol)`

| Getter | Quantity |
|---|---|
| `level_pressure` | level pressure [Pa] |
| `level_temperature` | level temperature [K] |

### Surface and top-of-atmosphere: inputs

| Getter | Quantity | Shape |
|---|---|---|
| `surface_temperature` | surface temperature [K] | `(ncol,)` |
| `surface_emissivity` | longwave surface emissivity | `(nbnd, ncol)` |
| `latitude` | column latitude [degrees] | `(ncol,)`, or `nothing` |
| `cos_zenith` | cosine of the solar zenith angle | `(ncol,)` |
| `toa_flux` | top-of-atmosphere solar flux [W/m²] | `(ncol,)` |
| `direct_sw_surface_albedo`, `diffuse_sw_surface_albedo` | shortwave surface albedo (direct / diffuse) | `(nbnd, ncol)` |
| `top_of_atmosphere_lw_flux_dn` | prescribed incident longwave flux [W/m²] | `(ngpt, ncol)`, or `nothing` |
| `top_of_atmosphere_diffuse_sw_flux_dn` | prescribed incident diffuse shortwave flux [W/m²] | `(ngpt, ncol)`, or `nothing` |

### Fluxes: outputs, `(nlev, ncol)` [W/m²]

| Getter | Quantity |
|---|---|
| `lw_flux_up`, `lw_flux_dn`, `lw_flux_net` | longwave up / down / net |
| `sw_flux_up`, `sw_flux_dn`, `sw_flux_net` | shortwave up / down / net |
| `sw_direct_flux_dn` | shortwave direct-beam downward |
| `net_flux` | combined longwave + shortwave net flux |
| [`heating_rate`](@ref RRTMGP.heating_rate) | radiative heating rate [K/s], `(nlay, ncol)` (computed on demand into a fresh array) |

The clear-sky counterparts (`AllSkyRadiationWithClearSkyDiagnostics` only)
mirror these: `clear_lw_flux_up`/`clear_lw_flux_dn`/`clear_lw_flux`,
`clear_sw_flux_up`/`clear_sw_flux_dn`/`clear_sw_direct_flux_dn`/
`clear_sw_flux`, and `clear_net_flux`. Per-band fluxes
(`spectral_{lw,sw}_flux_{up,dn,net}`, `(nlev, ncol, nbnd)`) are covered on the
[API](@ref) page; all six are views into the retained per-band buffers.

### Clouds and aerosols

| Getter | Quantity | Shape |
|---|---|---|
| `cloud_liquid_effective_radius`, `cloud_ice_effective_radius` | cloud effective radii [µm] | `(nlay, ncol)` |
| `cloud_liquid_water_path`, `cloud_ice_water_path` | cloud water paths [g/m²] | `(nlay, ncol)` |
| `cloud_fraction` | cloud fraction | `(nlay, ncol)` |
| `lw_cloud_cover`, `sw_cloud_cover` | diagnosed column cloud cover | `(ncol,)` |
| [`aerosol_radius`](@ref RRTMGP.aerosol_radius)`(s, name)`, [`aerosol_column_mass_density`](@ref RRTMGP.aerosol_column_mass_density)`(s, name)` | per-species aerosol size / column mass | `(nlay, ncol)` |
| `aod_sw_extinction`, `aod_sw_scattering` | shortwave aerosol optical depth | `(ncol,)` |

### Gases, coordinates, and configuration

| Getter | Quantity |
|---|---|
| [`volume_mixing_ratio`](@ref RRTMGP.volume_mixing_ratio)`(s, name)` | gas volume mixing ratio: a `(nlay, ncol)` view for `"h2o"`/`"o3"`; a host scalar for well-mixed gases with the default global-mean (`VmrGM`) storage. With full per-layer `Vmr` storage, every gas is a `(nlay, ncol)` view. Write with [`set_volume_mixing_ratio!`](@ref RRTMGP.set_volume_mixing_ratio!)`(s, name, value)` |
| `center_z`, `face_z` | layer-center / level altitudes [m], or `nothing` if not provided |
| `deep_atmosphere_inverse_scaling` | deep-atmosphere flux scaling `(nlev, ncol)`, or `nothing` |
| [`radiation_method`](@ref RRTMGP.radiation_method) | the solver's radiation method |
| `isothermal_boundary_layer` | whether the internal boundary layer is present |
| `optical_thickness_parameter` | gray optical-thickness parameters, or `nothing` (non-gray) |

## Docstrings

```@docs
RRTMGP.layer_pressure
RRTMGP.layer_temperature
RRTMGP.layer_relative_humidity
RRTMGP.level_pressure
RRTMGP.level_temperature
RRTMGP.surface_temperature
RRTMGP.surface_emissivity
RRTMGP.latitude
RRTMGP.cos_zenith
RRTMGP.toa_flux
RRTMGP.direct_sw_surface_albedo
RRTMGP.diffuse_sw_surface_albedo
RRTMGP.top_of_atmosphere_lw_flux_dn
RRTMGP.top_of_atmosphere_diffuse_sw_flux_dn
RRTMGP.lw_flux_up
RRTMGP.lw_flux_dn
RRTMGP.lw_flux_net
RRTMGP.sw_flux_up
RRTMGP.sw_flux_dn
RRTMGP.sw_flux_net
RRTMGP.sw_direct_flux_dn
RRTMGP.net_flux
RRTMGP.clear_lw_flux_up
RRTMGP.clear_lw_flux_dn
RRTMGP.clear_lw_flux
RRTMGP.clear_sw_flux_up
RRTMGP.clear_sw_flux_dn
RRTMGP.clear_sw_direct_flux_dn
RRTMGP.clear_sw_flux
RRTMGP.clear_net_flux
RRTMGP.cloud_liquid_effective_radius
RRTMGP.cloud_ice_effective_radius
RRTMGP.cloud_liquid_water_path
RRTMGP.cloud_ice_water_path
RRTMGP.cloud_fraction
RRTMGP.lw_cloud_cover
RRTMGP.sw_cloud_cover
RRTMGP.aod_sw_extinction
RRTMGP.aod_sw_scattering
RRTMGP.center_z
RRTMGP.face_z
RRTMGP.deep_atmosphere_inverse_scaling
RRTMGP.isothermal_boundary_layer
```
