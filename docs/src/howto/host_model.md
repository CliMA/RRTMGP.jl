# How to drive RRTMGP from a host model

This guide shows how to embed RRTMGP.jl in a host model (a GCM, a single-column
model, or any code with its own time loop): construct the solver once, then
exchange data through the [getter contract](../getters.md) and call
[`update_fluxes!`](@ref RRTMGP.update_fluxes!) every radiation step. The
[Functional core](../functional_core.md) page explains what these pieces are;
this page shows the order in which to assemble them.

## 1. Construct the grid parameters

[`RRTMGPGridParams`](@ref RRTMGP.RRTMGPGridParams) fixes the float type, the
compute device (via [ClimaComms](https://github.com/CliMA/ClimaComms.jl)), and
the grid size:

```julia
import ClimaComms
ClimaComms.@import_required_backends
using RRTMGP

context = ClimaComms.context()
grid_params = RRTMGP.RRTMGPGridParams(
    Float64;
    context,
    domain_nlay = 63,                  # your model's layer count
    ncol = ncolumns,
    isothermal_boundary_layer = true,  # add a buffer layer to the model top
)
```

`domain_nlay` is the number of physical layers. With
`isothermal_boundary_layer = true`, RRTMGP internally adds one extra layer
between your model top and the pressure minimum of the gas-optics tables, and
every getter masks it off (you read and write arrays sized to your domain).

## 2. Choose a radiation method and build the lookup tables

The method selects the physics and which lookup tables are loaded (load
NCDatasets first; the spectral methods live in a package extension):

```julia
using NCDatasets

method = RRTMGP.AllSkyRadiationWithClearSkyDiagnostics(
    true,  # aerosol_radiation
    false, # reset_rng_seed (McICA reproducibility)
)
lookups = RRTMGP.lookup_tables(grid_params, method)
```

Build the tables once and pass them to every solver you construct; see [How to
cache the lookup tables](lookup_cache.md) to skip the NetCDF read across
sessions.

## 3. Build the state and boundary conditions

The host owns the [`AtmosphericState`](@ref
RRTMGP.AtmosphericStates.AtmosphericState) (pressures, temperatures, mixing
ratios, optional [`CloudState`](@ref RRTMGP.AtmosphericStates.CloudState) and
[`AerosolState`](@ref RRTMGP.AtmosphericStates.AerosolState); see
[How to add clouds and aerosols](clouds_aerosols.md)) and the boundary
conditions ([`LwBCs`](@ref RRTMGP.BCs.LwBCs), [`SwBCs`](@ref RRTMGP.BCs.SwBCs)).
All arrays live on the device described by `grid_params`, with the vertical
dimension first and columns last. `test/read_all_sky_with_aerosols.jl` is a
worked example of filling these from data.

The boundary conditions close the column at both ends. At the top, the
insolation enters through the incident solar flux `toa_flux` [W/m²] and the
cosine of the solar zenith angle `cos_zenith`, both `(ncol,)`; the host computes
them from its solar and orbital code. The CliMA package
[Insolation.jl](https://github.com/CliMA/Insolation.jl) provides both as
functions of latitude, longitude, and time, from the orbital parameters. At the
bottom, the longwave surface emissivity `sfc_emis` and the shortwave direct and
diffuse albedos are `(nbnd, ncol)`, i.e., **per band**, so the surface can be
spectrally selective: vegetation, for example, absorbs strongly in the visible
and reflects in the near-infrared, and snow does the reverse. For a spectrally
flat surface, fill the band dimension with one value:

```julia
FT = Float64
DA = ClimaComms.array_type(ClimaComms.device(context))
(; nbnd_lw, nbnd_sw) = lookups

sfc_emis = fill!(DA{FT}(undef, nbnd_lw, ncolumns), FT(0.98))
bcs_lw = RRTMGP.BCs.LwBCs(sfc_emis, nothing)          # nothing: no prescribed
                                                      # incident longwave flux

cos_zenith  = DA{FT}(undef, ncolumns)                 # filled each step (below)
toa_flux    = DA{FT}(undef, ncolumns)
alb_direct  = fill!(DA{FT}(undef, nbnd_sw, ncolumns), FT(0.06))
alb_diffuse = copy(alb_direct)
bcs_sw = RRTMGP.BCs.SwBCs(cos_zenith, toa_flux, alb_direct, nothing,
                          alb_diffuse)
```

The two `nothing` slots are prescribed incident fluxes at the top (longwave and
diffuse shortwave), used for offline intercomparison cases; a host model
normally leaves them off, so the top boundary sees only the direct solar beam.
After construction, the host updates all of these through the getters
(`cos_zenith(solver)`, `toa_flux(solver)`, `surface_emissivity(solver)`,
`direct_sw_surface_albedo(solver)`, `diffuse_sw_surface_albedo(solver)`) as in
step 5, so a time-varying sun or a changing surface needs no reconstruction.

## 4. Construct the solver

```julia
solver = RRTMGP.RRTMGPSolver(
    grid_params, method, params, bcs_lw, bcs_sw, as;
    lookups,
    interpolation = RRTMGP.ArithmeticMean(),  # fill level values from layers
    # center_z, face_z,              # required for BestFit/HydrostaticBottom
    # deep_atmosphere_inverse_scaling,        # (nlev, ncol) metric factor
)
```

If your model provides level (face) pressures and temperatures directly, keep
the default `interpolation = NoInterpolation()`. Otherwise pick an
[`AbstractInterpolation`](@ref RRTMGP.AbstractInterpolation) and
[`AbstractBottomExtrapolation`](@ref RRTMGP.AbstractBottomExtrapolation);
z-based schemes need `center_z`/`face_z` at construction.

## 5. The radiation step

Write the current model state through the getters (they are writable views into
solver-owned device memory), then solve:

```julia
RRTMGP.layer_temperature(solver) .= T_from_host
RRTMGP.surface_temperature(solver) .= T_sfc
RRTMGP.volume_mixing_ratio(solver, "h2o") .= vmr_h2o
RRTMGP.cos_zenith(solver) .= μ₀

RRTMGP.update_fluxes!(solver)   # allocation-free; returns nothing

F_net = RRTMGP.net_flux(solver) # (nlev, ncol) domain-masked view
# heating_rate allocates; prefer the net-flux divergence in inner loops
hr    = RRTMGP.heating_rate(solver)
```

`update_fluxes!` runs the full cascade: level interpolation, boundary-layer
fill, clipping to the lookup tables' valid range, column amounts, longwave and
shortwave solves, and the combined net flux. To inspect the prepared state
without solving, call [`prepare_atmosphere!`](@ref RRTMGP.prepare_atmosphere!).

Two host responsibilities to remember (details in the [getter
contract](../getters.md)):

- **Relative humidity** feeds RH-dependent aerosol optics; update
  `layer_relative_humidity(solver)` (or call `compute_relative_humidity!`) if
  needed.
- **Reading well-mixed gases** via `volume_mixing_ratio(solver, "co2")` returns
  a host scalar; do this during setup.

## 6. Per-band (spectral) fluxes

The getters above return broadband fluxes. Some components need the spectral
decomposition: an atmospheric-chemistry scheme whose photolysis rates depend on
the ultraviolet and visible fluxes, or a land model that partitions the
radiation absorbed by a canopy into photosynthetically active (visible) and
near-infrared bands. For these, the solver can retain per-band fluxes. Request
them at construction and read them like the broadband getters, with the band as
a third dimension:

```julia
solver = RRTMGP.RRTMGPSolver(
    grid_params, method, params, bcs_lw, bcs_sw, as;
    lookups,
    spectral_fluxes = true,
)
RRTMGP.update_fluxes!(solver)

F_sw = RRTMGP.spectral_sw_flux_dn(solver)    # (nlev, ncol, nbnd_sw) view
wn = RRTMGP.sw_band_bounds(solver)   # (2, nbnd_sw) wavenumber edges [cm⁻¹]
```

The band edges from [`sw_band_bounds`](@ref RRTMGP.sw_band_bounds) /
[`lw_band_bounds`](@ref RRTMGP.lw_band_bounds) identify each band's spectral
interval, so the host can sum the surface downwelling flux over the bands its
chemistry or canopy scheme needs. See [How to get per-band (spectral)
fluxes](spectral_fluxes.md) for the full recipe and the supported
configurations.

## 7. Debugging: opt-in input validation

During development, enable range checks on every solve:

```julia
RRTMGP.check_values[] = true   # validate_inputs runs inside update_fluxes!
```

Out-of-range inputs (negative pressures, `cos_zenith` outside [-1, 1], albedos
outside [0, 1], negative mixing ratios, NaNs) raise an error naming the
offending field. The checks reduce over device arrays, so keep the toggle off in
production to maintain peak performance.
