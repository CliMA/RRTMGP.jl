# How to get per-band (spectral) fluxes

By default, RRTMGP retains only broadband fluxes. For applications that need
the spectral decomposition, the solver can retain per-band fluxes. Examples
include atmospheric chemistry, where photolysis rates depend on the
ultraviolet and visible fluxes; land models, where the radiation a canopy
absorbs divides into photosynthetically active (visible) and near-infrared
bands with very different vegetation albedos; band-by-band forcing
diagnostics; and satellite-channel proxies.

## Request per-band fluxes at construction

```julia
solver = RRTMGP.RRTMGPSolver(
    grid_params, method, params, bcs_lw, bcs_sw, as;
    lookups,
    spectral_fluxes = true,
)
RRTMGP.update_fluxes!(solver)
```

This is supported for the spectral (non-gray) methods with two-stream optics
in both bands; other configurations raise an informative error at
construction. The per-band buffers add `(nlev, ncol, n_bnd)` arrays per band
set, so they are opt-in.

## Read them

```julia
F_lw = RRTMGP.spectral_lw_flux_up(solver)   # (nlev, ncol, nbnd_lw) view
F_sw = RRTMGP.spectral_sw_flux_dn(solver)   # (nlev, ncol, nbnd_sw) view
```

Band `b`'s slice `F_lw[:, :, b]` has the same layout as the broadband getters,
and the bands sum to the broadband flux:

```julia
sum(RRTMGP.spectral_lw_flux_up(solver); dims = 3) ≈ RRTMGP.lw_flux_up(solver)
```

The `spectral_*_flux_net` getters are views into a retained per-band net-flux
buffer, updated on every solve like the `up`/`dn` buffers.

## Identify the bands

[`lw_band_bounds`](@ref RRTMGP.lw_band_bounds) and
[`sw_band_bounds`](@ref RRTMGP.sw_band_bounds) return the `(2, n_bnd)`
wavenumber edges (cm⁻¹) of each band:

```julia
wn = RRTMGP.lw_band_bounds(solver)
wn[:, 1]   # lower/upper wavenumber of longwave band 1
```

The RRTMGP longwave tables have 16 bands and the shortwave tables 14,
following [pincus2019](@citet).
