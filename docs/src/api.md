# API

RRTMGP has an API for creating various types of solvers, and accessing data
passed to it. Hosts read and write a solver's data through named getters; their
uniform layout/domain-masking/writability contract — the RRTMGP–ClimaCore
decoupling mechanism — is spelled out under [The getter contract](@ref).

```@meta
CurrentModule = RRTMGP
```

## Radiation modes

```@docs
RRTMGP.AbstractRRTMGPMethod
RRTMGP.AllSkyRadiation
RRTMGP.AllSkyRadiationWithClearSkyDiagnostics
RRTMGP.GrayRadiation
RRTMGP.ClearSkyRadiation
```

## Grid parameters

```@docs
RRTMGP.RRTMGPGridParams
```

## RRTMGPSolver

```@docs
RRTMGP.RRTMGPSolver
RRTMGP.radiation_method
RRTMGP.optical_thickness_parameter
```

## Lookup tables

```@docs
RRTMGP.lookup_tables
```

## Computing fluxes

```@docs
RRTMGP.update_fluxes!
RRTMGP.update_sw_fluxes!
RRTMGP.update_lw_fluxes!
RRTMGP.update_net_fluxes!
```

## Input validation

```@docs
RRTMGP.check_values
RRTMGP.validate_inputs
```

## Numerical policy

```@docs
RRTMGP.Numerics
RRTMGP.Numerics.k_min
RRTMGP.Numerics.τ_thresh
RRTMGP.Numerics.resonance_window
RRTMGP.Numerics.μ₀_min
RRTMGP.Numerics.pow_fast
```

## Spectrally-resolved fluxes

Optional per-band fluxes, enabled with `spectral_fluxes = true` when constructing the
[`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver). The [`spectral_lw_flux_up`](@ref
RRTMGP.spectral_lw_flux_up) docstring covers the full `spectral_{lw,sw}_flux_{up,dn,net}`
family.

```@docs
RRTMGP.spectral_lw_flux_up
RRTMGP.lw_band_bounds
RRTMGP.sw_band_bounds
RRTMGP.Fluxes.FluxBand
```

## Grid adaptation

```@docs
RRTMGP.AbstractInterpolation
RRTMGP.AbstractBottomExtrapolation
RRTMGP.interpolate_levels!
RRTMGP.add_isothermal_boundary_layer!
RRTMGP.clip!
RRTMGP.update_concentrations!
RRTMGP.get_p_min
```

## Aerosol properties

```@docs
RRTMGP.aerosol_radius
RRTMGP.aerosol_column_mass_density
RRTMGP.aerosol_index_map
RRTMGP.canonical_aerosol_name
RRTMGP.aerosol_index
```

## Volume mixing ratios

```@docs
RRTMGP.volume_mixing_ratio
VolumeMixingRatios.VolumeMixingRatioGlobalMean
```

## Standalone API

```@docs
RRTMGP.default_parameters
RRTMGP.solve_gray
RRTMGP.heating_rate
```

## Helpers

```@docs
RRTMGP.gas_names_sw
RRTMGP.aerosol_names
```
