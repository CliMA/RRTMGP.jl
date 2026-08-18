# API

This page lists the public API: the radiation methods, the solver and its
constructors, the flux update, and the named getters through which hosts read
and write a solver's data. The getters' layout, domain-masking, and writability
contract (the mechanism that decouples RRTMGP from ClimaCore) is spelled out
under [The getter contract](@ref).

## Versioning and API stability

`RRTMGP.PUBLIC_NAMES` lists the public API. On Julia 1.11 and later, the same
list is declared with `public`, so `Base.ispublic` agrees with it. Names on the
list change only with the package version; anything else reachable through the
module may change in a patch release.

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
RRTMGP.LookupBundle
RRTMGP.save_lookup_tables
RRTMGP.load_lookup_tables
```

## Computing fluxes

```@docs
RRTMGP.update_fluxes!
RRTMGP.prepare_atmosphere!
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

## Spectrally resolved fluxes

Optional per-band fluxes, enabled with `spectral_fluxes = true` when
constructing the [`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver). Summing a getter
over its band dimension recovers the corresponding broadband flux; see
[Get per-band (spectral) fluxes](howto/spectral_fluxes.md).

```@docs
RRTMGP.spectral_lw_flux_up
RRTMGP.spectral_lw_flux_dn
RRTMGP.spectral_lw_flux_net
RRTMGP.spectral_sw_flux_up
RRTMGP.spectral_sw_flux_dn
RRTMGP.spectral_sw_flux_net
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
RRTMGP.get_t_min
RRTMGP.get_t_max
```

## Level interpolation schemes

The schemes, and the two functions that define one, are derived and compared on
the [Level interpolation and extrapolation](@ref) page.

```@docs
RRTMGP.NoInterpolation
RRTMGP.ArithmeticMean
RRTMGP.GeometricMean
RRTMGP.UniformZ
RRTMGP.UniformP
RRTMGP.BestFit
RRTMGP.SameAsInterpolation
RRTMGP.UseSurfaceTempAtBottom
RRTMGP.HydrostaticBottom
RRTMGP.requires_z
RRTMGP.interp!
RRTMGP.extrap!
RRTMGP.uniform_z_p
RRTMGP.best_fit_p
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
RRTMGP.set_volume_mixing_ratio!
VolumeMixingRatios.VolumeMixingRatioGlobalMean
```

## Standalone API

```@docs
RRTMGP.default_parameters
RRTMGP.AtmosphereProfile
RRTMGP.standard_atmosphere
RRTMGP.solve
RRTMGP.RadiationOutput
RRTMGP.solve_gray
RRTMGP.heating_rate
```

## Helpers

```@docs
RRTMGP.gas_names_sw
RRTMGP.aerosol_names
```
