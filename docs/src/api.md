# API

RRTMGP has an API for creating various types of solvers, and accessing data
passed to it.

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
```

## RRTMGPData

```@docs
RRTMGP.RRTMGPData
RRTMGP.NVCData
RRTMGP.VCData
RRTMGP.NCData
RRTMGP.NData
RRTMGP.set_domain!
RRTMGP.set_cols!
RRTMGP.domain_view
```

## Lookup tables

```@docs
RRTMGP.lookup_tables
```

## Volume Mixing Ratio

```@docs
Vmrs.VolumeMixingRatioGlobalMean
```

## Computing fluxes

```@docs
RRTMGP.update_sw_fluxes!
RRTMGP.update_lw_fluxes!
```

## Aerosol properties

```@docs
RRTMGP.aero_radius
RRTMGP.aero_column_mass_density
```

## Volume mixing ratios

```@docs
RRTMGP.volume_mixing_ratio
```

## Helpers

```@docs
RRTMGP.gas_names_sw
RRTMGP.aerosol_names
```

## Internals

```@docs
RRTMGP.AbstractIndexOrder
RRTMGP.VCOrder
RRTMGP.NOrder
RRTMGP.NCOrder
RRTMGP.NVCOrder
```
