# API

RRTMGP has an API for creating various types of solvers, and accessing data
passed to it.

```@meta
CurrentModule = RRTMGP
```

First, users can construct a settings object using `RRTMGPSolver`

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
RRTMGP.radius_dust3
RRTMGP.radius_ss4
RRTMGP.radius_dust4
RRTMGP.radius_dust5
RRTMGP.radius_ss3
RRTMGP.radius_dust2
RRTMGP.radius_ss5
RRTMGP.radius_dust1
RRTMGP.radius_ss1
RRTMGP.radius_ss2
```

```@docs
RRTMGP.column_mass_density_dust5
RRTMGP.column_mass_density_dust2
RRTMGP.column_mass_density_ss5
RRTMGP.column_mass_density_dust3
RRTMGP.column_mass_density_dust4
RRTMGP.column_mass_density_ss3
RRTMGP.column_mass_density_ss1
RRTMGP.column_mass_density_dust1
RRTMGP.column_mass_density_ss2
RRTMGP.column_mass_density_ss4
```

## Volume mixing ratios

```@docs
RRTMGP.volume_mixing_ratio_h2o_frgn
RRTMGP.volume_mixing_ratio_ch4
RRTMGP.volume_mixing_ratio_co2
RRTMGP.volume_mixing_ratio_n2
RRTMGP.volume_mixing_ratio_hfc125
RRTMGP.volume_mixing_ratio_hfc32
RRTMGP.volume_mixing_ratio_o3
RRTMGP.volume_mixing_ratio_cfc11
RRTMGP.volume_mixing_ratio_h2o
RRTMGP.volume_mixing_ratio_hfc23
RRTMGP.volume_mixing_ratio_hfc134a
RRTMGP.volume_mixing_ratio_ccl4
RRTMGP.volume_mixing_ratio_o2
RRTMGP.volume_mixing_ratio_hfc143a
RRTMGP.volume_mixing_ratio_cfc12
RRTMGP.volume_mixing_ratio_no2
RRTMGP.volume_mixing_ratio_co
RRTMGP.volume_mixing_ratio_cf4
RRTMGP.volume_mixing_ratio_n2o
RRTMGP.volume_mixing_ratio_cfc22
RRTMGP.volume_mixing_ratio_h2o_self
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
