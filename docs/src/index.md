# RRTMGP.jl

Julia implementation of [RTE-RRTMGP](https://github.com/earth-system-radiation/rte-rrtmgp).

## Code structure

RRTMGP is fundamentally split into two parts:

 - **Optics** computes optical properties and source functions given atmospheric conditions (e.g., pressure, temperature, gas concentrations).
 - **RTE** computes radiative fluxes given optical properties and source functions.

 ## Authors

`RRTMGP.jl` is being developed by the [Climate Modeling Alliance](https://clima.caltech.edu/).
