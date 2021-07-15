# RRTMGP.jl

Julia implementation of Rapid and accurate Radiative Transfer Model for General Circulation Models in Parallel (RRTMGP).

## Code structure

RRTMGP is fundamentally split into two parts:

 - Optics** Computes optical depth given atmospheric conditions (pressure, temperature, gas concentrations, grid)
 - **RTE** Computes fluxes given optical depth

