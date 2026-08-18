# RRTMGP.jl

RRTMGP.jl is a Julia implementation of the radiative transfer solver RTE
(Radiative Transfer for Energetics) and the RRTMGP (RRTM for General circulation
model applications—Parallel) correlated-``k`` gas optics ([pincus2019](@citet)),
based on the reference Fortran implementation
[rte-rrtmgp](https://github.com/earth-system-radiation/rte-rrtmgp). It computes
longwave and shortwave fluxes and heating rates for clear, cloudy, and
aerosol-laden atmospheres on CPUs and GPUs, and it includes an analytic
gray-radiation mode for idealized studies. It is the radiation scheme of the
[CliMA](https://clima.caltech.edu/) Earth System Model.

## Quick start

Build an idealized atmospheric column, then solve it. Gray radiation uses
analytic formulas and needs no input data:

```@example index
using RRTMGP
profile = RRTMGP.standard_atmosphere(Float64; kind = :tropical)
out = RRTMGP.solve(profile; method = RRTMGP.GrayRadiation())  # gray optics
out.heating_rate[1, 1]                # heating rate at bottom layer [K/s]
```

The full gas optics reads lookup tables, downloaded on first use; loading
NCDatasets activates it:

```@example index
using NCDatasets                                # activates the gas optics
out = RRTMGP.solve(profile)                     # clear-sky correlated-k
out.lw_flux_up[end, 1]                                # OLR [W/m²]
```

## Code structure

RRTMGP is split into two parts:

- **Optics** computes optical properties and source functions given atmospheric
  conditions (pressure, temperature, gas concentrations, clouds, aerosols).
- **RTE** computes radiative fluxes given optical properties and source
  functions.

On top of this functional core sit the [`RRTMGPSolver`](@ref
RRTMGP.RRTMGPSolver), which host models drive through the [getter
contract](getters.md), and standalone entry points (`solve_gray`,
`standard_atmosphere`, `solve`) for single-column work.

## How the documentation is organized

- **Tutorials** (learning-oriented, runnable end to end): [A first radiation
  calculation](tutorials/getting_started.md) and [Radiative-convective
  equilibrium](tutorials/manabe_rce.md) (Manabe's classic climate-sensitivity
  experiment). Start here if you want to learn about radiative transfer and how
  to run this code.
- **How-to guides** (task-oriented recipes): [driving RRTMGP from a host
  model](howto/host_model.md), [adding clouds and
  aerosols](howto/clouds_aerosols.md), [running on GPUs](howto/gpu.md),
  [caching the lookup tables](howto/lookup_cache.md), [per-band
  fluxes](howto/spectral_fluxes.md), and [the validated test
  problems](Example.md). Start here if you are coupling RRTMGP to a climate
  model.
- **Explanation** (the concepts and the math): [the functional
  core](functional_core.md), the [RTE solvers](RTE.md), the [optics](Optics.md),
  [level interpolation](interpolation.md),
  [Float32 and Float64](precision.md), and the [Fortran and paper
  concordance](concordance.md) for readers coming from rte-rrtmgp.
- **Reference**: the [API](api.md), the [getter contract](getters.md), and
  per-module reference pages.

## Authors

`RRTMGP.jl` is developed by the [Climate Modeling
Alliance](https://clima.caltech.edu/). It is based on RTE+RRTMGP and its
reference [Fortran
implementation](https://github.com/earth-system-radiation/rte-rrtmgp), developed
by Robert Pincus, Eli Mlawer, and Jennifer Delamere ([pincus2019](@citet)).
