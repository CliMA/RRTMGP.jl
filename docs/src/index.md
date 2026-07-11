# RRTMGP.jl

RRTMGP.jl is a Julia implementation of the radiative transfer solver RTE and
the RRTMGP correlated-``k`` gas optics ([pincus2019](@citet)), based on the
reference Fortran implementation
[rte-rrtmgp](https://github.com/earth-system-radiation/rte-rrtmgp). It
computes longwave and shortwave fluxes and heating rates for clear, cloudy,
and aerosol-laden atmospheres, as well as an analytic gray-radiation mode for
idealized studies, on CPUs and GPUs. It is the radiation scheme of the
[CliMA](https://clima.caltech.edu/) Earth System Model.

## Quick start

Build an idealized atmospheric column, then solve it: gray radiation uses
analytic formulas and needs no input data, while the full gas optics loads
NCDatasets to download and parse the lookup tables:

```@example index
using RRTMGP
profile = RRTMGP.standard_atmosphere(Float64; kind = :tropical)
out = RRTMGP.solve(profile; method = RRTMGP.GrayRadiation())  # gray optics
out.heating_rate[1, 1]                          # heating rate at bottom layer [K/s]
```

```@example index
using NCDatasets                                # activates the gas optics
out = RRTMGP.solve(profile)                     # clear-sky correlated-k
out.lw_up[end, 1]                                # OLR [W/m²]
```

## Code structure

RRTMGP is split into two parts:

- **Optics** computes optical properties and source functions given
  atmospheric conditions (pressure, temperature, gas concentrations, clouds,
  aerosols).
- **RTE** computes radiative fluxes given optical properties and source
  functions.

On top of this functional core sits the [`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver),
which host models drive through the [getter contract](getters.md), and a standalone
front door (`solve_gray`, `standard_atmosphere`, `solve`) for single-column work.

## How the documentation is organized

- **Tutorials** (learning-oriented, runnable end to end):
  [A first radiation calculation](tutorials/getting_started.md) and
  [Radiative-convective equilibrium](tutorials/manabe_rce.md) (Manabe's
  classic climate-sensitivity experiment). Start here if you want to learn
  about radiative transfer and how to run this code.
- **How-to guides** (task-oriented recipes):
  [driving RRTMGP from a host model](howto/host_model.md),
  [running on GPUs](howto/gpu.md),
  [caching the lookup tables](howto/lookup_cache.md),
  [per-band fluxes](howto/spectral_fluxes.md), and
  [the validated test problems](Example.md). Start here if you are wiring
  RRTMGP into a climate model.
- **Explanation** (the concepts and the math):
  [the functional core](functional_core.md), the [RTE solvers](RTE.md), the
  [optics](Optics.md), and the
  [Fortran and paper concordance](concordance.md) for readers coming from
  rte-rrtmgp.
- **Reference**: the [API](api.md) and the [getter contract](getters.md).

## Authors

`RRTMGP.jl` is being developed by the
[Climate Modeling Alliance](https://clima.caltech.edu/).
