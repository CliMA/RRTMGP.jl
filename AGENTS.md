# RRTMGP.jl Agent Guide

## Ecosystem Guidelines

Please refer to the shared CliMA agent index for ecosystem-wide rules regarding architecture, performance, code quality, infrastructure, and workflows:

- [docs/dev-guides/AGENTS.md](docs/dev-guides/AGENTS.md) — Shared CliMA agent guidelines.

> Shared guides live at `docs/dev-guides/` and are vendored from the canonical source:
> https://github.com/CliMA/DeveloperGuides. Edit shared guides there, not here.

## Repo-Specific Guidelines

RRTMGP.jl computes radiative fluxes and heating rates with the RTE + RRTMGP correlated-`k` method, on CPU and GPU and in either `Float32` or `Float64`. Before working here, read the docs that define its public contract and internals:

- [docs/src/getters.md](docs/src/getters.md) — the getter contract: the ClimaCore-free data-exchange interface hosts use to read and write the solver's buffers.
- [docs/src/functional_core.md](docs/src/functional_core.md) — the functional core (`solve_lw!`/`solve_sw!`) beneath the `RRTMGPSolver` + `update_fluxes!` API.
- [docs/src/precision.md](docs/src/precision.md) — the `Float32` numerics: cancellation-avoiding reformulations, `eps`-scaled guard constants (`RRTMGP.Numerics`), and the consistency harness.
- [docs/src/RTE.md](docs/src/RTE.md) and [docs/src/Optics.md](docs/src/Optics.md) — the radiative transfer equations and the optics.

## Local norms

- The public API is layered: core containers → `RRTMGPSolver` + getters + `update_fluxes!` → the standalone entry points (`solve`, `solve_gray`, `standard_atmosphere`). Keep this separation; hosts should exchange data only through the getters.
- `update_fluxes!` must stay allocation-free and type-stable on both devices (asserted in [test/standalone.jl](test/standalone.jl) and [test/all_sky_with_aerosols_utils.jl](test/all_sky_with_aerosols_utils.jl)). Read fluxes through the getters — do not `copy(getter(solver))` or `parent(getter(solver))`; both defeat the copy-free contract.
- Kernels are device-agnostic and generic in the float type: they index buffers `[gcol, ilev]` (`(ncol, nlev)`), while physical storage is column-first on the GPU and vertical-first on the CPU. Numerical guard constants belong in `RRTMGP.Numerics`, each scaled by `eps(FT)`.
- Run `julia --project=.dev .dev/climaformat.jl .` before committing (margin 80). This matches the CI formatter (pinned to Julia 1.10); a plain `format(".")` from a default JuliaFormatter install produces a different diff.
- Prefer `Pkg.test()` over manually `include`ing `test/runtests.jl`: test-only deps (e.g. `JET`) are loaded through the package test path.
- For a local docs build, `Pkg.develop(path = pwd())` in the `docs` environment first (the committed `docs/Manifest.toml` can be stale); restore it with `git checkout docs/Manifest.toml` before committing, never mid-build.
- Match existing style: explicit names, narrow imports, comments that explain why. Follow the CliMA documentation policy (explicit `# Fields` sections in struct docstrings, not `DocStringExtensions.FIELDS`).

## Self-correction

- If a design doc referenced above is discovered to be stale relative to the code, update it.
- If the user gives a correction about how work should be done in this repo, add it to `Local norms` (or another clearly labeled persistent section) so future sessions inherit it.
