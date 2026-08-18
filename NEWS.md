RRTMGP.jl Release Notes
========================

main
------

v0.22.3
-------
- Documentation cleanup and additions, including new how-to guides for interpolation and clouds/aerosols.
- Added isothermal layer option in standalone runs and RCE calculations to remove noise.
- RCE corrections.
- Added NOTICE file.

v0.22.2
-------
- New **`set_volume_mixing_ratio!(solver, name, value)`**, the write counterpart
  to `volume_mixing_ratio`. It updates a gas's volume mixing ratio in place,
  correctly for **any** gas under either storage backend, without the caller
  knowing the storage class: for spatially varying gases (and every gas under
  per-layer `Vmr` storage) `value` broadcasts over the `(nlay, ncol)` field; for
  a well-mixed gas under the default global-mean (`VmrGM`) storage it writes the
  scalar slot. This fills a gap for well-mixed gases (e.g. prescribed,
  time-varying CO₂): `volume_mixing_ratio` returns those as a **host scalar copy**,
  so `volume_mixing_ratio(solver, "co2") .= value` did not write back. Spatially
  varying gases (`"h2o"`, `"o3"`) were already writable through the getter's view
  and remain so; use whichever reads more clearly (the setter needs no
  storage-class knowledge; the getter view supports in-place read-modify-write).

v0.22.1
-------
- Documented and added a regression test for the **copy-free flux-getter
  contract**: the output flux getters (`net_flux`, `lw_flux_up`, `sw_flux_dn`,
  …) return single-level views of plain `(nlev, ncol)` presentation buffers, so
  a host wraps them with ClimaCore's `array2field` without copying and in the
  correct memory order. Hosts must not `copy(getter(solver))` before wrapping —
  that materializes a full field per getter per call, which for a diagnostics
  pass reading many fluxes each step caused large per-step allocations. No
  behavior change; getters and buffers are unchanged.

v0.22.0
-------
- Internal compute buffers now use device-dependent physical layouts behind a
  uniform column-first `(ncol, nlay/nlev)` indexing convention: optical
  properties (`OneScalar`, `TwoStream`), source functions (`SourceLWNoScat`,
  `SourceLW2Str`, `SourceSW2Str`), and the broadband flux buffers (`FluxLW`,
  `FluxSW`). On the GPU the storage is column-first, giving coalesced access
  (on an A100 at DYAMOND scale this cuts longwave kernel times by 25-36% and
  shortwave all-sky by 14-20%); on the CPU the storage stays vertical-first
  under a lazy `PermutedDimsArray` wrapper, preserving the stride-1 vertical
  sweeps (and CPU performance) while running the same kernel code.
  The Layer-2 flux getters (`lw_flux_up`, `net_flux`, ...) still present
  plain `(nlev, ncol)` views: they read `Fluxes.FluxPresentation` arrays that
  `update_fluxes!` fills from the compute buffers with one fused transposing
  copy per call (so getters remain `Array`-materializable, broadcastable, and
  reducible on the GPU). The opt-in per-band buffers (`FluxBand`) keep the
  host-facing `(nlev, ncol, n_bnd)` layout. **Breaking for Layer-1
  (functional-core) users**: code reading the workspace fields directly
  (`slv.flux.flux_up`, `op.τ`, ...) must swap index order to `[gcol, ilev]`
  or wrap the buffer in `Fluxes.lazy_transpose`. State inputs
  (`AtmosphericState` arrays) and all boundary-condition arrays are
  unchanged.
- New `TransposedStateCache`: column-first copies of the hot
  `AtmosphericState` arrays (layer pressure/temperature/dry-air column
  amount/relative humidity, level temperatures), refreshed by one
  `permutedims!` per spectral solve, so the gas-optics g-point loop reads
  coalesced memory on GPUs. On by default on the GPU (the RTE workspace
  constructors build one; the `RRTMGPSolver` shares a single cache between
  the longwave and shortwave workspaces) and off (`nothing`) on the CPU,
  where the state's vertical-first layout is already cache-friendly; pass
  `state_cache = nothing` to opt out, or pass the same cache to several
  workspaces to share the storage.
- Documentation reorganized along tutorial / how-to / explanation / reference
  lines: two new executable (Literate.jl) tutorials — "A first radiation
  calculation" and "Radiative-convective equilibrium", which reproduces
  Manabe's classic fixed-relative-humidity climate-sensitivity experiment with
  RRTMGP's clear-sky optics — plus four how-to guides (driving RRTMGP from a
  host model, running on GPUs, caching the lookup tables, per-band fluxes).
  The README now gives a Thermodynamics.jl-style overview with quick starts.
- Struct docstrings migrated off `DocStringExtensions.FIELDS` to explicit
  `# Fields` sections (per the CliMA documentation policy), so field
  descriptions and units render on the API pages without inline field
  docstrings.
- New unit tests for previously indirect paths: every interpolation /
  bottom-extrapolation scheme against closed forms (including a dry-adiabat
  round-trip through `interpolate_levels!`), solver-constructor and getter
  guard errors, `validate_inputs` field-by-field, the incident-longwave-flux
  boundary condition, pointwise deep-atmosphere flux scaling, heating-rate =
  net-flux divergence, gray `clip!`, aerosol-only `LookupBundle` round-trips,
  cloud-radiative-effect sign checks in the all-sky tests, and the McICA
  reproducibility tests now also run in `Float32`.
- New advisory GPU benchmark ratchet (`perf/benchmark_ratchet.jl` + a
  soft-fail Buildkite step): the DYAMOND-scale `solve_lw!`/`solve_sw!` medians
  are compared against per-GPU baselines committed under
  `perf/benchmark_baselines/`, flagging wall-time regressions beyond 20%
  without blocking merges. The three DYAMOND benchmark scripts are now
  includable (their sweeps run only when executed as scripts).
- New docs page "Fortran and paper concordance": tables mapping RRTMGP.jl
  names (containers, RTE/gas-optics/cloud/aerosol kernels) to the Fortran
  rte-rrtmgp `mo_*`/`ty_*` names and to the papers whose equations they
  implement, so Fortran-literate readers can navigate instantly.
- The standalone (Layer-3) path is complete: `standard_atmosphere(FT; kind)`
  builds an idealized clear-sky `AtmosphereProfile` (`:tropical`,
  `:midlatitude_summer`, `:subarctic_winter` — analytic two-segment
  temperature, exact hydrostatic pressure, idealized water vapor/ozone,
  present-day well-mixed gases), and `solve(profile; method)` solves it in one
  call — `ClearSkyRadiation` (default) or `GrayRadiation`. Both `solve` and
  `solve_gray` now return a documented `RadiationOutput` struct with stable
  field names (`lw_up` … `sw_direct_dn`, `net`, `heating_rate`, `solver`);
  `solve_gray` previously returned a `NamedTuple` with the same names, so
  field access is unchanged.
- New `prepare_atmosphere!(solver)`: runs the atmospheric-state preparation
  cascade (level interpolation, isothermal boundary layer, clipping, dry-air
  column amounts) without solving, so hosts can inspect the prepared state —
  the prepare/solve split. `update_fluxes!` is now literally
  `prepare_atmosphere!` followed by the three solve/combine steps.
- **Breaking:** `lookup_tables` now returns a typed [`LookupBundle`] instead of
  a nested `(; lookups, lu_kwargs)` `NamedTuple`: stable fields
  (`lookup_lw`, `lookup_sw`, cloud/aerosol tables, the name→index maps, and the
  band/gas counts), `nothing` where a table is absent. Code that indexed the
  NamedTuple (e.g. `bundle.lookups.lookup_lw`, `bundle.lu_kwargs.nbnd_lw`) drops
  one level (`bundle.lookup_lw`, `bundle.nbnd_lw`).
- New `save_lookup_tables(path, bundle)` / `load_lookup_tables(path, grid_params)`
  cache the lookup tables on disk (Julia `Serialization`; a same-version cache,
  not an interchange format), so the spectral methods can run without NCDatasets
  once a cache has been generated — e.g. for standalone/classroom use.
- The CPU and CUDA solver drivers now share device-agnostic per-(g-point,
  column) bodies (`*_gpt_col!` in `src/rte/`), removing the duplicated
  orchestration where backend asymmetries repeatedly crept in. As part of the
  unification, the broadband shortwave direct beam `sw_direct_flux_dn` is now
  accumulated at **every** level: it was previously meaningful only at the
  surface row (a documented limitation inherited from PR #550), with the rows
  above holding the first g-point's profile on CPU and uninitialized memory on
  GPU.
- **Breaking (with aliases):** renamed two modules for clarity — `Vmrs` →
  `VolumeMixingRatios` and `GrayUtils` → `GrayAtmosphere`. The old names remain
  as `const` aliases for one release and will be removed in 0.24.
- **Breaking:** renamed `aerosol_idx()` → `aerosol_index_map()`, resolving the
  one-character collision with `aerosol_index(name)` (which returns a single
  index rather than the whole map).
- Float32 accuracy overhaul of the RTE kernels: series-switch threshold for the
  longwave no-scattering source at eps^(1/4) (was 100·eps), exact `γ1 − γ2`
  identities in the two-stream `k`, expm1-based thin-layer factors, a consistent
  off-resonance evaluation of the shortwave direct reflectance/transmittance,
  exact non-cancelling delta-scaling forms, an addition-built direct-beam
  profile, a continuous gas-optics η interpolation at η = 1, and a restructured
  longwave two-stream source (exact `1 ∓ Rdif − Tdif` factorizations, so no term
  divides by τ and thin layers keep their real O(τ) emission instead of being
  zeroed below a threshold). Measured Float32↔Float64 broadband agreement
  improves 25–90×: every longwave path now sits at its ~2–5e-4 W/m²
  interpolation-noise floor (from 1.1e-2–3.4e-2). A new ratcheting f32↔f64
  consistency test (test/float32_consistency.jl) locks the gains, and the
  Float32 longwave no-scattering reference tolerances tighten from 0.05 to
  5e-3 W/m².
- New `RRTMGP.Numerics` module: every numerical guard constant (`k_min`,
  `τ_thresh`, `resonance_window`, `μ₀_min`) in one place with its derivation.
- New opt-in input validation: set `RRTMGP.check_values[] = true` to have
  `update_fluxes!` validate solver inputs (pressures/temperatures positive and
  finite, `cos_zenith ∈ [-1, 1]`, emissivity/albedos ∈ [0, 1], VMRs ≥ 0) via
  `validate_inputs`; off by default and allocation-free when off.
- Lookup loaders now assert the canonical gas slots the kernels hard-code
  (h2o → 1, o3 → 3) and that `temperature_Planck` is in Kelvin (the g128 file
  variants store an integer index, which would silently corrupt the Planck
  interpolation).
- **Breaking:** the public API is now centered on `RRTMGPSolver` plus a complete set of named
  getters (`layer_temperature`, `level_pressure`, `net_flux`, ...), built on the functional
  `solve_lw!`/`solve_sw!` core. Hosts exchange every input and output through the getters — a
  documented, ClimaCore-free data-exchange contract — and drive the solve with a single
  allocation-free, type-stable `update_fluxes!` that returns `nothing` (results are read back
  through the getters).
- **Breaking:** removed the deprecated positional-argument constructors that warned to use
  `RRTMGPGridParams`: the old signatures of `OneScalar`, `TwoStream`, `FluxLW`, `FluxSW`,
  `SourceLWNoScat`, `SourceLW2Str`, `SourceSW2Str`, `NoScatLWRTE`, `TwoStreamLWRTE`,
  `NoScatSWRTE`, `TwoStreamSWRTE`, and `VmrGM`, plus `source_func_longwave`,
  `source_func_shortwave`, and `init_vmr`. Construct these with `RRTMGPGridParams` (and
  `VolumeMixingRatioGlobalMean` for `VmrGM`).
- **Breaking:** `RRTMGPGridParams` now takes `domain_nlay` (the physical layer count) instead
  of `nlay`, and adds the isothermal boundary layer internally; the stored `nlay` field is the
  total. Callers no longer compute `domain_nlay + 1` themselves.
- **Breaking:** removed the internal `RRTMGPData{Order}` data-layout types and helpers
  (`NVCData`/`VCData`/`NCData`/`NData`, `set_cols!`/`set_domain!`/`domain_view`, and the
  index-order tags); the getter-based API does not use them.
- Lifted the grid-adaptation helpers out of ClimaAtmos into RRTMGP as separable, in-place
  functions over plain `(nlay, ncol)` array views: `interpolate_levels!`,
  `add_isothermal_boundary_layer!`, `clip!`, and `update_concentrations!`.
- `clip!` now also clamps the spectral solvers' layer/level temperatures into the
  valid range of the optics lookup tables — the tables' first and last reference
  temperatures (`lookup.t_ref_min`/`t_ref_max`, 160–355 K for the standard tables),
  read straight from the lookup table — moving the clamp that ClimaAtmos applied
  before every radiation call into RRTMGP's own input preparation. The gray path is
  deliberately not clamped: it uses no lookup tables, and idealized gray atmospheres
  can legitimately reach temperatures outside the lookup range.
- Added standalone entry points that need no NetCDF lookup tables — `solve_gray`,
  `default_parameters`, and `heating_rate` (in K/s) — for single-column/classroom use.
- Added optional spectrally-resolved (per-band) fluxes: construct with `spectral_fluxes = true`,
  then read `spectral_{lw,sw}_flux_{up,dn,net}` and the band edges `{lw,sw}_band_bounds`. Off by
  default, and the broadband path is byte-for-byte unchanged.
- Added a canonical aerosol species name map: the aerosol getters take canonical names
  (`sulfate`, `sea_salt1`, `black_carbon_rh`, ...), and hosts map their own short names
  (`so4`, `ss1`, ...) to these.
- Added `deep_atmosphere_inverse_scaling` (an `RRTMGPSolver` kwarg/getter): a `(nlev, ncol)`
  factor applied to the fluxes for deep-atmosphere geometry, where the host passes the
  multiplicative inverse of its metric scaling (the low-level `metric_scaling` argument is
  unchanged).
- `RRTMGPSolver` and its optics/source workspaces are now `Adapt`-able via custom overloads
  that preserve the packed-buffer view topology across device transfers, so checkpoint restores
  stay zero-allocation.
- Fixed `VolumeMixingRatioGlobalMean`, which threw an `UndefVarError` on every call.
- New documentation: "The functional core" and "The getter contract" pages.

### Migration

`RRTMGPGridParams` now takes the physical `domain_nlay` and adds the isothermal boundary layer
itself, so host adapters (ClimaAtmos/ClimaCoupler) no longer compute the total layer count:

```diff
- nlay = domain_nlay + Int(add_isothermal_boundary_layer)
- gp = RRTMGPGridParams(FT; context, nlay, ncol, isothermal_boundary_layer)
+ gp = RRTMGPGridParams(FT; context, domain_nlay, ncol, isothermal_boundary_layer)
+ nlay = gp.nlay        # total, if you still need it (e.g. for state allocation)
```

Results are read through the getters after the solve: `update_fluxes!(solver, seed)` returns
`nothing`, so use `net_flux(solver)` (a domain-masked view) or the raw `solver.net_flux_buffer`.

v0.21.9
------

- Add cloud cover from McICA
  PR [#599](https://github.com/CliMA/RRTMGP.jl/pull/599)

v0.21.7
------

- Fix edge case cos_zenith = 0 in shortwave solver
  PR [#595](https://github.com/CliMA/RRTMGP.jl/pull/595)

v0.21.0
------

- Aerosol optical depth was added to the `AerosolState`
  PR [#567](https://github.com/CliMA/RRTMGP.jl/pulls/567)

v0.20.1
------

- Clip effective radius by the look-up table range
PR [#568](https://github.com/CliMA/RRTMGP.jl/pull/568)

v0.20.0
------

- Add five size bins of dust and sea salt aerosols.
PR [#564](https://github.com/CliMA/RRTMGP.jl/pull/564)

v0.19.2
-----

- Update cloud optics to the latest version of rrtmgp-data.
PR [#562](https://github.com/CliMA/RRTMGP.jl/pull/562)
- Remove pade approximation. PR [#563](https://github.com/CliMA/RRTMGP.jl/pull/563)

v0.19.1
-----

### Bug fixes

#### Fix  `flux_dn_dir` for non-gray radiation

Prior to this release, `flux_dn_dir` was not correctly set in the two stream
case for non-gray radiation, leading to incorrect values (whatever was in the
memory at initialization). Now, the variable is correctly accumulated over for
every g-point. Note, however, that only the value at the surface (`[1, :]`) is
updated. PR [#550](https://github.com/CliMA/RRTMGP.jl/pull/550).

#### Fix aerosol lookup table

Starting with this release, RRTGMP.jl will use an aerosol look up table that is internally stored, as opposed
to downloading it from the `rrtgmp-data` repository. The reason for this change is that the data distributed
with `rrtgmp-data` contains an an error in the array ordering for the aerosol optics lookup table for sea-salt (‘aero_salt_tbl’).
This error was fixed in the internal table. `RRTGMP.jl` will revert to using `rrtgmp-data`
once the repository updates their tables. PR [#548](https://github.com/CliMA/RRTMGP.jl/pull/548/).
This new lookup table fixes an error in the array ordering for the aerosol optics
lookup table for the shortwave sea-salt data (‘aero_salt_tbl’).
PR [#548](https://github.com/CliMA/RRTMGP.jl/pull/548/).

v0.19.0
-----

- Compute aero_mask internally and store the array.
  ([#528](https://github.com/CliMA/RRTMGP.jl/pull/528))
- Support 1D interpolation on non-uniform grid and fix relative humidity interpolation.
  ([#527](https://github.com/CliMA/RRTMGP.jl/pull/527))

v0.18.0
-----

- Add support for multiple aerosol types ([#523](https://github.com/CliMA/RRTMGP.jl/pull/523))

v0.17.0
-----

- Add support for aerosol optics ([#510](https://github.com/CliMA/RRTMGP.jl/pull/510))

v0.16.0
------

- Fix undefined variable in `rte_sw_noscat_solve!` ([#504](https://github.com/CliMA/RRTMGP.jl/pull/504))
- Add support for OneScalar cloud optics. ([#505](https://github.com/CliMA/RRTMGP.jl/pull/505))
- Rename `rte_lw_noscat!` and simplify the input arguments ([#506](https://github.com/CliMA/RRTMGP.jl/pull/506))

v0.15.1
------

- Force optical thickness to be non-negative ([#502](https://github.com/CliMA/RRTMGP.jl/pull/502))

v0.15.0
------

- Solver struct has been split to allow for independent RTE solver schemes for longwave and shortwave problems ([#492](https://github.com/CliMA/RRTMGP.jl/pull/492))
- Simplify arguments for solve_lw! and solve_sw!. ([#493](https://github.com/CliMA/RRTMGP.jl/pull/493))
- Update Artifacts to use lookup tables and reference data from ([#495](https://github.com/CliMA/RRTMGP.jl/pull/495))
- Move AngularDiscretization to `NoScatLWRTE` ([#496](https://github.com/CliMA/RRTMGP.jl/pull/496))
- Update longwave secants and weights ([#498](https://github.com/CliMA/RRTMGP.jl/pull/498))

v0.14.0
------

- CUDA is now an extension. Some methods have changed.
  ([#485](https://github.com/CliMA/RRTMGP.jl/pull/485)).

v0.13.4
------

Identical to v0.13.2

v0.13.3
------

Broken release. Do not use

v0.13.2
------

v0.13.1
------

- Broadcast FT over Arrays, test with NCDatasets@0.14; update docs env ([#484](https://github.com/CliMA/RRTMGP.jl/pull/484))
- Update argument types for compute_col_gas! ([#470](https://github.com/CliMA/RRTMGP.jl/pull/470))

v0.13.0
------

- Remove inferable fields from AtmosphericStates ([#453](https://github.com/CliMA/RRTMGP.jl/pull/453))
- Add CloudState ([#454](https://github.com/CliMA/RRTMGP.jl/pull/454))
- Restructure datalayout in AtmosphericStates to enable coalesced memory ([#455](https://github.com/CliMA/RRTMGP.jl/pull/455))
- Update from CLIMAParameters to ClimaParams ([#456](https://github.com/CliMA/RRTMGP.jl/pull/456)), Adapt v0.4 ([#462](https://github.com/CliMA/RRTMGP.jl/pull/462))
