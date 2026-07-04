RRTMGP.jl Release Notes
========================

main
------
- **Breaking:** removed the deprecated positional-argument constructors (the ones that
  warned "use ... with `RRTMGPGridParams`"): the old signatures of `OneScalar`, `TwoStream`,
  `FluxLW`, `FluxSW`, `SourceLWNoScat`, `SourceLW2Str`, `SourceSW2Str`, `NoScatLWRTE`,
  `TwoStreamLWRTE`, `NoScatSWRTE`, `TwoStreamSWRTE`, and `VmrGM`, plus `source_func_longwave`,
  `source_func_shortwave`, and `init_vmr`. Construct these with `RRTMGPGridParams` (and
  `VolumeMixingRatioGlobalMean` for `VmrGM`).
- **Breaking:** removed the internal `RRTMGPData{Order}` data-layout types and helpers
  (`NVCData`/`VCData`/`NCData`/`NData`, `set_cols!`/`set_domain!`/`domain_view`, the index-order
  tags) — they were unused by the getter-based API. This also drops the `LinearAlgebra`
  dependency.
- Removed `@fastmath` from the longwave kernels (kept for shortwave), so the non-scattering
  longwave flux matches CPU/GPU to the tight test tolerance again; the shortwave performance
  gain is retained.
- **Breaking:** the public API is now centered on `RRTMGPSolver` plus named getters
  (`layer_temperature`, `level_pressure`, `net_flux`, ...) built on the functional
  `solve_lw!`/`solve_sw!` core, with the grid-adaptation helpers (`interpolate_levels!`,
  `add_isothermal_boundary_layer!`, `clip!`, `update_concentrations!`) lifted in from
  ClimaAtmos. Standalone entry points `solve_gray`/`default_parameters`/`heating_rate`
  need no NetCDF lookup tables.
- **Breaking:** aerosol names must now be canonical (`sulfate`, `sea_salt1`,
  `black_carbon_rh`, ...); the host-specific short-name aliases (`so4`, `ss1`, ...) and
  `aerosol_aliases()` were removed — hosts map their own names to the canonical ones.
- **Breaking:** the Layer-2 deep-atmosphere getter/kwarg was renamed
  `deep_atmosphere_scaling` → `deep_atmosphere_inverse_scaling` (the caller passes the
  multiplicative inverse; the low-level `metric_scaling` argument is unchanged).
- Added optional spectrally-resolved (per-band) fluxes: construct with
  `spectral_fluxes = true`, then read `spectral_{lw,sw}_flux_{up,dn,net}` and the band
  edges `{lw,sw}_band_bounds`.
- `RRTMGPSolver` and its optics/source workspaces are now `Adapt`-able via custom
  overloads that preserve the packed-buffer view topology across device transfers, so
  checkpoint restores stay zero-allocation.
- **Breaking:** `RRTMGPGridParams` now takes `domain_nlay` (the physical layer count)
  instead of `nlay`, and adds the isothermal boundary layer internally; the stored `nlay`
  field is the total. Callers no longer compute `domain_nlay + 1` themselves.
- **Breaking:** `update_fluxes!`/`update_net_fluxes!` now return `nothing` (previously
  `net_flux(s)`); read results via the getters after the call. This keeps them
  allocation-free and type-stable.
- **Breaking:** the `RRTMGPSolver` output-buffer fields were renamed `net_flux`/
  `clear_net_flux` → `net_flux_buffer`/`clear_net_flux_buffer`, to disambiguate from the
  `net_flux(s)`/`clear_net_flux(s)` getters (the field is the full boundary-extended
  buffer; the getter is the domain-masked view).
- Flux buffers are stored `(nlev, ncol)`. (An intermediate `(ncol, nlev)` layout was tried
  for GPU coalescence, then reverted; getters and fields are `(nlev, ncol)` as before —
  relevant only if you tracked the intermediate branch.)
- The getter data-exchange contract (layout, domain masking, writability, "views into
  solver-owned buffers") is now documented on the "The getter contract" docs page.
- Wrong-mode getters now raise informative errors instead of failing with
  `getproperty(::Nothing, ...)` or a `NamedTuple` field error: the `clear_*` flux getters
  (require `AllSkyRadiationWithClearSkyDiagnostics`), `{lw,sw}_band_bounds` (require a
  spectral method), and the aerosol getters (require `aerosol_radiation = true`).
  `clear_net_flux` previously returned `nothing` on other solvers; it now errors too.
- `volume_mixing_ratio(s, "h2o_self")`/`"h2o_frgn"` with the global-mean (`VmrGM`) storage
  now return the water-vapor field, matching the solver kernels (they share the h2o slot);
  previously they silently returned the unused well-mixed slot (`0`).
- Fixed `Vmrs.VolumeMixingRatioGlobalMean` (threw an `UndefVarError`) and
  `domain_view(grid_params, ::RRTMGPData{VCOrder})`, which skipped the boundary-layer mask
  (the `Bool`-first form was unaffected).
- The `RRTMGPSolver` constructor now checks up front that `center_z`/`face_z` are provided
  when the interpolation/extrapolation scheme requires altitudes (`BestFit`,
  `HydrostaticBottom`), instead of a `MethodError` at solve time.

### Migration

Updating a host adapter (ClimaAtmos/ClimaCoupler) for the changes above:

```diff
  # Grid params: pass the physical domain layer count; RRTMGP adds the boundary layer.
- nlay = domain_nlay + Int(add_isothermal_boundary_layer)
- gp = RRTMGPGridParams(FT; context, nlay, ncol, isothermal_boundary_layer)
+ gp = RRTMGPGridParams(FT; context, domain_nlay, ncol, isothermal_boundary_layer)
+ nlay = gp.nlay        # total, if you still need it (e.g. for state allocation)

  # update_fluxes! no longer returns the flux — read it from a getter afterward.
- F = update_fluxes!(solver, seed)
+ update_fluxes!(solver, seed)
+ F = net_flux(solver)

  # The net_flux struct field was renamed; prefer the getter.
- solver.net_flux            # full (nlev, ncol) buffer
+ net_flux(solver)           # domain-masked view (recommended), or
+ solver.net_flux_buffer     # the raw full buffer

  # Deep-atmosphere scaling kwarg renamed (pass the multiplicative inverse).
- deep_atmosphere_scaling = ms
+ deep_atmosphere_inverse_scaling = ms
```

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

v0.12.0
------
- Started changelog
