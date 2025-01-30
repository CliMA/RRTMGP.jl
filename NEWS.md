RRTMGP.jl Release Notes
========================

main
------

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
