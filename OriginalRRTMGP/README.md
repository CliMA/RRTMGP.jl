# Regenerating the Fortran reference data

These scripts rebuild the reference [rte-rrtmgp](https://github.com/earth-system-radiation/rte-rrtmgp)
Fortran implementation and run its example drivers, to regenerate the reference
flux files that RRTMGP.jl's tests compare against (see `test/reference_files.jl`
and the [test problems](../docs/src/Example.md) page). They are a developer
convenience, not part of the package: nothing in `src/` or the test suite calls
them, and CI does not run them.

They assume a sibling checkout of rte-rrtmgp,

```
SomeDirectory/RRTMGP.jl/       # this repository
SomeDirectory/rte-rrtmgp/      # the Fortran reference
```

and a working `gfortran` and NetCDF (C and Fortran) toolchain. Run them from
this directory.

- `ConfigRRTMGP.jl` — sets the NetCDF paths, compiler flags, and rte-rrtmgp
  directory environment variables. Included by the others; not run on its own.
- `CleanRRTMGP.jl` — removes build artifacts and generated output from a prior
  rte-rrtmgp build (`make clean` plus the stray NetCDF outputs and templates).
- `AutoRunRRTMGP.jl` — builds rte-rrtmgp and runs its RFMIP clear-sky and/or
  all-sky example drivers, producing the reference flux files. Toggle the cases
  with the `run_clear_sky` / `run_all_sky` flags at the top.
- `TarNewReferenceData.jl` — gathers every `.nc` file under the rte-rrtmgp tree
  into a `rte-rrtmgp-data` copy and compresses it to `rte-rrtmgp-data.tar.gz`,
  the archive published as the test-data artifact.
