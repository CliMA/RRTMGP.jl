# How to cache the lookup tables

The spectral (non-gray) radiation methods need the RRTMGP gas/cloud/aerosol
lookup tables, which are parsed from NetCDF artifacts. Caching the tables avoids
re-reading the NetCDF files within a session and allows running without
NCDatasets in constrained environments.

## Reuse within a session

[`lookup_tables`](@ref RRTMGP.lookup_tables) returns a typed
[`LookupBundle`](@ref RRTMGP.LookupBundle); build it once and pass it to every
solver (and to every [`solve`](@ref RRTMGP.solve) call):

```julia
using RRTMGP, NCDatasets

lookups = RRTMGP.lookup_tables(grid_params, method)
solver_a = RRTMGP.RRTMGPSolver(grid_params, method, params, bcs_lw, bcs_sw,
                               as_a; lookups)
out = RRTMGP.solve(profile; lookups)
```

## Cache across sessions

[`save_lookup_tables`](@ref RRTMGP.save_lookup_tables) serializes a bundle to
disk; [`load_lookup_tables`](@ref RRTMGP.load_lookup_tables) restores it onto
the device described by your grid parameters:

```julia
# once, in an environment with NCDatasets:
using RRTMGP, NCDatasets
lookups = RRTMGP.lookup_tables(grid_params, RRTMGP.ClearSkyRadiation(false))
RRTMGP.save_lookup_tables("clearsky_lookups.jls", lookups)

# later, e.g., in a classroom environment with only RRTMGP:
using RRTMGP
lookups = RRTMGP.load_lookup_tables("clearsky_lookups.jls", grid_params)
out = RRTMGP.solve(profile; lookups)
```

!!! warning "A cache, not an interchange format"
    The file uses Julia's `Serialization` stdlib and is tied to the Julia
    version and package layout that wrote it. Regenerate it when environments
    change; the NetCDF artifacts remain the authoritative data.
