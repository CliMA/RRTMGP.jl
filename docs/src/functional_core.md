# The functional core

RRTMGP's API is layered, and this page documents its center: a small **functional core** that
you can drive directly — in a class, in a standalone single-column experiment, or as the
foundation the host-convenience [`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver) is built on. The
core needs no ClimaCore; a single column is simply `ncol = 1`.

## The model: three inputs and one verb

A radiative-transfer solve is split into explicitly-owned pieces:

1. **Immutable lookup tables** — the correlated-k gas/cloud/aerosol optics data, loaded once
   and shared across columns and timesteps. (Gray radiation needs none.)
2. **A caller-owned atmospheric state** — [`AtmosphericState`](@ref
   RRTMGP.AtmosphericStates.AtmosphericState) (or [`GrayAtmosphericState`](@ref
   RRTMGP.AtmosphericStates.GrayAtmosphericState)) holding the inputs — layer/level pressure
   and temperature, gas concentrations, cloud and aerosol properties — as plain device arrays
   laid out `(nlay, ncol)` at layer centers and `(nlev, ncol)` at level faces.
3. **A preallocated per-band workspace** — one of [`NoScatLWRTE`](@ref RRTMGP.RTE.NoScatLWRTE) /
   [`TwoStreamLWRTE`](@ref RRTMGP.RTE.TwoStreamLWRTE) / [`NoScatSWRTE`](@ref
   RRTMGP.RTE.NoScatSWRTE) / [`TwoStreamSWRTE`](@ref RRTMGP.RTE.TwoStreamSWRTE) — which owns the
   optics scratch and the output flux buffers, so a solve allocates nothing.

The verb is [`solve_lw!`](@ref RRTMGP.RTESolver.solve_lw!) / [`solve_sw!`](@ref
RRTMGP.RTESolver.solve_sw!). It reads the state, writes the workspace's flux buffers in place,
and you read the results straight out:

```julia
solve_lw!(slv_lw, state, lookup_lw)   # or solve_lw!(slv_lw, state) for gray
F = slv_lw.flux.flux_net              # (nlev, ncol) [W/m²]
```

This is what RRTMGP's own tests drive, and what [`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver)
wraps behind [`update_fluxes!`](@ref RRTMGP.update_fluxes!) and the named getters.

## Gray radiation (no lookup tables)

The gray path needs no NetCDF data, so it runs in a fresh session after `using RRTMGP`. It
builds an analytic gray profile, longwave and shortwave two-stream workspaces, solves, and
reads the net fluxes:

```julia
using RRTMGP
using RRTMGP.RTE: TwoStreamLWRTE, TwoStreamSWRTE
using RRTMGP.RTESolver: solve_lw!, solve_sw!
using RRTMGP.AtmosphericStates: setup_gray_as_pr_grid, GrayOpticalThicknessOGorman2008
import ClimaComms
ClimaComms.@import_required_backends

FT = Float64
context = ClimaComms.context()
DA = ClimaComms.array_type(ClimaComms.device(context))
nlay, ncol = 60, 1

grid_params = RRTMGP.RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
params = RRTMGP.default_parameters(FT)

# (2) caller-owned state: an analytic gray pressure/temperature profile
lat = DA{FT}([0])   # equator
as = setup_gray_as_pr_grid(context, nlay, lat, FT(1e5), FT(9e3),
                           GrayOpticalThicknessOGorman2008(FT), params, DA)

# (3) preallocated workspaces (own the output flux buffers)
sfc_emis = fill!(DA{FT}(undef, 1, ncol), FT(1))
slv_lw = TwoStreamLWRTE(grid_params; params, sfc_emis, inc_flux = nothing)

cos_zenith = fill!(DA{FT}(undef, ncol), FT(0.5))
toa_flux = fill!(DA{FT}(undef, ncol), FT(1361))
albedo = fill!(DA{FT}(undef, 1, ncol), FT(0.2))
slv_sw = TwoStreamSWRTE(grid_params; cos_zenith, toa_flux,
                        sfc_alb_direct = albedo, inc_flux_diffuse = nothing,
                        sfc_alb_diffuse = copy(albedo))

# the verb: solve, then read the fluxes out of the workspaces
solve_lw!(slv_lw, as)
solve_sw!(slv_sw, as)

lw_net = slv_lw.flux.flux_net   # (nlev, ncol) [W/m²]
sw_net = slv_sw.flux.flux_net
net = lw_net .+ sw_net
```

The one-liner [`solve_gray`](@ref RRTMGP.solve_gray) wraps precisely this.

## Clear-sky and all-sky radiation (with lookup tables)

The spectral path has the same shape; the differences are that the state is built from data
and the lookup tables are passed to the solve. Loading the tables needs `NCDatasets`:

```julia
using NCDatasets, RRTMGP
using RRTMGP.LookUpTables: LookUpLW, LookUpSW
using RRTMGP.RTE: TwoStreamLWRTE, TwoStreamSWRTE
using RRTMGP.RTESolver: solve_lw!, solve_sw!

# load the correlated-k lookup tables from the packaged artifacts
lookup_lw, idx_gases = NCDataset(lw_lookup_file) do ds
    LookUpLW(ds, FT, DA)
end
lookup_sw, _ = NCDataset(sw_lookup_file) do ds
    LookUpSW(ds, FT, DA)
end

# build the AtmosphericState from your data (pressures, temperatures, VMRs, ...)
as = ...   # e.g. setup_clear_sky_as in test/read_clear_sky.jl

grid_params = RRTMGP.RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
slv_lw = TwoStreamLWRTE(grid_params; params, sfc_emis, inc_flux = nothing)
slv_sw = TwoStreamSWRTE(grid_params; cos_zenith, toa_flux,
                        sfc_alb_direct, inc_flux_diffuse = nothing, sfc_alb_diffuse)

# pass the lookup tables to the solve (add cloud/aerosol lookups and a
# metric_scaling argument for cloudy or deep-atmosphere runs)
solve_lw!(slv_lw, as, lookup_lw)
solve_sw!(slv_sw, as, lookup_sw)
lw_net = slv_lw.flux.flux_net
```

See the full clear-sky driver in
[`test/clear_sky_utils.jl`](https://github.com/CliMA/RRTMGP.jl/blob/main/test/clear_sky_utils.jl)
for the complete state construction, and the all-sky/aerosol drivers alongside it for the
cloud and aerosol lookups.

## Optional grid adaptation

`solve_*!` assumes a complete, valid state. If a host supplies only layer-center values, or
wants the isothermal boundary layer, or must clip unphysical inputs, RRTMGP provides
separable, in-place helpers — [`interpolate_levels!`](@ref RRTMGP.interpolate_levels!),
[`add_isothermal_boundary_layer!`](@ref RRTMGP.add_isothermal_boundary_layer!),
[`clip!`](@ref RRTMGP.clip!), and [`update_concentrations!`](@ref
RRTMGP.update_concentrations!) — each usable on its own and none required by `solve_*!`.
(Relative humidity is a host responsibility; see [`update_concentrations!`](@ref
RRTMGP.update_concentrations!).)

## From the core to the host aggregate

For a host that runs radiation every step, [`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver) bundles
the two workspaces, the state, and the lookups, and drives the whole pipeline with one
allocation-free [`update_fluxes!`](@ref RRTMGP.update_fluxes!). Inputs and outputs are then
exchanged through named getters following the uniform contract on [The getter contract](@ref)
page. See the [API](@ref) reference for the full solver and getter listing.
