# The functional core

RRTMGP's API is layered, and this page documents its center: a **functional
core** that you can use in a class, in a standalone single-column experiment,
or as the foundation the [`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver) is built
upon. If you want results before concepts, start with the
[tutorials](tutorials/getting_started.md) and return here for what `solve`
assembles internally.

## The model: three inputs and one action

A radiative-transfer solve is split into separately-owned pieces:

1. **Immutable lookup tables**: the correlated-k gas/cloud/aerosol optics data,
   loaded once and shared across columns and timesteps. (Gray radiation uses
   analytic formulas.)
2. **A caller-owned atmospheric state**: [`AtmosphericState`](@ref
   RRTMGP.AtmosphericStates.AtmosphericState) (or [`GrayAtmosphericState`](@ref
   RRTMGP.AtmosphericStates.GrayAtmosphericState)) holding the inputs
   (layer/level pressure and temperature, gas concentrations, cloud and aerosol
   properties) as device arrays laid out `(nlay, ncol)` at layer centers and
   `(nlev, ncol)` at level faces.
3. **A preallocated per-band workspace**: one of [`NoScatLWRTE`](@ref
   RRTMGP.RTE.NoScatLWRTE) / [`TwoStreamLWRTE`](@ref RRTMGP.RTE.TwoStreamLWRTE)
   / [`NoScatSWRTE`](@ref RRTMGP.RTE.NoScatSWRTE) / [`TwoStreamSWRTE`](@ref
   RRTMGP.RTE.TwoStreamSWRTE) which owns the optics scratch and the output
   flux buffers, operating in place.

The action is a single function call: [`solve_lw!`](@ref RRTMGP.RTESolver.solve_lw!) /
[`solve_sw!`](@ref RRTMGP.RTESolver.solve_sw!) reads the state and writes the
workspace's flux buffers in place, and you can read the results:

```julia
solve_lw!(slv_lw, state, lookup_lw)   # or solve_lw!(slv_lw, state) for gray
F = slv_lw.flux.flux_net              # (ncol, nlev) [W/m²]
```

Note the axis order: the raw flux buffers are stored column-first, `(ncol,
nlev)`, so that neighboring GPU threads (one per column) access consecutive
memory. Only the Layer-2 getters (`net_flux`, `lw_flux_up`, ...) present the
vertical-first `(nlev, ncol)` orientation, from presentation copies that
`update_fluxes!` refreshes at the end of each call.

This is what RRTMGP's own tests drive, and what [`RRTMGPSolver`](@ref
RRTMGP.RRTMGPSolver) wraps behind [`update_fluxes!`](@ref RRTMGP.update_fluxes!)
and the named getters.

## Gray radiation

The gray path uses analytic formulas and runs after `using RRTMGP`. It
builds a gray profile, longwave and shortwave two-stream workspaces, solves, and
reads the net fluxes:

```@example functional_core
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

# solve, then read the fluxes out of the workspaces
solve_lw!(slv_lw, as)
solve_sw!(slv_sw, as)

lw_net = slv_lw.flux.flux_net   # (ncol, nlev) [W/m²]
sw_net = slv_sw.flux.flux_net
net = lw_net .+ sw_net
net[1, end]                     # net flux at the top of the atmosphere [W/m²]
```

The one-liner [`solve_gray`](@ref RRTMGP.solve_gray) wraps this.

## Clear-sky and all-sky radiation (with lookup tables)

The spectral path has the same shape as the gray path; the differences are that
the state is built from data and the lookup tables are passed to the solve.
Loading the tables needs `NCDatasets`:

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
as = ...   # e.g., setup_clear_sky_as in test/read_clear_sky.jl

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
for the complete state construction, and the all-sky/aerosol drivers alongside
it for the cloud and aerosol lookups.

## Optional grid adaptation

`solve_*!` assumes a complete, valid state. If a host supplies only
layer-center values, or wants the isothermal boundary layer (an extra buffer
layer above the model top; see [Boundary conditions](RTE.md) on the RTE page),
or must clip unphysical inputs, RRTMGP provides separable, in-place helpers: [`interpolate_levels!`](@ref
RRTMGP.interpolate_levels!), [`add_isothermal_boundary_layer!`](@ref
RRTMGP.add_isothermal_boundary_layer!), [`clip!`](@ref RRTMGP.clip!), and
[`update_concentrations!`](@ref RRTMGP.update_concentrations!), each usable
independently. (Relative humidity is a host responsibility; see
[`update_concentrations!`](@ref RRTMGP.update_concentrations!).)

## From the core to the host aggregate

For a host that runs radiation every step or every ``n`` steps (as is common
in atmosphere models because of the relatively large computational expense of
radiation), [`RRTMGPSolver`](@ref RRTMGP.RRTMGPSolver) bundles the two
workspaces, the state, and the lookups, and drives the pipeline with
[`update_fluxes!`](@ref RRTMGP.update_fluxes!).
Inputs and outputs are exchanged through named getters following the uniform
contract on [The getter contract](@ref) page. See the [API](@ref) reference for
the solver and getter listing.
