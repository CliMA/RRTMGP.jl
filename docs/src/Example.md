# Example

The quickest way to run RRTMGP is the standalone gray-atmosphere one-liner, which
needs no NetCDF data:

```julia
julia> using RRTMGP

julia> out = RRTMGP.solve_gray(Float64; nlay = 60, ncol = 1);

julia> out.net          # net flux at each level [W/m²], (nlev, ncol)
julia> out.heating_rate # heating rate at each layer [K/s], (nlay, ncol)
```

To assemble the pieces yourself — the atmospheric state, the RTE workspaces, and
`solve_lw!`/`solve_sw!` — see [The functional core](@ref), which walks through the
same gray problem step by step and sketches the clear-sky (lookup-table) path.

Beyond that, the test suite doubles as a set of complete, validated examples.
Each driver below builds its states from reference NetCDF data (downloaded
automatically as artifacts), runs the solvers, and compares against reference
results. Run them from the repository root with the `test` project
(`julia --project=test`).

## Gray radiation

`test/gray_atm.jl` exercises longwave and shortwave gray radiation with both the
non-scattering and two-stream solvers. For longwave-only gray radiation an
analytical radiative-equilibrium solution exists; `gray_atmos_lw_equil` integrates
to equilibrium and compares against it:

```julia
julia> include("test/gray_atm_utils.jl");

julia> gray_atmos_lw_equil(ClimaComms.context(), TwoStreamLWRTE, Float64);
Test Passed
```

Here is the vertical profile of temperature (`T_ex_lev`) in radiative equilibrium:

![](assets/gray_lw_T.png)

`gray_atmos_sw_test` computes shortwave-only gray fluxes and compares to the exact
solution:

```julia
julia> gray_atmos_sw_test(ClimaComms.context(), TwoStreamSWRTE, Float64, 1);
Test Passed
```

Here is the vertical profile of the downward shortwave radiative flux (`flux_dn_dir`):

![](assets/gray_sw_flux_dn.png)

## Gas optics (clear sky)

`test/clear_sky.jl` runs RRTMGP for the RFMIP clear-sky atmosphere states and
compares the results to reference data (running the file drives both longwave
solver types and the two-stream shortwave solver at `ncol = 250`):

```julia
julia> include("test/clear_sky.jl")
```

Here are the vertical profiles of downward longwave (`flux_dn_lw`) and shortwave
(`flux_dn_sw`) fluxes for the first column:

![](assets/clear_sky_lw_flux_dn.png) ![](assets/clear_sky_sw_flux_dn.png)

## Cloud and aerosol optics (all sky)

`test/all_sky_with_aerosols.jl` runs RRTMGP for all-sky atmosphere states with
idealized clouds (uniform condensate and particle size in the troposphere, McICA
cloud sampling) plus MERRA aerosols, and compares the results to reference data:

```julia
julia> include("test/all_sky_with_aerosols.jl")
```

Here are the vertical profiles of downward longwave (`flux_dn_lw`) and shortwave
(`flux_dn_sw`) fluxes for the first column:

![](assets/all_sky_lw_flux_dn.png) ![](assets/all_sky_sw_flux_dn.png)
