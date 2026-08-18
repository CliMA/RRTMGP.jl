# How to add clouds and aerosols

The all-sky radiation methods extend the clear-sky gas optics with cloud and
aerosol optics. This is the production configuration of the CliMA atmosphere.
The steps below assume the setup of
[How to drive RRTMGP from a host model](host_model.md) and cover only what
clouds and aerosols add: the [`CloudState`](@ref
RRTMGP.AtmosphericStates.CloudState) and [`AerosolState`](@ref
RRTMGP.AtmosphericStates.AerosolState), the McICA sampling options, and the
aerosol species ordering. The underlying optics (lookup tables, McICA
cloud-overlap sampling, MERRA aerosol species) are described on the
[Optics](../Optics.md) page.

## 1. Choose an all-sky method

```julia
using RRTMGP, NCDatasets

method = RRTMGP.AllSkyRadiationWithClearSkyDiagnostics(
    true,  # aerosol_radiation
    false, # reset_rng_seed (see step 5)
)
lookups = RRTMGP.lookup_tables(grid_params, method)
```

[`AllSkyRadiation`](@ref RRTMGP.AllSkyRadiation) computes all-sky fluxes;
[`AllSkyRadiationWithClearSkyDiagnostics`](@ref
RRTMGP.AllSkyRadiationWithClearSkyDiagnostics) additionally retains clear-sky
fluxes from the same state (the difference is the cloud radiative effect).
With `aerosol_radiation = true`, [`lookup_tables`](@ref RRTMGP.lookup_tables)
loads the MERRA aerosol tables alongside the gas and cloud tables.

## 2. Build the cloud state

The host owns the [`CloudState`](@ref RRTMGP.AtmosphericStates.CloudState) and
fills its physical inputs: effective radii [µm], water paths [g/m²], and cloud
fraction, all `(nlay, ncol)` on the device. Layers without cloud carry zeros.
The two `Bool` mask arrays are only allocated by the host; RRTMGP fills them
by McICA sampling on every solve. The two optional cloud-cover arrays,
`(ncol,)`, receive the sampled column cloud cover as a diagnostic (readable
through [`lw_cloud_cover`](@ref RRTMGP.lw_cloud_cover) /
[`sw_cloud_cover`](@ref RRTMGP.sw_cloud_cover)); a 9-argument constructor
without them is available if the diagnostic is not needed.

```julia
using RRTMGP.AtmosphericStates: CloudState, MaxRandomOverlap

zeros_lay() = DA(zeros(FT, nlay, ncol))
cloud_state = CloudState(
    zeros_lay(),                  # cld_r_eff_liq [µm]
    zeros_lay(),                  # cld_r_eff_ice [µm]
    zeros_lay(),                  # cld_path_liq [g/m²]
    zeros_lay(),                  # cld_path_ice [g/m²]
    zeros_lay(),                  # cld_frac in [0, 1]
    DA{FT}(undef, ncol),          # cld_cover_sw (diagnostic)
    DA{FT}(undef, ncol),          # cld_cover_lw (diagnostic)
    DA{Bool}(undef, nlay, ncol),  # mask_lw (built each solve)
    DA{Bool}(undef, nlay, ncol),  # mask_sw (built each solve)
    MaxRandomOverlap(),           # cloud-overlap assumption
    2,                            # ice_rgh: 1 none, 2 medium, 3 rough
)
```

The effective radii are clamped to the cloud lookup table ranges (liquid
2.5–21.5 µm, ice 5–90 µm; see [Optics](../Optics.md)), so values outside them
are used at the nearest table edge without warning. Note that the ice tables
are tabulated against diameter in the NetCDF file and halved to radii on
load. `ice_rgh` selects the ice-particle surface roughness of the
[yang2013](@citet) ice optics.

## 3. Build the aerosol state

The host owns the [`AerosolState`](@ref RRTMGP.AtmosphericStates.AerosolState)
and fills, per species and layer, the aerosol column mass `aero_mass` [kg/m²]
and, for dust and sea salt, the particle radius `aero_size` [µm], both
`(n_aero, nlay, ncol)`. The mask is allocated only (RRTMGP rebuilds it from
`aero_mass` on every solve), and the two `(ncol,)` optical-depth arrays are
diagnostics the solver fills (readable through
[`aod_sw_extinction`](@ref RRTMGP.aod_sw_extinction) /
[`aod_sw_scattering`](@ref RRTMGP.aod_sw_scattering)).

```julia
using RRTMGP.AtmosphericStates: AerosolState

n_aero = length(RRTMGP.aerosol_names())   # 15 MERRA species
aerosol_state = AerosolState(
    DA{FT}(undef, ncol),                  # aod_sw_ext (diagnostic)
    DA{FT}(undef, ncol),                  # aod_sw_sca (diagnostic)
    DA{Bool}(undef, nlay, ncol),          # aero_mask (built each solve)
    DA(zeros(FT, n_aero, nlay, ncol)),    # aero_size [µm]
    DA(zeros(FT, n_aero, nlay, ncol)),    # aero_mass [kg/m²]
)
```

!!! warning "Map aerosol species by name, not by position"
    The first index of `aero_mass` and `aero_size` selects the MERRA species,
    in an order fixed by RRTMGP. Filling these arrays positionally, from a
    host species list assumed to be in the same order, produces no error when
    the orders disagree: every species silently receives another species'
    optics. Resolve indices by name instead, through
    [`aerosol_index_map`](@ref RRTMGP.aerosol_index_map) (or
    [`aerosol_names`](@ref RRTMGP.aerosol_names), which lists the species in
    index order), or write through the named getters of step 4. RRTMGP's own
    ordering is locked by `test/aerosol_name_map.jl`; the host-side mapping is
    the host's responsibility.

The 15 species are `dust1` … `dust5` and `sea_salt1` … `sea_salt5` (five size
bins each), `sulfate`, and `black_carbon` / `organic_carbon`, each with a
hygroscopic `_rh` variant. The hygroscopic species (sea salt, sulfate, and the
`_rh` carbon variants) are interpolated in relative humidity, so keep
`layer_relative_humidity(solver)` current (see
[How to drive RRTMGP from a host model](host_model.md)).

## 4. Construct the solver and write the state each step

The cloud and aerosol states are the last two arguments of the
[`AtmosphericState`](@ref RRTMGP.AtmosphericStates.AtmosphericState):

```julia
as = RRTMGP.AtmosphericStates.AtmosphericState(
    lon, lat, layerdata, p_lev, t_lev, t_sfc, vmr,
    cloud_state, aerosol_state,
)
solver = RRTMGP.RRTMGPSolver(
    grid_params, method, params, bcs_lw, bcs_sw, as;
    lookups,
)
```

After construction, update the cloud and aerosol inputs through the getters,
like the rest of the state; the aerosol getters take the species name and
resolve the index internally:

```julia
RRTMGP.cloud_fraction(solver) .= cf                  # (nlay, ncol)
RRTMGP.cloud_liquid_water_path(solver) .= lwp        # [g/m²]
RRTMGP.cloud_liquid_effective_radius(solver) .= r_liq  # [µm]

RRTMGP.aerosol_column_mass_density(solver, "sulfate") .= so4_mass  # [kg/m²]
RRTMGP.aerosol_radius(solver, "dust1") .= dust1_radius             # [µm]

RRTMGP.update_fluxes!(solver)
```

`test/read_all_sky_with_aerosols.jl` is a complete worked example of filling
both states from data.

## 5. McICA sampling and reproducibility

With partial cloud fractions (`0 < cld_frac < 1`), the solvers sample the
cloud masks stochastically on every solve (McICA with maximum-random overlap;
see [Optics](../Optics.md)), drawing from the global random-number generator.
For reproducible fluxes, construct the method with `reset_rng_seed = true` and
pass a seed to the solve:

```julia
RRTMGP.update_fluxes!(solver, seed)   # reseeds the RNG before sampling
```

Without `reset_rng_seed = true`, the seed argument is ignored. On multiple
threads and on the GPU, the sampling of partially cloudy columns is
statistically but not bitwise reproducible; see the caveat in
[How to run on GPUs](gpu.md). Overcast (`cld_frac = 1`) and clear layers are
deterministic.

## 6. Read the diagnostics

With `AllSkyRadiationWithClearSkyDiagnostics`, the clear-sky counterparts of
every flux getter (`clear_lw_flux_up`, `clear_sw_flux_dn`, `clear_net_flux`,
…) hold the fluxes the same state would produce without clouds, so the cloud
radiative effect is, for example,
`RRTMGP.lw_flux_net(solver) .- RRTMGP.clear_lw_flux_net(solver)`. The sampled
column cloud covers and the 550 nm aerosol optical depths are available
through [`lw_cloud_cover`](@ref RRTMGP.lw_cloud_cover) /
[`sw_cloud_cover`](@ref RRTMGP.sw_cloud_cover) and
[`aod_sw_extinction`](@ref RRTMGP.aod_sw_extinction) /
[`aod_sw_scattering`](@ref RRTMGP.aod_sw_scattering).
