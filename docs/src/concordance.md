# Fortran and paper concordance

RRTMGP.jl is a Julia implementation of the algorithms in the reference Fortran
package [rte-rrtmgp](https://github.com/earth-system-radiation/rte-rrtmgp)
[pincus2019](@cite). If you know the Fortran code (or the papers it
implements), the tables below map its `mo_*` modules, `ty_*` derived types, and
kernel subroutines to their RRTMGP.jl counterparts, with the papers whose
equations each kernel implements.

Fortran names follow rte-rrtmgp v1.7; where a kernel was merged or renamed
across Fortran versions, both names are given. The shared vocabulary is the
same in both codes: `lay` for layer centers (`nlay`), `lev` for level faces
(`nlev = nlay + 1`), `gpt` for correlated-``k`` g-points, `bnd` for spectral
bands, `sfc` for surface, and the radiative-transfer symbols ``τ`` (optical
depth), ``ω₀``/`ssa` (single-scattering albedo), ``g`` (asymmetry parameter),
and ``μ₀`` (cosine of the solar zenith angle).

## Containers

| RRTMGP.jl | Fortran rte-rrtmgp | Papers / symbols |
|:----------|:-------------------|:-----------------|
| `AtmosphericState` (layer/level ``p``, ``T``; `vmr`; `col_dry`) | the `gas_optics()` arguments `play`, `plev`, `tlay`, `tlev`, `col_dry` | — |
| `GrayAtmosphericState` | (no gray path in the Fortran code) | [schneider2004](@cite), [frierson2006](@cite), [ogorman2008](@cite) |
| `VolumeMixingRatios.Vmr` / `VmrGM` | `mo_gas_concentrations`: `ty_gas_concs` (`set_vmr!`) | — |
| `Optics.OneScalar` (``τ`` only) | `mo_optical_props`: `ty_optical_props_1scl` | — |
| `Optics.TwoStream` (``τ``, ``ω₀``, ``g``) | `mo_optical_props`: `ty_optical_props_2str` | — |
| `SourceLWNoScat` (`lay_source`, `lev_source`, `sfc_source`), `SourceLW2Str` (`lev_source`, `sfc_source`) | `mo_source_functions`: `ty_source_func_lw` (older Fortran versions split `lev_source` into `lev_source_inc`/`dec`) | Planck source ``B`` |
| `SourceSW2Str` | internal arrays of `sw_solver_2stream` | — |
| `FluxLW`, `FluxSW` (`flux_up`, `flux_dn`, `flux_net`, `flux_dn_dir` — same field names) | `mo_fluxes`: `ty_fluxes_broadband` | ``F^↑``, ``F^↓`` |
| `Fluxes.FluxBand` (`spectral_fluxes = true`) | `extensions/mo_fluxes_byband`: `ty_fluxes_byband` | — |
| `LwBCs` (`sfc_emis`, `inc_flux`) | the `rte_lw` arguments `sfc_emis`, `inc_flux` | ``ε_{\mathrm{sfc}}`` |
| `SwBCs` (`cos_zenith`, `toa_flux`, `sfc_alb_direct`, `inc_flux_diffuse`, `sfc_alb_diffuse`) | the `rte_sw` arguments `mu0`, `inc_flux`, `sfc_alb_dir`, `sfc_alb_dif`, `inc_flux_dif` | ``μ₀``, ``α`` |
| `AngularDiscretization` (`gauss_Ds`, `gauss_wts`) | the `gauss_Ds`/`gauss_wts` tables in `mo_rte_solver_kernels` (`n_gauss_angles`) | Gauss–Jacobi-5 quadrature, [hogan2024](@cite) Table 1 |
| `LookupBundle`, `lookup_tables`, `LookUpLW`/`LookUpSW` | `mo_load_coefficients` loading a k-distribution into `mo_gas_optics_rrtmgp`: `ty_gas_optics_rrtmgp` | correlated-``k`` distribution, [pincus2019](@cite) |
| `LookUpCld` | `mo_load_cloud_coefficients` → `mo_cloud_optics_rrtmgp`: `ty_cloud_optics_rrtmgp` (LUT path) | — |
| `LookUpAerosolMerra` | `extensions/mo_aerosol_optics_rrtmgp_merra`: `ty_aerosol_optics_rrtmgp_merra` | MERRA-2 aerosol types |
| `CloudState`, `AerosolState` | the `cloud_optics()` / `aerosol_optics()` arguments (water paths, effective radii, relative humidity) | — |

## RTE solvers and kernels

The Fortran RTE kernels live in `rte-kernels/mo_rte_solver_kernels.F90`; the
RRTMGP.jl equivalents live in `src/rte/`, split by band and solver
(`longwave_noscat.jl`, `longwave_2stream.jl`, `shortwave_noscat.jl`,
`shortwave_2stream.jl`). In RRTMGP.jl the same per-(g-point, column) bodies
(`*_gpt_col!`) are driven by both the CPU and CUDA drivers.

| RRTMGP.jl | Fortran rte-rrtmgp | Papers / symbols |
|:----------|:-------------------|:-----------------|
| `RTESolver.solve_lw!` on a `NoScatLWRTE`/`TwoStreamLWRTE` workspace | `mo_rte_lw`: `rte_lw` | — |
| `RTESolver.solve_sw!` on a `NoScatSWRTE`/`TwoStreamSWRTE` workspace | `mo_rte_sw`: `rte_sw` | — |
| `rte_lw_noscat_solve!` → `lw_noscat_gpt_col!` | `lw_solver_noscat` (`lw_solver_noscat_oneangle`) | — |
| `lw_noscat_source_up` / `lw_noscat_source_dn` | `lw_source_noscat` | linear-in-``τ`` Planck source, [clough1992](@cite) Eq. 13 |
| `rte_lw_noscat_one_angle!` (transport sweep) | `lw_transport_noscat` | ``T = e^{-Dτ}`` at secant ``D`` |
| `rte_lw_2stream_solve!` → `lw_2stream_gpt_col!` | `lw_solver_2stream` | — |
| `lw_2stream_coeffs` | `lw_two_stream` | ``γ₁``, ``γ₂``, ``k = \sqrt{(γ₁-γ₂)(γ₁+γ₂)}``, ``R_{\mathrm{dif}}``/``T_{\mathrm{dif}}`` [meador1980](@cite) with the diffusivity ``D = 1.66`` of [fu1997](@cite) |
| longwave source part of `rte_lw_2stream!` | `lw_source_2str` | linear-in-``τ`` source for scattering atmospheres, [toon1989](@cite) |
| `rte_sw_noscat_solve!` → `sw_noscat_gpt_col!` (`rte_sw_noscat!`) | `sw_solver_noscat` | Beer–Lambert direct beam, ``e^{-τ/μ₀}`` |
| `rte_sw_2stream_solve!` → `sw_2stream_gpt_col!` | `sw_solver_2stream` | — |
| `sw_2stream_coeffs` | `sw_two_stream` + `sw_source_2str` (merged as `sw_dif_and_source` in later versions) | ``γ₁``–``γ₄`` of the Zdunkowski PIFM scheme [zdunkowski1980](@cite); ``R_{\mathrm{dir}}``/``T_{\mathrm{dir}}`` from [meador1980](@cite) Eqs. 14–18 |
| direct-beam prepass inside `rte_sw_2stream!` | the direct-beam accumulation inside `sw_solver_2stream` | ``F^↓_{\mathrm{dir}} = μ₀ F_0 e^{-τ_{\mathrm{cum}}/μ₀}`` |
| adding pass inside `rte_lw_2stream!` / `rte_sw_2stream!` | `adding` | interface albedo recursion, [shonk2008](@cite) |
| `Optics.delta_scale` | `mo_optical_props_kernels`: `delta_scale_2str_kernel` | δ-scaling, [joseph1976](@cite) (RRTMGP.jl uses the algebraically exact forms, e.g., ``g' = g/(1+g)``) |
| `Optics.increment_2stream` | `mo_optical_props_kernels`: `increment_2stream_by_2stream` | — |
| `Fluxes.compute_net_flux!` | `mo_fluxes_broadband_kernels`: `net_broadband` | ``F_{\mathrm{net}} = F^↑ - F^↓`` |

## Gas optics

The Fortran gas-optics kernels live in
`rrtmgp-kernels/mo_gas_optics_rrtmgp_kernels.F90`; the RRTMGP.jl equivalents in
`src/optics/gas_optics.jl` and `src/optics/compute_optical_props.jl`.

| RRTMGP.jl | Fortran rte-rrtmgp | Papers / symbols |
|:----------|:-------------------|:-----------------|
| `Optics.compute_optical_props!` | `ty_gas_optics_rrtmgp`: `gas_optics()` | [pincus2019](@cite) |
| `compute_gas_optics` (per g-point and column) | `compute_tau_absorption` + `compute_Planck_source` | ``τ``, Planck fraction |
| `compute_interp_frac_temp` / `_press` / `_η` | `interpolation` | temperature/pressure/mixing-fraction interpolation weights (`jtemp`, `jpress`, `jeta`, `fmajor`, `fminor`); ``η`` is the binary-species mixing fraction |
| `compute_τ_minor` | `gas_optical_depths_minor` | minor-gas absorption |
| `compute_τ_rayleigh` | `compute_tau_rayleigh` | Rayleigh scattering ``τ`` |
| `interp1d_equispaced`, `interp2d`, `interp3d` | `interpolate1D`, `interpolate2D`, `interpolate3D` | — |
| `Optics.compute_col_gas!` | `ty_gas_optics_rrtmgp`: `get_col_dry` | dry-air column amount [molecules cm⁻²] |
| `Optics.compute_relative_humidity!` | (host responsibility in Fortran) | feeds the MERRA aerosol lookup |

## Cloud and aerosol optics

| RRTMGP.jl | Fortran rte-rrtmgp | Papers / symbols |
|:----------|:-------------------|:-----------------|
| `add_cloud_optics_1scalar!` / `add_cloud_optics_2stream!` with `compute_lookup_cld_liq_props` / `compute_lookup_cld_ice_props` | `ty_cloud_optics_rrtmgp`: `cloud_optics()` (lookup-table path) | liquid/ice ``τ``, ``ω₀``, ``g`` from water path and effective radius |
| `build_cloud_mask!` with `MaxRandomOverlap` | `extensions/mo_cloud_sampling`: `sampled_mask` (maximum-random overlap) | cloud-overlap sampling |
| `add_aerosol_optics_1scalar!` / `add_aerosol_optics_2stream!` with `compute_lookup_aerosol` | `ty_aerosol_optics_rrtmgp_merra`: `aerosol_optics()` | MERRA-2 aerosol optics by type, radius bin, and relative humidity |

## RRTMGP.jl-only layers

These represent Julia-specific additions; the closest analog is the example driver
programs (`examples/rfmip-clear-sky`, `examples/all-sky`), which each user of
the Fortran code re-writes.

| RRTMGP.jl | Role |
|:----------|:-----|
| `RRTMGPSolver`, `update_fluxes!`, `prepare_atmosphere!`, the named getters | Layer-2 host convenience: construct once, drive the functional core each radiation step (see [The getter contract](@ref)) |
| `standard_atmosphere`, `solve(profile)`, `RadiationOutput`, `solve_gray` | Layer-3 standalone front door (teaching, single-column experiments) |
| `GrayAtmosphere` module and the gray optics kernels | analytic gray-atmosphere radiation [schneider2004](@cite), [frierson2006](@cite), [ogorman2008](@cite) |
| `Numerics` module (`k_min`, `τ_thresh`, `resonance_window`, `μ₀_min`) | documented numerical guard constants (hard-coded literals in the Fortran kernels) |
| `LookupBundle` caching via `save_lookup_tables` / `load_lookup_tables` | NetCDF-free lookup reuse |
