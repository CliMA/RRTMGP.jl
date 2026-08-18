#####
##### Layer-2 getters: the host data-exchange contract
#####
#
# Named getters over an `RRTMGPSolver`. Every getter carrying a vertical dimension returns a
# domain-masked `view` into a solver-owned buffer (physical dimension first, `ncol` last), so
# a host reads and writes RRTMGP's buffers without touching the packed layout -- which is how
# RRTMGP stays free of ClimaCore. The full contract (layout, boundary-layer masking,
# writability, and the deliberate exceptions: the well-mixed `volume_mixing_ratio` scalar
# and the allocating `heating_rate`) is documented in `docs/src/getters.md`.

import ..AtmosphericStates: getview_p_lay, getview_t_lay, getview_rel_hum
import ..VolumeMixingRatios

#! format: off

_maybe_transpose(x) = transpose(x)
_maybe_transpose(x::Nothing) = x

# internal API methods
_lookup_tables(s::RRTMGPSolver) = s.lookups
_radiation_method(s::RRTMGPSolver) = s.radiation_method
_shortwave_solver(s::RRTMGPSolver) = s.sws
_longwave_solver(s::RRTMGPSolver) = s.lws
_atmospheric_state(s::RRTMGPSolver) = s.as

# Return the domain-only view of a vertical (dim-1) array, excluding the extra
# isothermal boundary layer/level when one is present. Every getter that carries a
# vertical dimension is wrapped in this, so callers see arrays sized to the physical
# domain rather than RRTMGP's internal boundary-extended arrays. `nothing` passes
# through (for optional fields such as `deep_atmosphere_inverse_scaling` and
# `center_z`/`face_z`). Internal code that needs the full (boundary-extended) array
# uses the raw fields, or `parent(getter(s))`.
#
# This ALWAYS returns a `view` (spanning the whole array when there is no boundary
# layer) rather than branching to the bare array, so both cases share one concrete
# return type. That type stability is what keeps the getters -- and the Layer-2
# `update_fluxes!`/`update_net_fluxes!` broadcasts built on them -- allocation-free and
# JET-clean (see the @allocated/@test_opt tests in test/standalone.jl). `Int(bl)` ∈ {0,1}.
_domain_view(s::RRTMGPSolver, x) = _domain_view(s.grid_params.isothermal_boundary_layer, x)
_domain_view(::Bool, ::Nothing) = nothing
_domain_view(isothermal_boundary_layer::Bool, x::AbstractArray) =
    view(x, 1:(size(x, 1) - Int(isothermal_boundary_layer)), :)
# 3D variant for per-band `(nlev, ncol, n_bnd)` fluxes (mask the vertical dim only).
_domain_view(isothermal_boundary_layer::Bool, x::AbstractArray{<:Any, 3}) =
    view(x, 1:(size(x, 1) - Int(isothermal_boundary_layer)), :, :)

# The clear-sky flux workspaces exist only on `AllSkyRadiationWithClearSkyDiagnostics`
# solvers; give the `clear_*` getters an informative error elsewhere (instead of a
# cryptic `getproperty(::Nothing, ...)`).
_require_clear_sky(::Nothing) = error(
    "clear-sky diagnostic fluxes were not retained; construct the `RRTMGPSolver` \
     with `AllSkyRadiationWithClearSkyDiagnostics`.",
)
_require_clear_sky(x) = x

# Potentially views. Each getter carries a one-line docstring for REPL help; the shared
# layout/masking/writability rules live in the getter contract (docs/src/getters.md).

"""
    top_of_atmosphere_lw_flux_dn(s::RRTMGPSolver)

Return the prescribed incident longwave flux [W/m²], `(ngpt, ncol)`, or `nothing`.
"""
top_of_atmosphere_lw_flux_dn(s::RRTMGPSolver)         = _maybe_transpose(s.lws.bcs.inc_flux)

"""
    top_of_atmosphere_diffuse_sw_flux_dn(s::RRTMGPSolver)

Return the prescribed incident diffuse shortwave flux [W/m²], `(ngpt, ncol)`, or `nothing`.
"""
top_of_atmosphere_diffuse_sw_flux_dn(s::RRTMGPSolver) = _maybe_transpose(s.sws.bcs.inc_flux_diffuse)

"""
    lw_flux_up(s::RRTMGPSolver)

Return the upward longwave flux [W/m²]: a domain-masked `(nlev, ncol)` view.
"""
lw_flux_up(s::RRTMGPSolver)                           = _domain_view(s, s.presented_flux_lw.flux_up)

"""
    lw_flux_dn(s::RRTMGPSolver)

Return the downward longwave flux [W/m²]: a domain-masked `(nlev, ncol)` view.
"""
lw_flux_dn(s::RRTMGPSolver)                           = _domain_view(s, s.presented_flux_lw.flux_dn)

"""
    lw_flux_net(s::RRTMGPSolver)

Return the net (up - down) longwave flux [W/m²]: a domain-masked `(nlev, ncol)` view.
"""
lw_flux_net(s::RRTMGPSolver)                          = _domain_view(s, s.presented_flux_lw.flux_net)

"""
    clear_lw_flux_up(s::RRTMGPSolver)

Return the clear-sky upward longwave flux [W/m²], `(nlev, ncol)`
(`AllSkyRadiationWithClearSkyDiagnostics` only).
"""
clear_lw_flux_up(s::RRTMGPSolver)                     = _domain_view(s, _require_clear_sky(s.clear_flux_lw).flux_up)

"""
    clear_lw_flux_dn(s::RRTMGPSolver)

Return the clear-sky downward longwave flux [W/m²], `(nlev, ncol)`
(`AllSkyRadiationWithClearSkyDiagnostics` only).
"""
clear_lw_flux_dn(s::RRTMGPSolver)                     = _domain_view(s, _require_clear_sky(s.clear_flux_lw).flux_dn)

"""
    clear_lw_flux(s::RRTMGPSolver)

Return the clear-sky net longwave flux [W/m²], `(nlev, ncol)`
(`AllSkyRadiationWithClearSkyDiagnostics` only).
"""
clear_lw_flux(s::RRTMGPSolver)                        = _domain_view(s, _require_clear_sky(s.clear_flux_lw).flux_net)

"""
    surface_emissivity(s::RRTMGPSolver)

Return the longwave surface emissivity [-]: a writable `(nbnd, ncol)` device array.
"""
surface_emissivity(s::RRTMGPSolver)                   = s.lws.bcs.sfc_emis

"""
    sw_flux_up(s::RRTMGPSolver)

Return the upward shortwave flux [W/m²]: a domain-masked `(nlev, ncol)` view.
"""
sw_flux_up(s::RRTMGPSolver)                           = _domain_view(s, s.presented_flux_sw.flux_up)

"""
    sw_flux_dn(s::RRTMGPSolver)

Return the downward shortwave flux [W/m²]: a domain-masked `(nlev, ncol)` view.
"""
sw_flux_dn(s::RRTMGPSolver)                           = _domain_view(s, s.presented_flux_sw.flux_dn)

"""
    sw_flux_net(s::RRTMGPSolver)

Return the net (up - down) shortwave flux [W/m²]: a domain-masked `(nlev, ncol)` view.
"""
sw_flux_net(s::RRTMGPSolver)                          = _domain_view(s, s.presented_flux_sw.flux_net)

"""
    sw_direct_flux_dn(s::RRTMGPSolver)

Return the direct-beam downward shortwave flux [W/m²]: a domain-masked `(nlev, ncol)` view.
"""
sw_direct_flux_dn(s::RRTMGPSolver)                    = _domain_view(s, s.presented_flux_sw.flux_dn_dir)

"""
    clear_sw_flux_up(s::RRTMGPSolver)

Return the clear-sky upward shortwave flux [W/m²], `(nlev, ncol)`
(`AllSkyRadiationWithClearSkyDiagnostics` only).
"""
clear_sw_flux_up(s::RRTMGPSolver)                     = _domain_view(s, _require_clear_sky(s.clear_flux_sw).flux_up)

"""
    clear_sw_flux_dn(s::RRTMGPSolver)

Return the clear-sky downward shortwave flux [W/m²], `(nlev, ncol)`
(`AllSkyRadiationWithClearSkyDiagnostics` only).
"""
clear_sw_flux_dn(s::RRTMGPSolver)                     = _domain_view(s, _require_clear_sky(s.clear_flux_sw).flux_dn)

"""
    clear_sw_direct_flux_dn(s::RRTMGPSolver)

Return the clear-sky direct-beam downward shortwave flux [W/m²], `(nlev, ncol)`
(`AllSkyRadiationWithClearSkyDiagnostics` only).
"""
clear_sw_direct_flux_dn(s::RRTMGPSolver)              = _domain_view(s, _require_clear_sky(s.clear_flux_sw).flux_dn_dir)

"""
    clear_sw_flux(s::RRTMGPSolver)

Return the clear-sky net shortwave flux [W/m²], `(nlev, ncol)`
(`AllSkyRadiationWithClearSkyDiagnostics` only).
"""
clear_sw_flux(s::RRTMGPSolver)                        = _domain_view(s, _require_clear_sky(s.clear_flux_sw).flux_net)
"""
    cloud_liquid_effective_radius(s::RRTMGPSolver)

Return the cloud liquid effective radius [µm]: a writable, domain-masked `(nlay, ncol)` view.
"""
cloud_liquid_effective_radius(s::RRTMGPSolver)        = _domain_view(s, s.as.cloud_state.cld_r_eff_liq)

"""
    cloud_ice_effective_radius(s::RRTMGPSolver)

Return the cloud ice effective radius [µm]: a writable, domain-masked `(nlay, ncol)` view.
"""
cloud_ice_effective_radius(s::RRTMGPSolver)           = _domain_view(s, s.as.cloud_state.cld_r_eff_ice)

"""
    cloud_liquid_water_path(s::RRTMGPSolver)

Return the cloud liquid water path [g/m²]: a writable, domain-masked `(nlay, ncol)` view.
"""
cloud_liquid_water_path(s::RRTMGPSolver)              = _domain_view(s, s.as.cloud_state.cld_path_liq)

"""
    cloud_ice_water_path(s::RRTMGPSolver)

Return the cloud ice water path [g/m²]: a writable, domain-masked `(nlay, ncol)` view.
"""
cloud_ice_water_path(s::RRTMGPSolver)                 = _domain_view(s, s.as.cloud_state.cld_path_ice)

"""
    cloud_fraction(s::RRTMGPSolver)

Return the cloud fraction [-]: a writable, domain-masked `(nlay, ncol)` view.
"""
cloud_fraction(s::RRTMGPSolver)                       = _domain_view(s, s.as.cloud_state.cld_frac)

"""
    sw_cloud_cover(s::RRTMGPSolver)

Return the column cloud cover [-] diagnosed from the shortwave McICA sampling, `(ncol,)`.
"""
sw_cloud_cover(s::RRTMGPSolver)                       = s.as.cloud_state.cld_cover_sw

"""
    lw_cloud_cover(s::RRTMGPSolver)

Return the column cloud cover [-] diagnosed from the longwave McICA sampling, `(ncol,)`.
"""
lw_cloud_cover(s::RRTMGPSolver)                       = s.as.cloud_state.cld_cover_lw

"""
    aod_sw_extinction(s::RRTMGPSolver)

Return the aerosol extinction optical depth [-] in the shortwave band containing 550 nm,
`(ncol,)`.
"""
aod_sw_extinction(s::RRTMGPSolver)                    = s.as.aerosol_state.aod_sw_ext

"""
    aod_sw_scattering(s::RRTMGPSolver)

Return the aerosol scattering optical depth [-] in the shortwave band containing 550 nm,
`(ncol,)`.
"""
aod_sw_scattering(s::RRTMGPSolver)                    = s.as.aerosol_state.aod_sw_sca

"""
    center_z(s::RRTMGPSolver)

Return the layer-center altitudes [m], `(nlay, ncol)`, or `nothing` if not provided.
"""
center_z(s::RRTMGPSolver)                             = _domain_view(s, s.center_z)

"""
    face_z(s::RRTMGPSolver)

Return the level (face) altitudes [m], `(nlev, ncol)`, or `nothing` if not provided.
"""
face_z(s::RRTMGPSolver)                               = _domain_view(s, s.face_z)
"""
    cos_zenith(s::RRTMGPSolver)

Return the cosine of the solar zenith angle [-]: a writable `(ncol,)` device array.
"""
cos_zenith(s::RRTMGPSolver)                           = s.sws.bcs.cos_zenith

"""
    toa_flux(s::RRTMGPSolver)

Return the incident top-of-atmosphere solar flux [W/m²]: a writable `(ncol,)` device array.
"""
toa_flux(s::RRTMGPSolver)                             = s.sws.bcs.toa_flux

"""
    direct_sw_surface_albedo(s::RRTMGPSolver)

Return the direct-beam shortwave surface albedo [-]: a writable `(nbnd, ncol)` device array.
"""
direct_sw_surface_albedo(s::RRTMGPSolver)             = s.sws.bcs.sfc_alb_direct

"""
    diffuse_sw_surface_albedo(s::RRTMGPSolver)

Return the diffuse shortwave surface albedo [-]: a writable `(nbnd, ncol)` device array.
"""
diffuse_sw_surface_albedo(s::RRTMGPSolver)            = s.sws.bcs.sfc_alb_diffuse

"""
    latitude(s::RRTMGPSolver)

Return the column latitudes [degrees], `(ncol,)`, or `nothing` if not provided.
"""
latitude(s::RRTMGPSolver)                             = s.as.lat

"""
    surface_temperature(s::RRTMGPSolver)

Return the surface temperature [K]: a writable `(ncol,)` device array.
"""
surface_temperature(s::RRTMGPSolver)                  = s.as.t_sfc

"""
    layer_pressure(s::RRTMGPSolver)

Return the layer-center pressure [Pa]: a writable, domain-masked `(nlay, ncol)` view.
"""
layer_pressure(s::RRTMGPSolver)                       = _domain_view(s, getview_p_lay(s.as))

"""
    layer_temperature(s::RRTMGPSolver)

Return the layer-center temperature [K]: a writable, domain-masked `(nlay, ncol)` view.
"""
layer_temperature(s::RRTMGPSolver)                    = _domain_view(s, getview_t_lay(s.as))

"""
    layer_relative_humidity(s::RRTMGPSolver)

Return the layer-center relative humidity [-]: a writable, domain-masked `(nlay, ncol)`
view. Feeds the humidity-dependent aerosol optics; keeping it current is the host's job.
"""
layer_relative_humidity(s::RRTMGPSolver)              = _domain_view(s, getview_rel_hum(s.as))

"""
    level_pressure(s::RRTMGPSolver)

Return the level (face) pressure [Pa]: a writable, domain-masked `(nlev, ncol)` view.
"""
level_pressure(s::RRTMGPSolver)                       = _domain_view(s, s.as.p_lev)

"""
    level_temperature(s::RRTMGPSolver)

Return the level (face) temperature [K]: a writable, domain-masked `(nlev, ncol)` view.
"""
level_temperature(s::RRTMGPSolver)                    = _domain_view(s, s.as.t_lev)

"""
    net_flux(s::RRTMGPSolver)

Return the combined longwave + shortwave net flux [W/m²]: a domain-masked `(nlev, ncol)`
view.
"""
net_flux(s::RRTMGPSolver)                             = _domain_view(s, s.net_flux_buffer)

"""
    clear_net_flux(s::RRTMGPSolver)

Return the clear-sky combined net flux [W/m²], `(nlev, ncol)`
(`AllSkyRadiationWithClearSkyDiagnostics` only).
"""
clear_net_flux(s::RRTMGPSolver)                       = _domain_view(s, _require_clear_sky(s.clear_net_flux_buffer))

"""
    deep_atmosphere_inverse_scaling(s::RRTMGPSolver)

Return the deep-atmosphere flux scaling [-], `(nlev, ncol)`, or `nothing` if not provided.
"""
deep_atmosphere_inverse_scaling(s::RRTMGPSolver)      = _domain_view(s, s.deep_atmosphere_inverse_scaling)
#! format: on

"""
    radiation_method(s::RRTMGPSolver)

Return the radiation method the solver was constructed with: `GrayRadiation`,
`ClearSkyRadiation`, `AllSkyRadiation`, or `AllSkyRadiationWithClearSkyDiagnostics`.
"""
radiation_method(s::RRTMGPSolver) = _radiation_method(s)

# Gray radiation carries no spectral lookup tables; raise an informative error for
# getters that need them (band bounds, gas indices).
function _require_spectral_lookups(s::RRTMGPSolver)
    lookups = _lookup_tables(s)
    isnothing(lookups.lookup_lw) && error(
        "spectral lookup tables are not present on this solver (gray radiation); \
         construct it with a spectral radiation method (e.g., `ClearSkyRadiation`).",
    )
    return lookups
end

# The aerosol index map exists only when the solver was built with aerosol optics.
function _require_aerosol_lookups(s::RRTMGPSolver)
    lookups = _lookup_tables(s)
    isnothing(lookups.idx_aerosol_sw) && error(
        "aerosol lookup tables are not present on this solver; construct it with a \
         radiation method with `aerosol_radiation = true`.",
    )
    return lookups
end

# Internal helpers for the optional per-band flux buffers.
_require_band_flux(::Nothing) = error(
    "spectral fluxes were not retained; construct the `RRTMGPSolver` with `spectral_fluxes = true`.",
)
_require_band_flux(band::Fluxes.FluxBand) = band
_solver_band_flux(ws) =
    hasproperty(ws, :band_flux) ? _require_band_flux(ws.band_flux) :
    error("spectral fluxes require a two-stream, non-gray solver.")

# The six per-band getters share one contract: a domain-masked
# `(nlev, ncol, n_bnd)` view of a retained band buffer, available only when the
# solver was built with `spectral_fluxes = true`. Each carries its own
# docstring so REPL help works for all of them; the shared detail is on the
# "Get per-band (spectral) fluxes" how-to page.

"""
    spectral_lw_flux_up(s::RRTMGPSolver)

Return the per-band upward longwave flux [W/m²]: a domain-masked
`(nlev, ncol, n_bnd)` view. Requires `spectral_fluxes = true`; see
[`lw_band_bounds`](@ref) for each band's wavenumber range.
"""
spectral_lw_flux_up(s::RRTMGPSolver) =
    _domain_view(s, _solver_band_flux(s.lws).flux_up)

"""
    spectral_lw_flux_dn(s::RRTMGPSolver)

Return the per-band downward longwave flux [W/m²]: a domain-masked
`(nlev, ncol, n_bnd)` view. Requires `spectral_fluxes = true`; see
[`lw_band_bounds`](@ref) for each band's wavenumber range.
"""
spectral_lw_flux_dn(s::RRTMGPSolver) =
    _domain_view(s, _solver_band_flux(s.lws).flux_dn)

"""
    spectral_lw_flux_net(s::RRTMGPSolver)

Return the per-band net (up - down) longwave flux [W/m²]: a domain-masked
`(nlev, ncol, n_bnd)` view. Requires `spectral_fluxes = true`; see
[`lw_band_bounds`](@ref) for each band's wavenumber range.
"""
spectral_lw_flux_net(s::RRTMGPSolver) =
    _domain_view(s, _solver_band_flux(s.lws).flux_net)

"""
    spectral_sw_flux_up(s::RRTMGPSolver)

Return the per-band upward shortwave flux [W/m²]: a domain-masked
`(nlev, ncol, n_bnd)` view. Requires `spectral_fluxes = true`; see
[`sw_band_bounds`](@ref) for each band's wavenumber range.
"""
spectral_sw_flux_up(s::RRTMGPSolver) =
    _domain_view(s, _solver_band_flux(s.sws).flux_up)

"""
    spectral_sw_flux_dn(s::RRTMGPSolver)

Return the per-band downward shortwave flux [W/m²]: a domain-masked
`(nlev, ncol, n_bnd)` view. Requires `spectral_fluxes = true`; see
[`sw_band_bounds`](@ref) for each band's wavenumber range.
"""
spectral_sw_flux_dn(s::RRTMGPSolver) =
    _domain_view(s, _solver_band_flux(s.sws).flux_dn)

"""
    spectral_sw_flux_net(s::RRTMGPSolver)

Return the per-band net (up - down) shortwave flux [W/m²]: a domain-masked
`(nlev, ncol, n_bnd)` view. Requires `spectral_fluxes = true`; see
[`sw_band_bounds`](@ref) for each band's wavenumber range.
"""
spectral_sw_flux_net(s::RRTMGPSolver) =
    _domain_view(s, _solver_band_flux(s.sws).flux_net)

"""
    lw_band_bounds(s::RRTMGPSolver)

Return the `(2, n_bnd)` lower/upper wavenumber edges [cm⁻¹] of the longwave bands,
identifying the spectral range of each band in the per-band fluxes (e.g.,
[`spectral_lw_flux_up`](@ref)).
"""
lw_band_bounds(s::RRTMGPSolver) =
    _require_spectral_lookups(s).lookup_lw.band_data.bnd_lims_wn

"""
    sw_band_bounds(s::RRTMGPSolver)

Return the `(2, n_bnd)` lower/upper wavenumber edges [cm⁻¹] of the shortwave bands,
identifying the spectral range of each band in
[the per-band fluxes](@ref RRTMGP.spectral_lw_flux_up) (`spectral_sw_flux_up`
and companions).
"""
sw_band_bounds(s::RRTMGPSolver) =
    _require_spectral_lookups(s).lookup_sw.band_data.bnd_lims_wn

"""
    optical_thickness_parameter(s::RRTMGPSolver)

For a gray-radiation solver, return the gray optical-thickness parameters (an
`AbstractGrayOpticalThickness`, e.g., `GrayOpticalThicknessOGorman2008`) the atmospheric
state was built with; `nothing` for lookup-table (non-gray) radiation, which has no such
parameter.
"""
optical_thickness_parameter(s::RRTMGPSolver) =
    hasproperty(_atmospheric_state(s), :otp) ? _atmospheric_state(s).otp : nothing

"""
    isothermal_boundary_layer(s::RRTMGPSolver)

Return whether the solver carries an internal isothermal boundary layer above the
physical domain (which every getter masks off).
"""
isothermal_boundary_layer(s::RRTMGPSolver) = s.grid_params.isothermal_boundary_layer

"""
    aerosol_radius(s::RRTMGPSolver, name::AbstractString)

Return the aerosol radius for the given aerosol name.

Available names are:
$(aerosol_names_docs())
"""
function aerosol_radius(s::RRTMGPSolver, name::AbstractString)
    lookups = _require_aerosol_lookups(s) # informative error before field access
    return _domain_view(s, view(
        s.as.aerosol_state.aero_size,
        lookups.idx_aerosol_sw[canonical_aerosol_name(name)],
        :,
        :
    ))
end

"""
    aerosol_column_mass_density(s::RRTMGPSolver, name::AbstractString)

Return the aerosol column mass density [kg/m²] for the given aerosol name.

Available names are:
$(aerosol_names_docs())
"""
function aerosol_column_mass_density(s::RRTMGPSolver, name::AbstractString)
    lookups = _require_aerosol_lookups(s) # informative error before field access
    return _domain_view(s, view(
        s.as.aerosol_state.aero_mass,
        lookups.idx_aerosol_sw[canonical_aerosol_name(name)],
        :,
        :
    ))
end

"""
    gas_names_sw()

Return a vector containing the gas names in the shortwave lookup tables.

This should be the same list of gases returned from the following code,
and is tested in `lookup_tables`
```julia
function gas_names_sw_from_artifacts()
   artifact(t, b, n) =
       NC.Dataset(RRTMGP.ArtifactPaths.get_lookup_filename(t, b)) do ds
           getproperty(RRTMGP.LookUpTables, n)(ds, Float64, Array)
       end
   _, idx_gases_sw = artifact(:gas, :sw, :LookUpSW)
   return keys(idx_gases_sw)
end
```
"""
gas_names_sw() = [
    "h2o",
    "cfc11",
    "h2o_self",
    "co2",
    "cfc12",
    "hfc134a",
    "cfc22",
    "ch4",
    "hfc23",
    "ccl4",
    "hfc143a",
    "co",
    "no2",
    "n2",
    "o2",
    "o3",
    "h2o_frgn",
    "hfc32",
    "n2o",
    "cf4",
    "hfc125"
]
gas_names_sw_docs() = map(x->"`$x`", gas_names_sw())

"""
    volume_mixing_ratio(s::RRTMGPSolver, name::AbstractString)

Return the volume mixing ratio for gas `name`. `"h2o"` and `"o3"` vary by layer
and column (a domain-masked `(nlay, ncol)` view); with the global-mean `VmrGM`
storage every other (well-mixed) gas is a single scalar.

A well-mixed scalar is returned as a host `Number` copied off the device, so each call triggers
a device→host synchronization on the GPU. Read these during setup to maintain performance.

Available names are:
$(gas_names_sw_docs())
"""
volume_mixing_ratio(s::RRTMGPSolver, name::AbstractString) =
    _volume_mixing_ratio(s, _atmospheric_state(s).vmr, name)

function _volume_mixing_ratio(
    s::RRTMGPSolver,
    vmr::VolumeMixingRatios.VmrGM,
    name::AbstractString,
)
    name == "h2o" && return _domain_view(s, vmr.vmr_h2o)
    name == "o3" && return _domain_view(s, vmr.vmr_o3)
    idx = _lookup_tables(s).idx_gases_sw[name]
    # The water-vapor continuum pseudo-gases ("h2o_self"/"h2o_frgn") share the h2o
    # slot: the loader maps them to `idx_h2o` and `get_vmr` maps index 1 to `vmr_h2o`,
    # so return the h2o field rather than the (unused) well-mixed slot 1.
    idx == 1 && return _domain_view(s, vmr.vmr_h2o)
    # A well-mixed gas is a single value; return it as a host scalar. (A 0-d device
    # view would scalar-index the GPU when read.)
    return only(Array(view(vmr.vmr, idx:idx)))
end
_volume_mixing_ratio(s::RRTMGPSolver, vmr::VolumeMixingRatios.Vmr, name::AbstractString) =
    _domain_view(s, view(vmr.vmr, _lookup_tables(s).idx_gases_sw[name], :, :))

"""
    set_volume_mixing_ratio!(s::RRTMGPSolver, name::AbstractString, value)

Set the volume mixing ratio for gas `name` to `value`, returning `value`.

This is the write counterpart to [`volume_mixing_ratio`](@ref). For `"h2o"`/`"o3"`
(and, with per-layer `Vmr` storage, every gas), `value` broadcasts over the
`(nlay, ncol)` field, so it may be a scalar or a domain-sized array. For a
well-mixed gas with the default global-mean (`VmrGM`) storage, the mixing ratio
is a single scalar; pass a scalar `value`.

Unlike `volume_mixing_ratio`, this is the supported way to update a well-mixed
gas: that getter returns a read-only host copy, so `volume_mixing_ratio(s, name) .= value`
does not write back. Use this to update time-varying trace gases (e.g. prescribed CO₂).

Available names are:
$(gas_names_sw_docs())
"""
set_volume_mixing_ratio!(s::RRTMGPSolver, name::AbstractString, value) =
    _set_volume_mixing_ratio!(s, _atmospheric_state(s).vmr, name, value)

function _set_volume_mixing_ratio!(
    s::RRTMGPSolver,
    vmr::VolumeMixingRatios.VmrGM,
    name::AbstractString,
    value,
)
    name == "h2o" && return (_domain_view(s, vmr.vmr_h2o) .= value; value)
    name == "o3" && return (_domain_view(s, vmr.vmr_o3) .= value; value)
    idx = _lookup_tables(s).idx_gases_sw[name]
    # The water-vapor continuum pseudo-gases share the h2o slot (see the getter);
    # route them to `vmr_h2o` rather than the unused well-mixed slot 1.
    idx == 1 && return (_domain_view(s, vmr.vmr_h2o) .= value; value)
    # A well-mixed gas is a single value. `fill!` over a 1-element view writes the
    # scalar without scalar-indexing the GPU (a bare `vmr.vmr[idx] = value` would).
    fill!(view(vmr.vmr, idx:idx), value)
    return value
end
_set_volume_mixing_ratio!(
    s::RRTMGPSolver,
    vmr::VolumeMixingRatios.Vmr,
    name::AbstractString,
    value,
) = (_domain_view(s, view(vmr.vmr, _lookup_tables(s).idx_gases_sw[name], :, :)) .= value; value)
