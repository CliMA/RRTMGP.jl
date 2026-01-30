#! format: off

maybe_transpose(x) = transpose(x)
maybe_transpose(x::Nothing) = x

# TODO: this is not user-facing yet in that `lookup_tables` still returns a
# NamedTuple, and we should probably make a dedicated struct for this.

# internal API methods
_lookup_tables(s::RRTMGPSolver) = s.lookups
_radiation_method(s::RRTMGPSolver) = s.radiation_method
_shortwave_solver(s::RRTMGPSolver) = s.sws
_longwave_solver(s::RRTMGPSolver) = s.lws
_atmospheric_state(s::RRTMGPSolver) = s.as

# Potentially views
top_of_atmosphere_lw_flux_dn(s::RRTMGPSolver)         = maybe_transpose(s.lws.bcs.inc_flux)
top_of_atmosphere_diffuse_sw_flux_dn(s::RRTMGPSolver) = maybe_transpose(s.sws.bcs.inc_flux_diffuse)
lw_flux_up(s::RRTMGPSolver)                           = s.lws.flux.flux_up
lw_flux_dn(s::RRTMGPSolver)                           = s.lws.flux.flux_dn
lw_flux_net(s::RRTMGPSolver)                          = s.lws.flux.flux_net
clear_lw_flux_up(s::RRTMGPSolver)                     = s.clear_flux_lw.flux_up
clear_lw_flux_dn(s::RRTMGPSolver)                     = s.clear_flux_lw.flux_dn
clear_lw_flux(s::RRTMGPSolver)                        = s.clear_flux_lw.flux_net
surface_emissivity(s::RRTMGPSolver)                   = s.lws.bcs.sfc_emis
sw_flux_up(s::RRTMGPSolver)                           = s.sws.flux.flux_up
sw_flux_dn(s::RRTMGPSolver)                           = s.sws.flux.flux_dn
sw_flux_net(s::RRTMGPSolver)                          = s.sws.flux.flux_net
sw_direct_flux_dn(s::RRTMGPSolver)                    = s.sws.flux.flux_dn_dir
clear_sw_flux_up(s::RRTMGPSolver)                     = s.clear_flux_sw.flux_up
clear_sw_flux_dn(s::RRTMGPSolver)                     = s.clear_flux_sw.flux_dn
clear_sw_direct_flux_dn(s::RRTMGPSolver)              = s.clear_flux_sw.flux_dn_dir
clear_sw_flux(s::RRTMGPSolver)                        = s.clear_flux_sw.flux_net
cloud_liquid_effective_radius(s::RRTMGPSolver)        = s.as.cloud_state.cld_r_eff_liq
cloud_ice_effective_radius(s::RRTMGPSolver)           = s.as.cloud_state.cld_r_eff_ice
cloud_liquid_water_path(s::RRTMGPSolver)              = s.as.cloud_state.cld_path_liq
cloud_ice_water_path(s::RRTMGPSolver)                 = s.as.cloud_state.cld_path_ice
cloud_fraction(s::RRTMGPSolver)                       = s.as.cloud_state.cld_frac
aod_sw_extinction(s::RRTMGPSolver)                    = s.as.aerosol_state.aod_sw_ext
aod_sw_scattering(s::RRTMGPSolver)                    = s.as.aerosol_state.aod_sw_sca
center_z(s::RRTMGPSolver)                             = s.center_z
face_z(s::RRTMGPSolver)                               = s.face_z
cos_zenith(s::RRTMGPSolver)                           = s.sws.bcs.cos_zenith
toa_flux(s::RRTMGPSolver)                             = s.sws.bcs.toa_flux
direct_sw_surface_albedo(s::RRTMGPSolver)             = s.sws.bcs.sfc_alb_direct
diffuse_sw_surface_albedo(s::RRTMGPSolver)            = s.sws.bcs.sfc_alb_diffuse
latitude(s::RRTMGPSolver)                             = s.as.lat
surface_temperature(s::RRTMGPSolver)                  = s.as.t_sfc
pressure(s::RRTMGPSolver)                             = getview_p_lay(s.as)
temperature(s::RRTMGPSolver)                          = getview_t_lay(s.as)
relative_humidity(s::RRTMGPSolver)                    = getview_rel_hum(s.as)

#! format: on
