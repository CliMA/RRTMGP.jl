#####
##### Deprecated names
#####
#
# Renames made while settling the public API for v1. Each old name forwards to
# its replacement with a deprecation warning; remove them in the next breaking
# release. `Base.@deprecate` cannot carry a docstring, so the replacements are
# named here:
#
#   clear_lw_flux                        -> clear_lw_flux_net
#   clear_sw_flux                        -> clear_sw_flux_net
#       (matching lw_flux_net / sw_flux_net)
#   toa_flux                             -> toa_sw_flux_dn
#   top_of_atmosphere_lw_flux_dn         -> toa_lw_flux_dn
#   top_of_atmosphere_diffuse_sw_flux_dn -> toa_diffuse_sw_flux_dn
#       (one `toa_` prefix for the three incident-flux getters)

Base.@deprecate clear_lw_flux(s::RRTMGPSolver) clear_lw_flux_net(s) false
Base.@deprecate clear_sw_flux(s::RRTMGPSolver) clear_sw_flux_net(s) false
Base.@deprecate toa_flux(s::RRTMGPSolver) toa_sw_flux_dn(s) false
Base.@deprecate top_of_atmosphere_lw_flux_dn(s::RRTMGPSolver) toa_lw_flux_dn(s) false
Base.@deprecate top_of_atmosphere_diffuse_sw_flux_dn(s::RRTMGPSolver) toa_diffuse_sw_flux_dn(
    s,
) false
