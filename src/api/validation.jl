#####
##### Optional input validation (opt-in, in the spirit of rte-rrtmgp's
##### mo_rte_config check_values)
#####

"""
    check_values

Global toggle for input validation: when enabled with
`RRTMGP.check_values[] = true`, [`update_fluxes!`](@ref) calls
[`validate_inputs`](@ref) before each solve. Off by default — the checks
reduce over device arrays, so they are intended for development and debugging,
not for production time-stepping (with the toggle off they cost one branch and
keep `update_fluxes!` allocation-free).
"""
const check_values = Ref(false)

@noinline _validation_error(name::Symbol) = error(
    "RRTMGP input validation failed: `$name` contains values outside its \
     physical range (or non-finite values). Inspect the corresponding getter \
     on the solver.",
)

_all(f::F, x::AbstractArray) where {F} = mapreduce(f, &, x; init = true)
_check(f::F, x::AbstractArray, name::Symbol) where {F} =
    _all(f, x) || _validation_error(name)
_check(f, x::Nothing, name::Symbol) = true

_in_unit_interval(v) = isfinite(v) & (zero(v) ≤ v ≤ one(v))
_in_cos_range(v) = isfinite(v) & (-one(v) ≤ v ≤ one(v))
_positive_finite(v) = isfinite(v) & (v > zero(v))
_nonnegative_finite(v) = isfinite(v) & (v ≥ zero(v))

_check_vmr(::Nothing) = true
_check_vmr(vmr::VolumeMixingRatios.VmrGM) = begin
    _check(_nonnegative_finite, vmr.vmr_h2o, :vmr_h2o)
    _check(_nonnegative_finite, vmr.vmr_o3, :vmr_o3)
    _check(_nonnegative_finite, vmr.vmr, :vmr)
end
_check_vmr(vmr::VolumeMixingRatios.Vmr) = _check(_nonnegative_finite, vmr.vmr, :vmr)

"""
    validate_inputs(s::RRTMGPSolver)

Validate the solver's host-provided inputs against their physical ranges,
raising an informative error on the first violation:

- level/layer pressures and temperatures: positive and finite,
- `cos_zenith ∈ [-1, 1]`, `toa_flux ≥ 0`,
- surface emissivity and the two shortwave surface albedos ∈ [0, 1],
- gas volume mixing ratios ≥ 0 (spectral states).

Runs automatically inside [`update_fluxes!`](@ref) when
[`check_values`](@ref)`[] = true`; can also be called directly at any time.
Reduces over the solver's device arrays (not allocation-free).
"""
function validate_inputs(s::RRTMGPSolver)
    as = _atmospheric_state(s)
    _check(_positive_finite, as.p_lev, :level_pressure)
    _check(_positive_finite, as.t_lev, :level_temperature)
    _check(_positive_finite, getview_p_lay(as), :layer_pressure)
    _check(_positive_finite, getview_t_lay(as), :layer_temperature)
    _check(_positive_finite, as.t_sfc, :surface_temperature)
    _check(_in_cos_range, s.sws.bcs.cos_zenith, :cos_zenith)
    _check(_nonnegative_finite, s.sws.bcs.toa_flux, :toa_flux)
    _check(_in_unit_interval, s.lws.bcs.sfc_emis, :surface_emissivity)
    _check(_in_unit_interval, s.sws.bcs.sfc_alb_direct, :direct_sw_surface_albedo)
    _check(_in_unit_interval, s.sws.bcs.sfc_alb_diffuse, :diffuse_sw_surface_albedo)
    hasproperty(as, :vmr) && _check_vmr(as.vmr)
    return nothing
end
