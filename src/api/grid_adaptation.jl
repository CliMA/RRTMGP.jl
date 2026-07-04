#####
##### Grid adaptation: prepare an AtmosphericState for the RTE solvers
#####
#
# TODO: Remove comment after ClimaAtmos update:
#
# Optional, separable steps lifted from the ClimaAtmos RRTMGP interface so they
# live with the radiation code and can be used standalone (without ClimaCore).
# Each operates in place on an `AtmosphericState` (or `GrayAtmosphericState`) and
# dispatches on the state / sub-state type, so it is a no-op where not
# applicable. The arithmetic is kept identical to the original to preserve
# bit-for-bit reproducibility of downstream fluxes.

import ..AtmosphericStates:
    AbstractAtmosphericState,
    AtmosphericState,
    GrayAtmosphericState,
    CloudState,
    AerosolState,
    getview_p_lay,
    getview_t_lay,
    getview_rel_hum,
    getview_col_dry
import ..Vmrs
import ..Optics

"""
    get_p_min(as, lookup_lw)

Return the minimum pressure supported by the radiation scheme: zero for gray
radiation, and the longwave lookup table's reference minimum pressure otherwise.
"""
get_p_min(as::GrayAtmosphericState, lookup_lw) = zero(eltype(as.p_lay))
get_p_min(as::AtmosphericState, lookup_lw) = lookup_lw.p_ref_min

#####
##### Level <-> layer interpolation
#####

"""
    interpolate_levels!(as, interpolation, bottom_extrapolation, params;
                        center_z = nothing, face_z = nothing,
                        isothermal_boundary_layer = false)

Fill the level (cell-face) pressures and temperatures of `as` by interpolating
and extrapolating its layer (cell-center) values: `interpolation` is used on the
interior faces and the top face, and `bottom_extrapolation` on the bottom face.
A no-op for `NoInterpolation` (the caller is assumed to have provided level
values directly). `center_z` / `face_z` are required only for z-based methods
(`BestFit`, `HydrostaticBottom`). Mutates `as` in place and returns it.
"""
interpolate_levels!(
    as::AbstractAtmosphericState,
    ::NoInterpolation,
    ::AbstractBottomExtrapolation,
    params;
    kwargs...,
) = as

function interpolate_levels!(
    as::AbstractAtmosphericState,
    interpolation::AbstractInterpolation,
    bottom_extrapolation::AbstractBottomExtrapolation,
    params;
    center_z = nothing,
    face_z = nothing,
    isothermal_boundary_layer::Bool = false,
)
    p_lev = as.p_lev
    t_lev = as.t_lev
    t_sfc = as.t_sfc
    p_lay = getview_p_lay(as)
    t_lay = getview_t_lay(as)
    nlay = size(p_lay, 1) - Int(isothermal_boundary_layer)
    if requires_z(interpolation) || requires_z(bottom_extrapolation)
        z_lay = parent(center_z)
        z_lev = parent(face_z)
    end
    mode = interpolation
    outs = requires_z(mode) ? (p_lev, t_lev, z_lev) : (p_lev, t_lev)
    ins = requires_z(mode) ? (p_lay, t_lay, z_lay) : (p_lay, t_lay)
    update_views(interp!, mode, outs, ins, (), 2:nlay, 1:(nlay - 1), 2:nlay)
    others = (t_sfc, params)
    update_views(extrap!, mode, outs, ins, others, nlay + 1, nlay, nlay - 1)
    mode =
        bottom_extrapolation isa SameAsInterpolation ? interpolation :
        bottom_extrapolation
    outs = requires_z(mode) ? (p_lev, t_lev, z_lev) : (p_lev, t_lev)
    ins = requires_z(mode) ? (p_lay, t_lay, z_lay) : (p_lay, t_lay)
    update_views(extrap!, mode, outs, ins, others, 1, 1, 2)
    return as
end

update_views(f, mode, outs, ins, others, out_range, in_range1, in_range2) = f(
    mode,
    map(out -> view(out, out_range, :), outs)...,
    map(in -> view(in, in_range1, :), ins)...,
    map(in -> view(in, in_range2, :), ins)...,
    others...,
)

#####
##### Isothermal boundary layer
#####

"""
    add_isothermal_boundary_layer!(as, p_min)

Fill the extra (top) isothermal layer of `as`: its top level pressure is set to
`p_min`, while its layer/level temperatures, relative humidity, volume mixing
ratios, cloud properties, and aerosol properties are copied from the layer below.
Mutates `as` in place and returns it.
"""
function add_isothermal_boundary_layer!(as::AtmosphericState, p_min)
    p_lev = as.p_lev
    t_lev = as.t_lev
    p_lay = getview_p_lay(as)
    t_lay = getview_t_lay(as)
    rel_hum = getview_rel_hum(as)
    @views p_lay[end, :] .= (p_lev[end - 1, :] .+ p_min) ./ 2
    @views p_lev[end, :] .= p_min
    @views t_lay[end, :] .= t_lev[end - 1, :]
    @views t_lev[end, :] .= t_lev[end - 1, :]
    @views rel_hum[end, :] .= rel_hum[end - 1, :]
    _extend_boundary_vmr!(as.vmr)
    _extend_boundary_cloud!(as.cloud_state)
    _extend_boundary_aerosol!(as.aerosol_state)
    return as
end

function add_isothermal_boundary_layer!(as::GrayAtmosphericState, p_min)
    p_lev = as.p_lev
    t_lev = as.t_lev
    p_lay = getview_p_lay(as)
    t_lay = getview_t_lay(as)
    @views p_lay[end, :] .= (p_lev[end - 1, :] .+ p_min) ./ 2
    @views p_lev[end, :] .= p_min
    @views t_lay[end, :] .= t_lev[end - 1, :]
    @views t_lev[end, :] .= t_lev[end - 1, :]
    return as
end

function _extend_boundary_vmr!(vmr::Vmrs.VmrGM)
    @views vmr.vmr_h2o[end, :] .= vmr.vmr_h2o[end - 1, :]
    @views vmr.vmr_o3[end, :] .= vmr.vmr_o3[end - 1, :]
    return nothing
end
function _extend_boundary_vmr!(vmr::Vmrs.Vmr)
    @views vmr.vmr[:, end, :] .= vmr.vmr[:, end - 1, :]
    return nothing
end

_extend_boundary_cloud!(::Nothing) = nothing
function _extend_boundary_cloud!(cloud_state::CloudState)
    @views cloud_state.cld_r_eff_liq[end, :] .=
        cloud_state.cld_r_eff_liq[end - 1, :]
    @views cloud_state.cld_r_eff_ice[end, :] .=
        cloud_state.cld_r_eff_ice[end - 1, :]
    @views cloud_state.cld_path_liq[end, :] .=
        cloud_state.cld_path_liq[end - 1, :]
    @views cloud_state.cld_path_ice[end, :] .=
        cloud_state.cld_path_ice[end - 1, :]
    @views cloud_state.cld_frac[end, :] .= cloud_state.cld_frac[end - 1, :]
    return nothing
end

_extend_boundary_aerosol!(::Nothing) = nothing
function _extend_boundary_aerosol!(aerosol_state::AerosolState)
    @views aerosol_state.aero_size[:, end, :] .=
        aerosol_state.aero_size[:, end - 1, :]
    @views aerosol_state.aero_mass[:, end, :] .=
        aerosol_state.aero_mass[:, end - 1, :]
    return nothing
end

#####
##### Clipping and derived concentrations
#####

# Water-vapor volume mixing ratio, whichever storage is used. `idx_h2o` is only
# needed for the (per-gas) `Vmr` storage; it is ignored for `VmrGM`.
_vmr_h2o(vmr::Vmrs.VmrGM, idx_h2o) = vmr.vmr_h2o
_vmr_h2o(vmr::Vmrs.Vmr, idx_h2o) = view(vmr.vmr, idx_h2o, :, :)

"""
    clip!(as, p_min[, idx_h2o])

Clip the layer/level pressures of `as` to be at least `p_min`, and (for non-gray
states) clip the water-vapor volume mixing ratio to be nonnegative. Mutates `as`
in place and returns it.
"""
function clip!(as::GrayAtmosphericState, p_min, idx_h2o = nothing)
    p_lev = as.p_lev
    p_lay = getview_p_lay(as)
    @. p_lay = max(p_lay, p_min)
    @. p_lev = max(p_lev, p_min)
    return as
end
function clip!(as::AtmosphericState, p_min, idx_h2o = nothing)
    vmr_h2o = _vmr_h2o(as.vmr, idx_h2o)
    @. vmr_h2o = max(vmr_h2o, 0)
    p_lev = as.p_lev
    p_lay = getview_p_lay(as)
    @. p_lay = max(p_lay, p_min)
    @. p_lev = max(p_lev, p_min)
    return as
end

"""
    update_concentrations!(as, params, device[, idx_h2o])

Compute the column dry-air amount of `as` from its level pressures and
water-vapor volume mixing ratio. A no-op for gray radiation. Mutates `as` in
place and returns it.

This updates only the dry-air column amount (`compute_col_gas!`); it does **not**
recompute relative humidity. Relative humidity is a host responsibility: a caller that
needs an up-to-date `layer_relative_humidity` (e.g. for RH-dependent aerosol optics) must
call `compute_relative_humidity!` itself after updating temperature, pressure, or humidity.
"""
update_concentrations!(
    as::GrayAtmosphericState,
    params,
    device,
    idx_h2o = nothing,
) = as
function update_concentrations!(
    as::AtmosphericState,
    params,
    device,
    idx_h2o = nothing,
)
    Optics.compute_col_gas!(
        device,
        as.p_lev,
        getview_col_dry(as),
        params,
        _vmr_h2o(as.vmr, idx_h2o),
        as.lat,
    )
    return as
end
