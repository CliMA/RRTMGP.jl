#####
##### Layer-2 orchestration: update_{lw,sw,net,}fluxes! over an RRTMGPSolver
#####

"""
    update_lw_fluxes!(s::RRTMGPSolver)

Update the longwave fluxes.
"""
update_lw_fluxes!(s::RRTMGPSolver) = update_lw_fluxes!(s, _radiation_method(s))

update_lw_fluxes!(s::RRTMGPSolver, ::GrayRadiation) = RTESolver.solve_lw!(
    _longwave_solver(s),
    _atmospheric_state(s),
    _deep_atmosphere_inverse_scaling(s),
)
update_lw_fluxes!(s::RRTMGPSolver, ::ClearSkyRadiation) = RTESolver.solve_lw!(
    _longwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookup_lw,
    nothing,
    _lookup_tables(s).lookup_lw_aero,
    _deep_atmosphere_inverse_scaling(s),
)
update_lw_fluxes!(s::RRTMGPSolver, ::AllSkyRadiation) = RTESolver.solve_lw!(
    _longwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookup_lw,
    _lookup_tables(s).lookup_lw_cld,
    _lookup_tables(s).lookup_lw_aero,
    _deep_atmosphere_inverse_scaling(s),
)
function update_lw_fluxes!(
    s::RRTMGPSolver,
    ::AllSkyRadiationWithClearSkyDiagnostics,
)
    as = _atmospheric_state(s)
    lookups = _lookup_tables(s)
    lw_solver = _longwave_solver(s)
    ms = _deep_atmosphere_inverse_scaling(s)
    RTESolver.solve_lw!(
        lw_solver,
        as,
        lookups.lookup_lw,
        nothing,
        lookups.lookup_lw_aero,
        ms,
    )
    parent(clear_lw_flux_up(s)) .= parent(lw_flux_up(s))
    parent(clear_lw_flux_dn(s)) .= parent(lw_flux_dn(s))
    parent(clear_lw_flux(s)) .= parent(lw_flux_net(s))
    RTESolver.solve_lw!(
        lw_solver,
        as,
        lookups.lookup_lw,
        lookups.lookup_lw_cld,
        lookups.lookup_lw_aero,
        ms,
    )
end

"""
    update_sw_fluxes!(s::RRTMGPSolver)

Update the shortwave fluxes.
"""
update_sw_fluxes!(s::RRTMGPSolver) = update_sw_fluxes!(s, _radiation_method(s))

update_sw_fluxes!(s::RRTMGPSolver, ::GrayRadiation) = RTESolver.solve_sw!(
    _shortwave_solver(s),
    _atmospheric_state(s),
    _deep_atmosphere_inverse_scaling(s),
)
update_sw_fluxes!(s::RRTMGPSolver, ::ClearSkyRadiation) = RTESolver.solve_sw!(
    _shortwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookup_sw,
    nothing,
    _lookup_tables(s).lookup_sw_aero,
    _deep_atmosphere_inverse_scaling(s),
)
update_sw_fluxes!(s::RRTMGPSolver, ::AllSkyRadiation) = RTESolver.solve_sw!(
    _shortwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookup_sw,
    _lookup_tables(s).lookup_sw_cld,
    _lookup_tables(s).lookup_sw_aero,
    _deep_atmosphere_inverse_scaling(s),
)
function update_sw_fluxes!(
    s::RRTMGPSolver,
    ::AllSkyRadiationWithClearSkyDiagnostics,
)
    lookups = _lookup_tables(s)
    sw_solver = _shortwave_solver(s)
    as = _atmospheric_state(s)
    ms = _deep_atmosphere_inverse_scaling(s)
    RTESolver.solve_sw!(
        sw_solver,
        as,
        lookups.lookup_sw,
        nothing,
        lookups.lookup_sw_aero,
        ms,
    )
    parent(clear_sw_flux_up(s)) .= parent(sw_flux_up(s))
    parent(clear_sw_flux_dn(s)) .= parent(sw_flux_dn(s))
    parent(clear_sw_direct_flux_dn(s)) .= parent(sw_direct_flux_dn(s))
    parent(clear_sw_flux(s)) .= parent(sw_flux_net(s))

    RTESolver.solve_sw!(
        sw_solver,
        as,
        lookups.lookup_sw,
        lookups.lookup_sw_cld,
        lookups.lookup_sw_aero,
        ms,
    )
end

#####
##### Orchestration
#####

_deep_atmosphere_inverse_scaling(s::RRTMGPSolver) =
    s.deep_atmosphere_inverse_scaling

# Longwave lookup table (used for p_min) and H2O index (used for col_dry and
# clipping). Both are `nothing` for gray radiation, which uses no lookup tables.
_lw_lookup(s::RRTMGPSolver) = _lw_lookup(s, _radiation_method(s))
_lw_lookup(::RRTMGPSolver, ::GrayRadiation) = nothing
_lw_lookup(s::RRTMGPSolver, ::AbstractRRTMGPMethod) =
    _lookup_tables(s).lookup_lw

_idx_h2o(s::RRTMGPSolver) = _idx_h2o(s, _radiation_method(s))
_idx_h2o(::RRTMGPSolver, ::GrayRadiation) = nothing
_idx_h2o(s::RRTMGPSolver, ::AbstractRRTMGPMethod) =
    _lookup_tables(s).lookup_lw.idx_h2o

_maybe_reset_rng_seed!(::AbstractRRTMGPMethod, seedval) = nothing
function _maybe_reset_rng_seed!(
    rm::Union{AllSkyRadiation, AllSkyRadiationWithClearSkyDiagnostics},
    seedval,
)
    rm.reset_rng_seed && !isnothing(seedval) && Random.seed!(seedval)
    return nothing
end

"""
    update_net_fluxes!(s::RRTMGPSolver)

Combine the longwave and shortwave net fluxes into `net_flux(s)` (and, for
`AllSkyRadiationWithClearSkyDiagnostics`, the clear-sky pair into
`clear_net_flux(s)`).
"""
update_net_fluxes!(s::RRTMGPSolver) =
    update_net_fluxes!(s, _radiation_method(s))
function update_net_fluxes!(s::RRTMGPSolver, ::AbstractRRTMGPMethod)
    # Operate on the raw (boundary-extended) buffers, not `parent(getter(s))`: the
    # getters allocate a `view` per call that would then just be stripped by `parent`.
    s.net_flux_buffer .= s.lws.flux.flux_net .+ s.sws.flux.flux_net
    _update_band_net_fluxes!(s.lws)
    _update_band_net_fluxes!(s.sws)
    return nothing
end
function update_net_fluxes!(
    s::RRTMGPSolver,
    ::AllSkyRadiationWithClearSkyDiagnostics,
)
    s.net_flux_buffer .= s.lws.flux.flux_net .+ s.sws.flux.flux_net
    s.clear_net_flux_buffer .=
        s.clear_flux_lw.flux_net .+ s.clear_flux_sw.flux_net
    _update_band_net_fluxes!(s.lws)
    _update_band_net_fluxes!(s.sws)
    return nothing
end

_update_band_net_fluxes!(ws) = _update_band_net_fluxes!(hasproperty(ws, :band_flux) ? ws.band_flux : nothing)
_update_band_net_fluxes!(::Nothing) = nothing
function _update_band_net_fluxes!(band::Fluxes.FluxBand)
    band.flux_net .= band.flux_up .- band.flux_dn
    return nothing
end

"""
    update_fluxes!(s::RRTMGPSolver, seedval = nothing)

Run the radiation update: prepare the atmospheric state (interpolate levels,
add the isothermal boundary layer, clip pressures/temperatures/humidity to the
range the optics support, and compute concentrations), solve the
longwave and shortwave problems (applying `deep_atmosphere_inverse_scaling` if present), and
combine them into the net flux. Mutates `s` in place (its atmospheric state and
flux buffers) and returns `nothing` (read results via `net_flux(s)` and the
other flux getters). When the radiation method requests reproducible seeding,
`seedval` reseeds the RNG used for cloud sampling.

This is designed to be allocation-free and type-stable, which matters because a
host calls it every radiation step. CI asserts `@allocated == 0` and
`JET.@test_opt` for the gray Layer-2 aggregate and for the Layer-1
`solve_lw!`/`solve_sw!` kernels of the spectral modes (single-threaded CPU).

See also `prepare_atmosphere!`, `update_lw_fluxes!`, `update_sw_fluxes!`, and
`update_net_fluxes!`.
"""
function update_fluxes!(s::RRTMGPSolver, seedval = nothing)
    # opt-in input validation (see `check_values`/`validate_inputs`); a single
    # branch when off, so the zero-allocation contract is unaffected
    check_values[] && validate_inputs(s)
    _maybe_reset_rng_seed!(_radiation_method(s), seedval)
    prepare_atmosphere!(s)
    update_lw_fluxes!(s)
    update_sw_fluxes!(s)
    update_net_fluxes!(s)
    return nothing
end

"""
    prepare_atmosphere!(s::RRTMGPSolver)

Run the atmospheric-state preparation cascade without solving: interpolate
level pressures/temperatures from layer values (per the solver's
`interpolation`/`bottom_extrapolation` configuration; a no-op for
`NoInterpolation`), fill the isothermal boundary layer (when configured), clip
unphysical inputs (pressures below the lookup tables' minimum, negative water
vapor, and — for non-gray optics — temperatures outside the lookup tables'
range), and compute the dry-air column amounts. Mutates the solver's
atmospheric state in place and returns `nothing`.

[`update_fluxes!`](@ref) calls this before solving; call it to inspect
the prepared state (interpolated level values, column amounts) without running
the radiative transfer (the prepare/solve split familiar from other radiation
drivers). Relative humidity is managed by the host (see [`update_concentrations!`](@ref)).
"""
function prepare_atmosphere!(s::RRTMGPSolver)
    as = _atmospheric_state(s)
    lookup_lw = _lw_lookup(s)
    p_min = get_p_min(as, lookup_lw)
    interpolate_levels!(
        as,
        s.interpolation,
        s.bottom_extrapolation,
        s.params;
        center_z = s.center_z,
        face_z = s.face_z,
        isothermal_boundary_layer = s.grid_params.isothermal_boundary_layer,
    )
    s.grid_params.isothermal_boundary_layer &&
        add_isothermal_boundary_layer!(as, p_min)
    clip!(
        as,
        p_min,
        _idx_h2o(s);
        t_min = get_t_min(as, lookup_lw),
        t_max = get_t_max(as, lookup_lw),
    )
    update_concentrations!(
        as,
        s.params,
        ClimaComms.device(s.grid_params),
        _idx_h2o(s),
    )
    return nothing
end
