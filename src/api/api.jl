using ..Fluxes
using ..Sources
using ..AtmosphericStates: AtmosphericState
using ..RTE
using ..RTESolver
using ..BCs
using ..Optics
using ClimaComms

"""
    lookup_tables(::AbstractRRTMGPMethod, device, FloatType)

Return a NamedTuple, containing
 - RRTMGP lookup tables (varies by radiation mode)
 - nbnd_lw
 - nbnd_sw
 - ngas_lw (absent for GrayRadiation)
 - ngas_sw (absent for GrayRadiation)

This is an extension that requires NCDatasets be loaded prior to using RRTMGP.

TODO:
 - We should add type annotations for the data read from NC files as this will
   improve inference and the return type of `lookup_tables`.
"""
function lookup_tables end

"""
	RRTMGPSolver

A RRTMGP model that contains all RRTMGP structs,
and provides getter methods to retrieve data.
"""
struct RRTMGPSolver{S, RM, P, SWS, LWS, AS <: AtmosphericState, LU, CFLW, CFSW, V <: Union{AbstractArray, Nothing}}
    grid_params::S
    radiation_method::RM
    params::P
    sws::SWS
    lws::LWS
    as::AS
    lookups::LU
    clear_flux_lw::CFLW
    clear_flux_sw::CFSW
    center_z::V
    face_z::V
end

function RRTMGPSolver(
    grid_params::RRTMGPGridParams,
    radiation_method::AbstractRRTMGPMethod,
    params::RP.ARP,
    bcs_lw::BCs.LwBCs,
    bcs_sw::SwBCs,
    as::AtmosphericState;
    op_lw::Optics.AbstractOpticalProps = Optics.TwoStream(grid_params),
    op_sw::Optics.AbstractOpticalProps = Optics.TwoStream(grid_params),
    center_z = nothing,
    face_z = nothing,
)
    (; context) = grid_params

    src_lw = op_lw isa OneScalar ? SourceLWNoScat(grid_params; params) : SourceLW2Str(grid_params; params)
    src_sw = op_sw isa OneScalar ? nothing : SourceSW2Str(grid_params)

    fluxb_lw = radiation_method isa GrayRadiation ? nothing : Fluxes.FluxLW(grid_params)
    fluxb_sw = radiation_method isa GrayRadiation ? nothing : Fluxes.FluxSW(grid_params)

    flux_lw = Fluxes.FluxLW(grid_params)
    flux_sw = Fluxes.FluxSW(grid_params)
    if radiation_method isa AllSkyRadiationWithClearSkyDiagnostics
        clear_flux_lw = Fluxes.FluxLW(grid_params)
        clear_flux_sw = Fluxes.FluxSW(grid_params)
    else
        clear_flux_lw = nothing
        clear_flux_sw = nothing
    end

    sws = RTE.TwoStreamSWRTE(
        context,
        op_sw,
        src_sw, # inferrable
        bcs_sw, # has data
        fluxb_sw, # inferable, need radiation_method
        flux_sw, # views attached
    )
    lws = RTE.TwoStreamLWRTE(
        context,
        op_lw,
        src_lw, # inferrable
        bcs_lw, # has data
        fluxb_lw, # inferable, need radiation_method
        flux_lw, # views attached
    )
    lookups = lookup_tables(grid_params, radiation_method)
    return RRTMGPSolver(
        grid_params,
        radiation_method,
        params,
        sws,
        lws,
        as,
        lookups,
        clear_flux_lw,
        clear_flux_sw,
        center_z,
        face_z,
    )
end

"""
    update_lw_fluxes!(s::RRTMGPSolver)

Updates the longwave fluxes.
"""
update_lw_fluxes!(s::RRTMGPSolver) = update_lw_fluxes!(s, _radiation_method(s))

update_lw_fluxes!(s::RRTMGPSolver, ::GrayRadiation) = RTESolver.solve_lw!(_longwave_solver(s), _atmospheric_state(s))
update_lw_fluxes!(s::RRTMGPSolver, ::ClearSkyRadiation) = RTESolver.solve_lw!(
    _longwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookups.lookup_lw,
    nothing,
    _lookup_tables(s).lookups.lookup_lw_aero,
)
update_lw_fluxes!(s::RRTMGPSolver, ::AllSkyRadiation) = RTESolver.solve_lw!(
    _longwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookups.lookup_lw,
    _lookup_tables(s).lookups.lookup_lw_cld,
    _lookup_tables(s).lookups.lookup_lw_aero,
)
function update_lw_fluxes!(s::RRTMGPSolver, ::AllSkyRadiationWithClearSkyDiagnostics)
    as = _atmospheric_state(s)
    (; lookups) = _lookup_tables(s)
    lw_solver = _longwave_solver(s)
    RTESolver.solve_lw!(lw_solver, as, lookups.lookup_lw, nothing, lookups.lookup_lw_aero)
    parent(clear_lw_flux_up(s)) .= parent(lw_flux_up(s))
    parent(clear_lw_flux_dn(s)) .= parent(lw_flux_dn(s))
    parent(clear_lw_flux(s)) .= parent(lw_flux_net(s))
    RTESolver.solve_lw!(lw_solver, as, lookups.lookup_lw, lookups.lookup_lw_cld, lookups.lookup_lw_aero)
end

"""
    update_sw_fluxes!(s::RRTMGPSolver)

Updates the shortwave fluxes.
"""
update_sw_fluxes!(s::RRTMGPSolver) = update_sw_fluxes!(s, _radiation_method(s))

update_sw_fluxes!(s::RRTMGPSolver, ::GrayRadiation) = RTESolver.solve_sw!(_shortwave_solver(s), _atmospheric_state(s))
update_sw_fluxes!(s::RRTMGPSolver, ::ClearSkyRadiation) = RTESolver.solve_sw!(
    _shortwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookups.lookup_sw,
    nothing,
    _lookup_tables(s).lookups.lookup_sw_aero,
)
update_sw_fluxes!(s::RRTMGPSolver, ::AllSkyRadiation) = RTESolver.solve_sw!(
    _shortwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookups.lookup_sw,
    _lookup_tables(s).lookups.lookup_sw_cld,
    _lookup_tables(s).lookups.lookup_sw_aero,
)
function update_sw_fluxes!(s::RRTMGPSolver, ::AllSkyRadiationWithClearSkyDiagnostics)
    (; lookups) = _lookup_tables(s)
    sw_solver = _shortwave_solver(s)
    as = _atmospheric_state(s)
    RTESolver.solve_sw!(sw_solver, as, lookups.lookup_sw, nothing, lookups.lookup_sw_aero)
    parent(clear_sw_flux_up(s)) .= parent(sw_flux_up(s))
    parent(clear_sw_flux_dn(s)) .= parent(sw_flux_dn(s))
    parent(clear_sw_direct_flux_dn(s)) .= parent(sw_direct_flux_dn(s))
    parent(clear_sw_flux(s)) .= parent(sw_flux_net(s))

    RTESolver.solve_sw!(sw_solver, as, lookups.lookup_sw, lookups.lookup_sw_cld, lookups.lookup_sw_aero)
end
