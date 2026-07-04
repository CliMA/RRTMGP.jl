using ..Fluxes
using ..Sources
using ..AtmosphericStates: AtmosphericState, AbstractAtmosphericState
using ..RTE
using ..RTESolver
using ..BCs
using ..Optics
using ClimaComms
import Adapt
import Random

"""
    lookup_tables(grid_params::RRTMGPGridParams, radiation_method::AbstractRRTMGPMethod)

Build the lookup tables for `radiation_method`, returning a `NamedTuple`
`(; lookups, lu_kwargs)`:
 - `lookups`: the RRTMGP lookup tables and gas/aerosol index maps (varies by
   radiation mode; empty for `GrayRadiation`),
 - `lu_kwargs`: band/gas counts (`nbnd_lw`, `nbnd_sw`, and — spectral modes
   only — `ngas_lw`, `ngas_sw`).

Treat the return value as an **opaque token**: build it once and pass it back to
the `RRTMGPSolver` constructor via `lookups = ...` to avoid a second NetCDF read.
Its internal layout is not API and may become a typed struct in a future release.

The spectral (non-gray) methods are provided by an extension: load NCDatasets
(`using NCDatasets`) first.

TODO:
 - We should add type annotations for the data read from NC files as this will
   improve inference and the return type of `lookup_tables`.
"""
function lookup_tables end

# Gray radiation uses no lookup tables, so this method needs no NCDatasets (the
# spectral methods are provided by the NCDatasets extension).
lookup_tables(::RRTMGPGridParams, ::GrayRadiation) =
    (; lookups = (;), lu_kwargs = (; nbnd_lw = 1, nbnd_sw = 1))

# Non-gray radiation needs lookup tables loaded from NetCDF; surface a clear,
# actionable error if the NCDatasets extension has not been loaded.
#
# FUTURE IMPROVEMENT:
# To improve the standalone experience, RRTMGP.jl should eventually provide serialized Julia-native 
# lookup tables (e.g., via JLD2 or Serialization binary arrays in Artifacts). This would allow loading 
# real-gas lookup tables completely standalone without requiring NCDatasets.jl or C-library wrappers.
_check_lookup_support(::GrayRadiation) = nothing
function _check_lookup_support(radiation_method::AbstractRRTMGPMethod)
    if isnothing(Base.get_extension(@__MODULE__, :RRTMGPNCDatasetsExt))
        error(
            "RRTMGP's $(nameof(typeof(radiation_method))) mode requires lookup \
             tables; run `using NCDatasets` before constructing the solver.",
        )
    end
    return nothing
end

"""
    RRTMGPSolver

Aggregate bundling everything needed to run an RRTMGP radiation calculation. It
holds the radiation configuration, the atmospheric state, the longwave and
shortwave solvers, the lookup tables, and the output flux buffers, and exposes
getter methods (e.g. `layer_temperature`, `net_flux`) to read and write its data.
Construct it with the `RRTMGPSolver` constructor and drive it with
`update_fluxes!`.

# Fields
- `grid_params`: grid and device configuration (`RRTMGPGridParams`).
- `radiation_method`: the radiation method (gray, clear-sky, or all-sky).
- `interpolation`: scheme for filling level values from layer values.
- `bottom_extrapolation`: scheme for the bottom-level value.
- `params`: RRTMGP physical parameters.
- `sws`: shortwave RTE solver and its flux/scratch buffers.
- `lws`: longwave RTE solver and its flux/scratch buffers.
- `as`: the atmospheric state (solver inputs).
- `lookups`: the `(; lookups, lu_kwargs)` bundle from `lookup_tables` (for gray
  radiation the tables are empty and only the band counts are carried).
- `clear_flux_lw`: clear-sky longwave fluxes, or `nothing`.
- `clear_flux_sw`: clear-sky shortwave fluxes, or `nothing`.
- `center_z`: layer-center altitudes [m], or `nothing`.
- `face_z`: level (face) altitudes [m], or `nothing`.
- `deep_atmosphere_inverse_scaling`: `(nlev, ncol)` factor multiplied directly into the
  fluxes for deep-atmosphere geometric scaling, or `nothing`.
- `net_flux_buffer`: combined longwave + shortwave net flux at each level [W/m²], the full
  boundary-extended `(nlev, ncol)` buffer (read the domain-masked view via `net_flux(s)`).
- `clear_net_flux_buffer`: combined clear-sky net-flux buffer, or `nothing`.

# Constructor
    RRTMGPSolver(grid_params, radiation_method, params, bcs_lw, bcs_sw, as; <keyword arguments>)

Build a solver from the grid parameters, radiation method, physical `params`, the longwave
and shortwave boundary conditions, and the atmospheric state. Keyword arguments:
- `op_lw`, `op_sw`: longwave/shortwave optics (`OneScalar` or `TwoStream`); default `TwoStream`.
- `center_z`, `face_z`: layer-center and level altitudes [m], needed only for z-based
  interpolation; default `nothing`.
- `interpolation`: scheme for filling level values from layer values; default `NoInterpolation`.
- `bottom_extrapolation`: scheme for the bottom-level value; default `SameAsInterpolation`.
- `deep_atmosphere_inverse_scaling`: a `(nlev, ncol)` array multiplied directly into the
  fluxes for deep-atmosphere geometric scaling — the host supplies the multiplicative
  inverse of its metric scaling, hence the name — or `nothing` (default) for the
  shallow-atmosphere approximation.
- `lookups`: prebuilt lookup tables to reuse (avoids a second NetCDF read), or `nothing`
  (default) to build them internally.
- `spectral_fluxes`: if `true`, also retain per-band fluxes (two-stream, non-gray only);
  default `false`.
"""
struct RRTMGPSolver{
    S,
    RM,
    IN,
    BE,
    P,
    SWS,
    LWS,
    AS <: AbstractAtmosphericState,
    LU,
    CFLW,
    CFSW,
    ZC <: Union{AbstractArray, Nothing},
    ZF <: Union{AbstractArray, Nothing},
    MS <: Union{AbstractArray, Nothing},
    NF,
    CNF,
}
    grid_params::S
    radiation_method::RM
    interpolation::IN
    bottom_extrapolation::BE
    params::P
    sws::SWS
    lws::LWS
    as::AS
    lookups::LU
    clear_flux_lw::CFLW
    clear_flux_sw::CFSW
    center_z::ZC
    face_z::ZF
    deep_atmosphere_inverse_scaling::MS
    net_flux_buffer::NF
    clear_net_flux_buffer::CNF
end
Adapt.@adapt_structure RRTMGPSolver

function RRTMGPSolver(
    grid_params::RRTMGPGridParams,
    radiation_method::AbstractRRTMGPMethod,
    params::RP.ARP,
    bcs_lw::BCs.LwBCs,
    bcs_sw::SwBCs,
    as::AbstractAtmosphericState;
    op_lw::Optics.AbstractOpticalProps = Optics.TwoStream(grid_params),
    op_sw::Optics.AbstractOpticalProps = Optics.TwoStream(grid_params),
    center_z = nothing,
    face_z = nothing,
    interpolation::AbstractInterpolation = NoInterpolation(),
    bottom_extrapolation::AbstractBottomExtrapolation = SameAsInterpolation(),
    deep_atmosphere_inverse_scaling = nothing,
    lookups = nothing,
    spectral_fluxes = false,
)
    (; context) = grid_params

    # Non-scattering shortwave optics are only meaningful for gray radiation;
    # spectral shortwave radiation requires scattering (the spectral `NoScatSWRTE`
    # solve is gas-only and cannot be driven by `update_sw_fluxes!`).
    if op_sw isa OneScalar && !(radiation_method isa GrayRadiation)
        error(
            "non-scattering shortwave optics (`op_sw::OneScalar`) are only supported \
             with `GrayRadiation`; spectral shortwave radiation requires scattering, \
             so pass `op_sw = TwoStream(grid_params)`.",
        )
    end

    # z-based interpolation/extrapolation reads `center_z`/`face_z` at solve time;
    # fail here with an actionable message rather than a `MethodError` mid-solve.
    if (requires_z(interpolation) || requires_z(bottom_extrapolation)) &&
       (isnothing(center_z) || isnothing(face_z))
        error(
            "`$(nameof(typeof(interpolation)))` interpolation / \
             `$(nameof(typeof(bottom_extrapolation)))` extrapolation requires \
             altitudes; pass `center_z` and `face_z` to the `RRTMGPSolver` \
             constructor.",
        )
    end

    # Build (or reuse) the lookup tables first: the band counts are needed for the
    # optional spectral buffers below, and a host that already built them can pass them
    # in (`lookups = ...`) to avoid a second, expensive NetCDF read.
    if isnothing(lookups)
        _check_lookup_support(radiation_method)
        lookups = lookup_tables(grid_params, radiation_method)
    end

    src_lw =
        op_lw isa OneScalar ? SourceLWNoScat(grid_params; params) :
        SourceLW2Str(grid_params; params)
    src_sw = op_sw isa OneScalar ? nothing : SourceSW2Str(grid_params)

    fluxb_lw =
        radiation_method isa GrayRadiation ? nothing :
        Fluxes.FluxLW(grid_params)
    fluxb_sw =
        radiation_method isa GrayRadiation ? nothing :
        Fluxes.FluxSW(grid_params)

    flux_lw = Fluxes.FluxLW(grid_params)
    flux_sw = Fluxes.FluxSW(grid_params)
    net_flux_buffer = similar(flux_lw.flux_net)
    if radiation_method isa AllSkyRadiationWithClearSkyDiagnostics
        clear_flux_lw = Fluxes.FluxLW(grid_params)
        clear_flux_sw = Fluxes.FluxSW(grid_params)
        clear_net_flux_buffer = similar(flux_lw.flux_net)
    else
        clear_flux_lw = nothing
        clear_flux_sw = nothing
        clear_net_flux_buffer = nothing
    end

    # Optional per-band (spectrally-resolved) flux buffers, allocated only on request.
    # Only meaningful for spectral radiation with two-stream optics (gray radiation is a
    # single band, and the non-scattering solvers are not wired for per-band retention).
    if spectral_fluxes
        radiation_method isa GrayRadiation && error(
            "spectral_fluxes = true is not supported for GrayRadiation (a single band).",
        )
        (op_lw isa OneScalar || op_sw isa OneScalar) && error(
            "spectral_fluxes = true requires two-stream optics for both bands.",
        )
        band_flux_lw = Fluxes.FluxBand(grid_params, lookups.lu_kwargs.nbnd_lw)
        band_flux_sw = Fluxes.FluxBand(grid_params, lookups.lu_kwargs.nbnd_sw)
    else
        band_flux_lw = nothing
        band_flux_sw = nothing
    end

    sws =
        op_sw isa OneScalar ?
        RTE.NoScatSWRTE(context, op_sw, bcs_sw, fluxb_sw, flux_sw) :
        RTE.TwoStreamSWRTE(
            context,
            op_sw,
            src_sw, # inferrable
            bcs_sw, # has data
            fluxb_sw, # inferable, need radiation_method
            flux_sw, # views attached
            band_flux_sw, # optional per-band fluxes (or `nothing`)
        )
    lws =
        op_lw isa OneScalar ?
        RTE.NoScatLWRTE(
            context,
            op_lw,
            src_lw, # inferrable
            bcs_lw, # has data
            fluxb_lw, # inferable, need radiation_method
            flux_lw, # views attached
            AngularDiscretizations.AngularDiscretization(grid_params, 1),
        ) :
        RTE.TwoStreamLWRTE(
            context,
            op_lw,
            src_lw, # inferrable
            bcs_lw, # has data
            fluxb_lw, # inferable, need radiation_method
            flux_lw, # views attached
            band_flux_lw, # optional per-band fluxes (or `nothing`)
        )
    return RRTMGPSolver(
        grid_params,
        radiation_method,
        interpolation,
        bottom_extrapolation,
        params,
        sws,
        lws,
        as,
        lookups,
        clear_flux_lw,
        clear_flux_sw,
        center_z,
        face_z,
        deep_atmosphere_inverse_scaling,
        net_flux_buffer,
        clear_net_flux_buffer,
    )
end

"""
    update_lw_fluxes!(s::RRTMGPSolver)

Updates the longwave fluxes.
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
    _lookup_tables(s).lookups.lookup_lw,
    nothing,
    _lookup_tables(s).lookups.lookup_lw_aero,
    _deep_atmosphere_inverse_scaling(s),
)
update_lw_fluxes!(s::RRTMGPSolver, ::AllSkyRadiation) = RTESolver.solve_lw!(
    _longwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookups.lookup_lw,
    _lookup_tables(s).lookups.lookup_lw_cld,
    _lookup_tables(s).lookups.lookup_lw_aero,
    _deep_atmosphere_inverse_scaling(s),
)
function update_lw_fluxes!(
    s::RRTMGPSolver,
    ::AllSkyRadiationWithClearSkyDiagnostics,
)
    as = _atmospheric_state(s)
    (; lookups) = _lookup_tables(s)
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

Updates the shortwave fluxes.
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
    _lookup_tables(s).lookups.lookup_sw,
    nothing,
    _lookup_tables(s).lookups.lookup_sw_aero,
    _deep_atmosphere_inverse_scaling(s),
)
update_sw_fluxes!(s::RRTMGPSolver, ::AllSkyRadiation) = RTESolver.solve_sw!(
    _shortwave_solver(s),
    _atmospheric_state(s),
    _lookup_tables(s).lookups.lookup_sw,
    _lookup_tables(s).lookups.lookup_sw_cld,
    _lookup_tables(s).lookups.lookup_sw_aero,
    _deep_atmosphere_inverse_scaling(s),
)
function update_sw_fluxes!(
    s::RRTMGPSolver,
    ::AllSkyRadiationWithClearSkyDiagnostics,
)
    (; lookups) = _lookup_tables(s)
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
    _lookup_tables(s).lookups.lookup_lw

_idx_h2o(s::RRTMGPSolver) = _idx_h2o(s, _radiation_method(s))
_idx_h2o(::RRTMGPSolver, ::GrayRadiation) = nothing
_idx_h2o(s::RRTMGPSolver, ::AbstractRRTMGPMethod) =
    _lookup_tables(s).lookups.lookup_lw.idx_h2o

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
    return nothing
end
function update_net_fluxes!(
    s::RRTMGPSolver,
    ::AllSkyRadiationWithClearSkyDiagnostics,
)
    s.net_flux_buffer .= s.lws.flux.flux_net .+ s.sws.flux.flux_net
    s.clear_net_flux_buffer .=
        s.clear_flux_lw.flux_net .+ s.clear_flux_sw.flux_net
    return nothing
end

"""
    update_fluxes!(s::RRTMGPSolver, seedval = nothing)

Run the full radiation update: prepare the atmospheric state (interpolate levels,
add the isothermal boundary layer, clip, and compute concentrations), solve the
longwave and shortwave problems (applying `deep_atmosphere_inverse_scaling` if present), and
combine them into the net flux. Mutates `s` in place — its atmospheric state and
flux buffers — and returns `nothing` (read results via `net_flux(s)` and the
other flux getters). When the radiation method requests reproducible seeding,
`seedval` reseeds the RNG used for cloud sampling.

This is designed to be allocation-free and type-stable, which matters because a
host calls it every radiation step. CI asserts `@allocated == 0` and
`JET.@test_opt` for the gray Layer-2 aggregate and for the Layer-1
`solve_lw!`/`solve_sw!` kernels of the spectral modes (single-threaded CPU).

See also `update_lw_fluxes!`, `update_sw_fluxes!`, and `update_net_fluxes!`.
"""
function update_fluxes!(s::RRTMGPSolver, seedval = nothing)
    _maybe_reset_rng_seed!(_radiation_method(s), seedval)
    as = _atmospheric_state(s)
    p_min = get_p_min(as, _lw_lookup(s))
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
    clip!(as, p_min, _idx_h2o(s))
    update_concentrations!(
        as,
        s.params,
        ClimaComms.device(s.grid_params),
        _idx_h2o(s),
    )
    update_lw_fluxes!(s)
    update_sw_fluxes!(s)
    update_net_fluxes!(s)
    return nothing
end
