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

Build the lookup tables for `radiation_method`, returning a
[`LookupBundle`](@ref) — the gas/cloud/aerosol lookup tables, the name→index
maps, and the band/gas counts. Build it once and pass it back to the
`RRTMGPSolver` constructor via `lookups = ...` to avoid a second NetCDF read,
or cache it on disk with [`save_lookup_tables`](@ref).

The spectral (non-gray) methods are provided by an extension: load NCDatasets
(`using NCDatasets`) first (or [`load_lookup_tables`](@ref) from a cache).
"""
function lookup_tables end

# Gray radiation uses no lookup tables, so this method needs no NCDatasets (the
# spectral methods are provided by the NCDatasets extension).
lookup_tables(::RRTMGPGridParams, ::GrayRadiation) =
    LookupBundle(; nbnd_lw = 1, nbnd_sw = 1)

# Non-gray radiation needs lookup tables loaded from NetCDF; surface a clear,
# actionable error if the NCDatasets extension has not been loaded.
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
- `lookups`: the [`LookupBundle`](@ref) from `lookup_tables` (for gray
  radiation the tables are `nothing` and only the band counts are carried).
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
        band_flux_lw = Fluxes.FluxBand(grid_params, lookups.nbnd_lw)
        band_flux_sw = Fluxes.FluxBand(grid_params, lookups.nbnd_sw)
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

