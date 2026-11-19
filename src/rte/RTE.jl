module RTE
using Adapt
import ClimaComms
using ..AngularDiscretizations
using ..AtmosphericStates
using ..Sources
using ..Fluxes
using ..Optics
import ..RRTMGPGridParams
using ..BCs

import ..Parameters as RP

export NoScatLWRTE, TwoStreamLWRTE, NoScatSWRTE, TwoStreamSWRTE
export LWGptGroups, SWGptGroups, n_gpt_groups

# Device-dependent default for the workspaces' `state_cache` keyword. On the
# GPU, a `TransposedStateCache` gives the g-point loop coalesced state reads
# across the warp's columns. On the CPU, the kernels read the
# `AtmosphericState` directly (`nothing`): its vertical-first layout is
# already cache-friendly for the per-column sweeps, and a transposed copy
# would only add work.
_default_state_cache(gp::RRTMGPGridParams) =
    _default_state_cache(ClimaComms.device(gp), gp)
_default_state_cache(::ClimaComms.AbstractCPUDevice, gp::RRTMGPGridParams) =
    nothing
_default_state_cache(::ClimaComms.AbstractDevice, gp::RRTMGPGridParams) =
    TransposedStateCache(gp)
# With the state in hand, the cache also takes column-first copies of the gas,
# cloud and aerosol inputs it carries.
_default_state_cache(gp::RRTMGPGridParams, as::AtmosphericState) =
    _default_state_cache(ClimaComms.device(gp), gp, as)
_default_state_cache(
    ::ClimaComms.AbstractCPUDevice,
    gp::RRTMGPGridParams,
    as::AtmosphericState,
) = nothing
_default_state_cache(
    ::ClimaComms.AbstractDevice,
    gp::RRTMGPGridParams,
    as::AtmosphericState,
) = TransposedStateCache(gp, as)

# How many pieces the fused g-point loop is split into. One is the serial loop
# every other solver uses; above one the solve widens its workspaces to
# `ncol * ngroups` and covers them in a single kernel launch, trading that
# memory for blocks the one wave cannot otherwise fill.
const MAX_NGPT_GROUPS = 16

# Wide columns the default group count aims for: enough blocks to cover an A100
# many waves over, which is where the imbalance between columns stops showing,
# and past which another group adds no parallelism the scheduler can use. Every
# widened workspace is linear in `ngroups`, so this also bounds their memory --
# the accumulators are `ncol * ngroups * nlev` by two planes in the longwave and
# three in the shortwave.
const TARGET_GPT_GROUP_COLUMNS = 393216

"""
    default_ngpt_groups(grid_params::RRTMGPGridParams)

How many g-point groups a two-stream solver widens to when the caller names
none: enough to reach `TARGET_GPT_GROUP_COLUMNS` wide columns, capped at
`MAX_NGPT_GROUPS`.

A column count that already fills the device gets one group, and so does a CPU,
whose driver solves the whole g-point loop per column and ignores the workspaces.
Both gates matter because the count rises as columns fall, so an ungated default
would give the smallest runs the widest workspaces. Pass `ngpt_groups` to
override in either direction.

Grouping changes the order in which a column's g-points are summed, and
floating-point addition is not associative, so grouped fluxes differ from a
serial solve at roundoff -- about 1e-6 relative in the configurations tested --
and the clear-sky half of a both-skies solve no longer equals a standalone
clear-sky solve exactly. Pass `ngpt_groups = 1` where bit-for-bit agreement is
required. The cloud-cover tally is a count of g-points and stays exact. What
grouping does to a sum is the whole of the difference: at one group the g-points
are taken in the serial order, and the solve is then bit-for-bit with the serial
one, which is what the tests assert.
"""
default_ngpt_groups(grid_params::RRTMGPGridParams) =
    ClimaComms.device(grid_params) isa ClimaComms.AbstractCPUDevice ? 1 :
    clamp(cld(TARGET_GPT_GROUP_COLUMNS, grid_params.ncol), 1, MAX_NGPT_GROUPS)

include("gpt_groups.jl")

"""
    NoScatLWRTE(grid_params::RRTMGPGridParams; params, sfc_emis, inc_flux,
                n_gauss_angles = 1)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a non-scattering longwave simulation.

# Fields
- `context`: ClimaComms context.
- `op`: Optical properties.
- `src`: Longwave source functions.
- `bcs`: Longwave boundary conditions.
- `fluxb`: Temporary storage for bandwise calculations.
- `flux`: Longwave fluxes.
- `angle_disc`: Angular discretization; `n_gauss_angles` sets the number of
  Gauss-Jacobi-5 quadrature angles (1-4) used for the hemispheric integral of
  the spectral solve. One angle is the diffusivity approximation; more angles
  resolve the angular integral better, at proportional cost. Gray radiation is
  single-angle and rejects anything but `1`.
- `state_cache`: [`TransposedStateCache`](@ref
  RRTMGP.AtmosphericStates.TransposedStateCache) with column-first copies of
  the hot state arrays (refreshed at each spectral solve), or `nothing`.
"""
struct NoScatLWRTE{
    C,
    OP <: OneScalar,
    SL <: SourceLWNoScat,
    BC <: LwBCs,
    FXBL,
    FXL <: FluxLW,
    AD,
    TSC,
}
    context::C
    op::OP
    src::SL
    bcs::BC
    fluxb::FXBL
    flux::FXL
    angle_disc::AD
    state_cache::TSC
end
Adapt.@adapt_structure NoScatLWRTE

function NoScatLWRTE(
    grid_params::RRTMGPGridParams;
    params,
    sfc_emis,
    inc_flux,
    state_cache = _default_state_cache(grid_params),
    n_gauss_angles::Int = 1,
)
    (; context) = grid_params
    op = OneScalar(grid_params)
    src = SourceLWNoScat(grid_params; params)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(grid_params)
    flux = FluxLW(grid_params)
    ad = AngularDiscretization(grid_params, n_gauss_angles)
    return NoScatLWRTE(context, op, src, bcs, fluxb, flux, ad, state_cache)
end

"""
    TwoStreamLWRTE(grid_params::RRTMGPGridParams; params, sfc_emis, inc_flux)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a `2-stream` longwave simulation.

# Fields
- `context`: ClimaComms context.
- `op`: Optical properties.
- `src`: Longwave source functions.
- `bcs`: Longwave boundary conditions.
- `fluxb`: Temporary storage for bandwise calculations.
- `flux`: Longwave fluxes.
- `band_flux`: Optional per-band longwave fluxes `(nlev, ncol, n_bnd)`, or `nothing`.
- `state_cache`: [`TransposedStateCache`](@ref
  RRTMGP.AtmosphericStates.TransposedStateCache) with column-first copies of
  the hot state arrays (refreshed at each spectral solve), or `nothing`.
"""
struct TwoStreamLWRTE{
    C,
    OP <: TwoStream,
    SL <: SourceLW2Str,
    BC <: LwBCs,
    FXBL,
    FXL <: FluxLW,
    FXBND,
    TSC,
    GG,
}
    context::C
    op::OP
    src::SL
    bcs::BC
    fluxb::FXBL
    flux::FXL
    band_flux::FXBND
    state_cache::TSC
    gpt_groups::GG
end
Adapt.@adapt_structure TwoStreamLWRTE
# By default no per-band (spectrally-resolved) fluxes are retained, and the
# g-point loop is not split.
TwoStreamLWRTE(context, op, src, bcs, fluxb, flux, band_flux) =
    TwoStreamLWRTE(context, op, src, bcs, fluxb, flux, band_flux, nothing)
TwoStreamLWRTE(context, op, src, bcs, fluxb, flux, band_flux, state_cache) =
    TwoStreamLWRTE(
        context,
        op,
        src,
        bcs,
        fluxb,
        flux,
        band_flux,
        state_cache,
        nothing,
    )

function TwoStreamLWRTE(
    grid_params::RRTMGPGridParams;
    params,
    sfc_emis,
    inc_flux,
    state_cache = _default_state_cache(grid_params),
    ngpt_groups::Int = default_ngpt_groups(grid_params),
)
    (; context) = grid_params
    op = TwoStream(grid_params)
    src = SourceLW2Str(grid_params; params)
    bcs = LwBCs(sfc_emis, inc_flux)
    fluxb = FluxLW(grid_params)
    flux = FluxLW(grid_params)
    # Only the fused both-skies solve reads these; the others ignore them.
    gpt_groups =
        ngpt_groups > 1 ?
        LWGptGroups(grid_params, ngpt_groups; params) : nothing
    return TwoStreamLWRTE(
        context,
        op,
        src,
        bcs,
        fluxb,
        flux,
        nothing,
        state_cache,
        gpt_groups,
    )
end

"""
    NoScatSWRTE(grid_params::RRTMGPGridParams; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a non-scattering shortwave simulation.

# Fields
- `context`: ClimaComms context.
- `op`: Optical properties.
- `bcs`: Shortwave boundary conditions.
- `fluxb`: Temporary storage for bandwise calculations.
- `flux`: Shortwave fluxes.
- `state_cache`: [`TransposedStateCache`](@ref
  RRTMGP.AtmosphericStates.TransposedStateCache) with column-first copies of
  the hot state arrays (refreshed at each spectral solve), or `nothing`.
"""
struct NoScatSWRTE{C, OP <: OneScalar, BC <: SwBCs, FXBS, FXS <: FluxSW, TSC}
    context::C
    op::OP
    bcs::BC
    fluxb::FXBS
    flux::FXS
    state_cache::TSC
end
Adapt.@adapt_structure NoScatSWRTE

function NoScatSWRTE(
    grid_params::RRTMGPGridParams;
    cos_zenith,
    toa_flux,
    sfc_alb_direct,
    inc_flux_diffuse,
    sfc_alb_diffuse,
    state_cache = _default_state_cache(grid_params),
)
    (; context) = grid_params
    op = OneScalar(grid_params)
    bcs = SwBCs(
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
    fluxb = FluxSW(grid_params)
    flux = FluxSW(grid_params)
    return NoScatSWRTE(context, op, bcs, fluxb, flux, state_cache)
end

"""
    TwoStreamSWRTE(grid_params::RRTMGPGridParams; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

A high-level RRTMGP data structure storing the optical
properties, sources, boundary conditions and fluxes
configurations for a `2-stream` shortwave simulation.

# Fields
- `context`: ClimaComms context.
- `op`: Optical properties.
- `src`: Shortwave source functions.
- `bcs`: Shortwave boundary conditions.
- `fluxb`: Temporary storage for bandwise calculations.
- `flux`: Shortwave fluxes.
- `band_flux`: Optional per-band shortwave fluxes `(nlev, ncol, n_bnd)`, or `nothing`.
- `state_cache`: [`TransposedStateCache`](@ref
  RRTMGP.AtmosphericStates.TransposedStateCache) with column-first copies of
  the hot state arrays (refreshed at each spectral solve), or `nothing`.
"""
struct TwoStreamSWRTE{
    C,
    OP <: TwoStream,
    SS,
    BC <: SwBCs,
    FXBS,
    FXS <: FluxSW,
    FXBND,
    TSC,
    GG,
}
    context::C
    op::OP
    src::SS
    bcs::BC
    fluxb::FXBS
    flux::FXS
    band_flux::FXBND
    state_cache::TSC
    gpt_groups::GG
end
Adapt.@adapt_structure TwoStreamSWRTE
# By default no per-band (spectrally-resolved) fluxes are retained.
TwoStreamSWRTE(context, op, src, bcs, fluxb, flux, band_flux) =
    TwoStreamSWRTE(context, op, src, bcs, fluxb, flux, band_flux, nothing)
TwoStreamSWRTE(context, op, src, bcs, fluxb, flux, band_flux, state_cache) =
    TwoStreamSWRTE(
        context,
        op,
        src,
        bcs,
        fluxb,
        flux,
        band_flux,
        state_cache,
        nothing,
    )

function TwoStreamSWRTE(
    grid_params::RRTMGPGridParams;
    cos_zenith,
    toa_flux,
    sfc_alb_direct,
    inc_flux_diffuse,
    sfc_alb_diffuse,
    state_cache = _default_state_cache(grid_params),
    ngpt_groups::Int = default_ngpt_groups(grid_params),
)
    (; context) = grid_params
    op = TwoStream(grid_params)
    src = SourceSW2Str(grid_params)
    bcs = SwBCs(
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
    fluxb = FluxSW(grid_params)
    flux = FluxSW(grid_params)
    # Only the fused both-skies solve reads these; the others ignore them.
    gpt_groups =
        ngpt_groups > 1 ?
        SWGptGroups(grid_params, ngpt_groups) : nothing
    return TwoStreamSWRTE(
        context,
        op,
        src,
        bcs,
        fluxb,
        flux,
        nothing,
        state_cache,
        gpt_groups,
    )
end


end
