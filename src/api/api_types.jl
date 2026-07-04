using ClimaComms


"""
    RRTMGPGridParams(
        FT;
        context::ClimaComms.AbstractCommsContext,
        domain_nlay::Int,
        ncol::Int,
        isothermal_boundary_layer::Bool = false,
    )

Grid parameters for RRTMGP, parametrized on the float type `FT`.

Specify `domain_nlay`, the number of layers in your **physical** grid. When
`isothermal_boundary_layer = true`, RRTMGP adds one extra layer at the top of the
domain *internally*, so the stored field `nlay` — the total number of layers RRTMGP
works with — is `domain_nlay + 1` (and just `domain_nlay` when the flag is `false`).
Callers never add the extra layer themselves, so the boundary-layer bookkeeping
cannot be got wrong; the getters return domain-sized arrays and `heating_rate`
reports on the physical grid.

# Keyword Arguments
 - `context`: the `ClimaComms` context.
 - `domain_nlay`: the number of layers in the physical domain.
 - `ncol`: the number of columns.
 - `isothermal_boundary_layer = false`: whether RRTMGP adds an isothermal boundary
    layer/level at the top of the domain.

Migrating from the pre-`RRTMGPGridParams` workspace constructors (`OneScalar(FT, ncol,
nlay, DA)`, `TwoStream(...)`, ...): pass the same physical layer count as `domain_nlay`,
and do not add the boundary layer yourself — RRTMGP does when `isothermal_boundary_layer =
true`. A caller that previously baked the extra layer into `nlay` must drop it here, or the
boundary layer is double-counted.
"""
struct RRTMGPGridParams{FT, C}
    context::C
    "total layer count RRTMGP works with = `domain_nlay + Int(isothermal_boundary_layer)`"
    nlay::Int
    ncol::Int
    isothermal_boundary_layer::Bool # for isothermal boundary layer
end
function RRTMGPGridParams(
    ::Type{FT};
    context::ClimaComms.AbstractCommsContext,
    domain_nlay::Int,
    ncol::Int,
    isothermal_boundary_layer::Bool = false,
) where {FT}
    # Add the isothermal boundary layer here so callers only ever pass the physical
    # `domain_nlay`; the stored `nlay` is the total and cannot be got wrong.
    nlay = domain_nlay + Int(isothermal_boundary_layer)
    return RRTMGPGridParams{FT, typeof(context)}(
        context,
        nlay,
        ncol,
        isothermal_boundary_layer,
    )
end

# Overload some methods for convenience:
Base.eltype(s::RRTMGPGridParams{FT}) where {FT} = FT
ClimaComms.device(s::RRTMGPGridParams) = ClimaComms.device(s.context)
ClimaComms.array_type(s::RRTMGPGridParams) =
    ClimaComms.array_type(ClimaComms.device(s))


"""
    AbstractRRTMGPMethod

An abstract type used for different radiation methods.

These subtypes are helpful for configuring lookup tables, and pre-configuring
caches for different radiation modes.
"""
abstract type AbstractRRTMGPMethod end

"""
    GrayRadiation()

Gray-gas radiation: a single-band ("gray") atmosphere whose optical thickness follows a
prescribed analytic profile (a `GrayOpticalThickness*` parameter set) instead of
correlated-k lookup tables. Because it needs no NetCDF data, this is the standalone/teaching
mode driven by [`solve_gray`](@ref RRTMGP.solve_gray). Takes no options.
"""
struct GrayRadiation <: AbstractRRTMGPMethod end

"""
    ClearSkyRadiation(aerosol_radiation::Bool)

Clear-sky spectral radiation: molecular (gas) absorption and emission from the correlated-k
lookup tables, with **no clouds**. Requires the lookup tables (load `NCDatasets`).

# Fields
 - `aerosol_radiation::Bool`: include aerosol optics.
"""
struct ClearSkyRadiation <: AbstractRRTMGPMethod
    aerosol_radiation::Bool
end

"""
    AllSkyRadiation(aerosol_radiation::Bool, reset_rng_seed::Bool)

All-sky spectral radiation: gas absorption plus **cloud optics**, with cloud overlap sampled
by McICA. Requires the lookup tables (load `NCDatasets`).

# Fields
 - `aerosol_radiation::Bool`: include aerosol optics.
 - `reset_rng_seed::Bool`: when `true`, `update_fluxes!(s, seedval)` reseeds the RNG with
   `seedval` before the solve (hosts typically pass the timestep number); with no `seedval`
   the flag has no effect. Because the McICA cloud sampler draws random numbers, reseeding
   makes a CPU run fully reproducible and restartable; disable it for production runs. On
   the GPU the sampler draws from the device RNG, which this does not seed — per-column
   McICA sampling is not reproducible there (see `build_cloud_mask!`).
"""
struct AllSkyRadiation <: AbstractRRTMGPMethod
    aerosol_radiation::Bool
    reset_rng_seed::Bool
end

"""
    AllSkyRadiationWithClearSkyDiagnostics(aerosol_radiation::Bool, reset_rng_seed::Bool)

Like [`AllSkyRadiation`](@ref), but each call also runs a parallel **cloud-free** (clear-sky)
solve. Those fluxes are exposed through the `clear_*` getters (e.g. `clear_net_flux`,
`clear_lw_flux_up`); differenced against the all-sky fluxes they give the cloud radiative
effect. Requires the lookup tables (load `NCDatasets`).

# Fields
 - `aerosol_radiation::Bool`: include aerosol optics.
 - `reset_rng_seed::Bool`: reseed the RNG from the `seedval` passed to `update_fluxes!`
   (see [`AllSkyRadiation`](@ref)).
"""
struct AllSkyRadiationWithClearSkyDiagnostics <: AbstractRRTMGPMethod
    aerosol_radiation::Bool
    reset_rng_seed::Bool
end
