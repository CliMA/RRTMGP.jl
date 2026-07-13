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
solve. Those fluxes are exposed through the `clear_*` getters (e.g., `clear_net_flux`,
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
