using ClimaComms


"""
    RRTMGPGridParams(
        FT;
        context::ClimaComms.AbstractCommsContext,
        nlay::Int,
        ncol::Int
    )

A struct containing grid parameters.
"""
struct RRTMGPGridParams{FT, C}
    context::C
    nlay::Int
    ncol::Int
    extra_layer::Bool # for isothermal boundary layer
end
function RRTMGPGridParams(
    ::Type{FT};
    context::ClimaComms.AbstractCommsContext,
    nlay::Int,
    ncol::Int,
    extra_layer::Bool = false,
) where {FT}
    return RRTMGPGridParams{FT, typeof(context)}(context, nlay, ncol, extra_layer)
end

# Overload some methods for convenience:
Base.eltype(s::RRTMGPGridParams{FT}) where {FT} = FT
ClimaComms.device(s::RRTMGPGridParams) = ClimaComms.device(s.context)
ClimaComms.array_type(s::RRTMGPGridParams) = ClimaComms.array_type(ClimaComms.device(s))


"""
    AbstractRRTMGPMethod

An abstract type used for different radiation methods.

These subtypes are helpful for configuring lookup tables, and pre-configuring
caches for different radiation modes.
"""
abstract type AbstractRRTMGPMethod end

"""
    GrayRadiation

The gray radiation mode, given:
"""
struct GrayRadiation <: AbstractRRTMGPMethod end

"""
    ClearSkyRadiation

The clear-sky radiation mode, given:

 - `aerosol_radiation` bool indicating to turn on aerosol radiation
"""
struct ClearSkyRadiation <: AbstractRRTMGPMethod
    aerosol_radiation::Bool
end

"""
    AllSkyRadiation

The all-sky radiation mode, given:

 - `aerosol_radiation` bool indicating to turn on aerosol radiation
 - `reset_rng_seed` Reset the RNG seed before calling RRTMGP to a known value
   (the timestep number). When modeling cloud optics, RRTMGP uses a random
   number generator. Resetting the seed every time RRTMGP is called to a
   deterministic value ensures that the simulation is fully reproducible and
   can be restarted in a reproducible way. Disable this option when running
   production runs.
"""
struct AllSkyRadiation <: AbstractRRTMGPMethod
    aerosol_radiation::Bool
    reset_rng_seed::Bool
end

"""
    AllSkyRadiationWithClearSkyDiagnostics

The all-sky radiation mode (with diagnostics), given:

 - `aerosol_radiation` bool indicating to turn on aerosol radiation
 - `reset_rng_seed` Reset the RNG seed before calling RRTMGP to a known value
   (the timestep number). When modeling cloud optics, RRTMGP uses a random
   number generator. Resetting the seed every time RRTMGP is called to a
   deterministic value ensures that the simulation is fully reproducible and
   can be restarted in a reproducible way. Disable this option when running
   production runs.
"""
struct AllSkyRadiationWithClearSkyDiagnostics <: AbstractRRTMGPMethod
    aerosol_radiation::Bool
    reset_rng_seed::Bool
end
