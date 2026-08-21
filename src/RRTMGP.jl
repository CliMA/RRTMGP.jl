module RRTMGP
using Artifacts

_get_artifact_path() = joinpath(artifact"rrtmgp-data", "rrtmgp-data-1.9")

import ClimaComms

# File-naming convention: a CamelCase file defines the module of the same name;
# a snake_case file is an include-fragment of its enclosing module. The include
# list below is ordered by architectural layer (see the "Functional core" and
# "getter contract" docs pages).

# ── Numerical policy, physical parameters, data artifacts
include("Numerics.jl")
import .Numerics
import .Numerics: pow_fast
include("Parameters.jl")
import .Parameters as RP
include("ArtifactPaths.jl")

# ── Configuration types and the canonical aerosol name map
include(joinpath("api", "grid_params.jl"))
include(joinpath("api", "radiation_methods.jl"))
include(joinpath("api", "aerosols.jl"))

# ── Layer 0: core containers — states, sources, optics, fluxes, BCs, lookups
include(joinpath("optics", "VolumeMixingRatios.jl"))
include(joinpath("optics", "Fluxes.jl"))
include(joinpath("optics", "LookUpTables.jl"))
include(joinpath("optics", "AngularDiscretizations.jl"))
include(joinpath("optics", "AtmosphericStates.jl"))
include(joinpath("optics", "Sources.jl"))
include(joinpath("optics", "Optics.jl"))
include(joinpath("optics", "GrayAtmosphere.jl"))
include(joinpath("optics", "BCs.jl"))

# ── Layer 1: RTE solver workspaces and the functional solve_lw!/solve_sw! core
include(joinpath("rte", "RTE.jl"))
include(joinpath("rte", "RTESolver.jl"))

# ── Layer-1 preparation: level interpolation and grid adaptation
include(joinpath("api", "interpolation.jl"))
include(joinpath("api", "grid_adaptation.jl"))

# ── Layer 2: the RRTMGPSolver aggregate, orchestration, getters, validation
include(joinpath("api", "lookup_bundle.jl"))
include(joinpath("api", "solver.jl"))
include(joinpath("api", "update_fluxes.jl"))
include(joinpath("api", "getters.jl"))
include(joinpath("api", "validation.jl"))

# ── Layer 3: standalone front door
include(joinpath("api", "atmosphere_profile.jl"))
include(joinpath("api", "standalone.jl"))

# ── The public API surface (see test/public_api.jl)
include(joinpath("api", "public.jl"))

# ── Warm up the standalone gray solve at load time
include("precompile.jl")

end # module
