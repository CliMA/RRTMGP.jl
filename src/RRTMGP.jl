module RRTMGP
using Artifacts

get_artifact_path() = joinpath(artifact"rrtmgp-data", "rrtmgp-data-1.9")

import ClimaComms
include("Numerics.jl")
import .Numerics
import .Numerics: pow_fast
include("Parameters.jl")
import .Parameters as RP

include("ArtifactPaths.jl")
include(joinpath("api", "api_types.jl"))
include(joinpath("api", "aerosols.jl"))
include(joinpath("optics", "Vmrs.jl"))
include(joinpath("optics", "LookUpTables.jl"))
include(joinpath("optics", "AngularDiscretizations.jl"))
include(joinpath("optics", "AtmosphericStates.jl"))
include(joinpath("optics", "Sources.jl"))
include(joinpath("optics", "Optics.jl"))
include(joinpath("optics", "Fluxes.jl"))
include(joinpath("optics", "GrayUtils.jl"))
include(joinpath("optics", "BCs.jl"))
include(joinpath("optics", "RTE.jl"))

include(joinpath("rte", "RTESolver.jl"))

include(joinpath("api", "interpolation.jl"))
include(joinpath("api", "grid_adaptation.jl"))

include(joinpath("api", "api.jl"))
include(joinpath("api", "getters.jl"))
include(joinpath("api", "validation.jl"))
include(joinpath("api", "standalone.jl"))

end # module
