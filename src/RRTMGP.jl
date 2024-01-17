module RRTMGP

using Artifacts
lookup_data() = artifact"rrtmgp-lookup-data"

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
pow_fast(x, y) = exp(y * log(x))

import ClimaComms
include("Parameters.jl")
import .Parameters as RP

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

if !isdefined(Base, :get_extension)
	include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # module
