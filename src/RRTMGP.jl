module RRTMGP

include("Device.jl")

include(joinpath("optics", "GrayAngularDiscretizations.jl"))
include(joinpath("optics", "GrayOptics.jl"))
include(joinpath("optics", "GraySources.jl"))
include(joinpath("optics", "GrayFluxes.jl"))
include(joinpath("optics", "GrayAtmosphericStates.jl"))
include(joinpath("optics", "GrayUtils.jl"))
include(joinpath("optics", "GrayBCs.jl"))
include(joinpath("optics", "GrayAtmos.jl"))

include(joinpath("rte", "GrayRTESolver.jl"))

end # module
