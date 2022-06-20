module RRTMGP

include("Device.jl")

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

end # module
