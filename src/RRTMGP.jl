module RRTMGP

include("Misc.jl")
include("FortranIntrinsics.jl")
include("Utilities.jl")

include(joinpath("rte","ArrayUtilities.jl"))
include(joinpath("rte","OpticalProps.jl"))
include(joinpath("rte","Fluxes.jl"))
include(joinpath("rte","SourceFunctions.jl"))
include(joinpath("extensions","CloudOptics.jl"))
include(joinpath("rte","RTESolver.jl"))

include(joinpath("rrtmgp","PhysicalConstants.jl"))
include(joinpath("rrtmgp","GasConcentrations.jl"))
include(joinpath("rrtmgp","GasOptics.jl"))
include(joinpath("examples","LoadCoefficients.jl"))
include(joinpath("examples","SimpleNetCDF.jl"))
include(joinpath("examples","RFMIPIO.jl"))
include(joinpath("examples","LoadCloudCoefficients.jl"))

end # module
