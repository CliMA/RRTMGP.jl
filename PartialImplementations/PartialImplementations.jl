using RRTMGP
using DocStringExtensions
using RRTMGP.OpticalProps
using RRTMGP.FortranIntrinsics
using RRTMGP.ArrayUtilities
using RRTMGP.GasOptics
using RRTMGP.GasConcentrations
using RRTMGP.RTESolver
using RRTMGP.Fluxes
using RRTMGP.SourceFunctions
using RRTMGP.CloudOptics

import RRTMGP.OpticalProps: validate!, subset_range!, delta_scale!, get_nmom

include(joinpath("..","ReadInputData","ReadInputs.jl"))
include(joinpath("..","ReadInputData","LoadCoefficients.jl"))
include(joinpath("..","ReadInputData","LoadCloudCoefficients.jl"))

include("OpticalPropsNStream.jl")
include("OpticalPropsNStream_kernels.jl")
include("ExtractSubset.jl")
