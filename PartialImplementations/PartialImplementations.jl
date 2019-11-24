using RRTMGP
using DocStringExtensions
using RRTMGP.OpticalProps
using RRTMGP.SimpleNetCDF
using RRTMGP.FortranIntrinsics
using RRTMGP.ArrayUtilities
using RRTMGP.GasOptics
using RRTMGP.GasConcentrations
using RRTMGP.RTESolver
using RRTMGP.Fluxes
using RRTMGP.LoadCoefficients
using RRTMGP.RFMIPIO
using RRTMGP.SourceFunctions
using RRTMGP.CloudOptics
using RRTMGP.LoadCloudCoefficients

import RRTMGP.OpticalProps: validate!, subset_range!, delta_scale!, get_nmom

include("OpticalPropsNStream.jl")
include("OpticalPropsNStream_kernels.jl")
include("ExtractSubset.jl")
