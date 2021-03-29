using Test
using Pkg.Artifacts
using RRTMGP
using RRTMGP.Device: array_type, array_device
using RRTMGP.Vmrs
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.AngularDiscretizations
using RRTMGP.RTE
using RRTMGP.RTESolver

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
# overriding some parameters to match with RRTMGP FORTRAN code
CLIMAParameters.Planet.grav(::EarthParameterSet) = 9.80665
CLIMAParameters.Planet.molmass_dryair(::EarthParameterSet) = 0.028964
CLIMAParameters.Planet.molmass_water(::EarthParameterSet) = 0.018016

using NCDatasets

I = Int
FT = Float64
DA = array_type()
@show DA

homefolder = "RRTMGPReferenceData"
#homefolder = artifact"RRTMGPReferenceData"
sw_cld_file = joinpath(
    homefolder,
    "lookup_tables",
    "cloud_optics",
    "rrtmgp-cloud-optics-coeffs-sw.nc",
)

# reading longwave lookup data
ds_sw_cld = Dataset(sw_cld_file, "r")

lookup_sw_cld = LookUpCld(ds_sw_cld, I, FT, DA)
#close(ds_sw_cloud)
