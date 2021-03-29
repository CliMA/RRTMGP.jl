using Test
using Pkg.Artifacts
using NCDatasets

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

include("reference_files.jl")
include("read_all_sky.jl")

I = Int
FT = Float64
DA = array_type()
max_threads = Int(256)
@show DA


lw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :lw) # lw lookup tables for gas optics
lw_cld_file = get_ref_filename(:lookup_tables, :cloudysky, λ = :lw) # lw cloud lookup tables
lw_input_file = get_ref_filename(:atmos_state, :cloudysky)      # all-sky atmos state
@show lw_file
@show lw_cld_file
@show lw_input_file

#reading longwave gas optics lookup data
ds_lw = Dataset(lw_file, "r")
lookup_lw, idx_gases = LookUpLW(ds_lw, I, FT, DA)
close(ds_lw)

# reading longwave cloud lookup data
ds_lw_cld = Dataset(lw_cld_file, "r")
lookup_lw_cld = LookUpCld(ds_lw_cld, I, FT, DA)
close(ds_lw_cld)

# reading input file 
ds_lw_in = Dataset(lw_input_file, "r")

setup_allsky_as(ds_lw_in, idx_gases, lookup_lw, FT, DA, max_threads)
