module RRTMGP

include("misc.jl")
include("fortran_intrinsics.jl")
include("Utilities.jl")

include(joinpath("rte","mo_util_array.jl"))
include(joinpath("rte","mo_optical_props.jl"))
include(joinpath("rte","mo_fluxes.jl"))
include(joinpath("rte","mo_source_functions.jl"))
include(joinpath("extensions","mo_fluxes_byband_kernels.jl"))
include(joinpath("extensions", "cloud_optics", "mo_cloud_optics.jl"))
# include(joinpath("extensions","mo_fluxes_byband.jl"))
include(joinpath("rte","RTESolver.jl"))

include(joinpath("rrtmgp","mo_rrtmgp_constants.jl"))
include(joinpath("rrtmgp","mo_gas_concentrations.jl"))
include(joinpath("rrtmgp","mo_gas_optics_rrtmgp.jl"))
include(joinpath("examples","mo_load_coefficients.jl"))
include(joinpath("examples","mo_simple_netcdf.jl"))
include(joinpath("examples","rfmip-clear-sky","mo_rfmip_io.jl"))
include(joinpath("examples","mo_load_cloud_coefficients.jl"))

end # module
