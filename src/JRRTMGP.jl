module JRRTMGP

include("fortran_intrinsics.jl")
include(joinpath("rte","mo_util_array.jl"))
include(joinpath("rte","kernels","mo_optical_props_kernels.jl"))
include(joinpath("rte","kernels","mo_fluxes_broadband_kernels.jl"))
include(joinpath("rte","mo_optical_props.jl"))
include(joinpath("rte","mo_fluxes.jl"))
include(joinpath("rte","mo_source_functions.jl"))

include(joinpath("rrtmgp","mo_util_string.jl"))
include(joinpath("rrtmgp","mo_rrtmgp_constants.jl"))
include(joinpath("rrtmgp","kernels","mo_reorder_kernels.jl"))
include(joinpath("rrtmgp","mo_util_reorder.jl"))
# include(joinpath("rrtmgp","mo_gas_optics.jl"))
include(joinpath("rrtmgp","kernels","mo_gas_optics_kernels.jl"))
include(joinpath("rrtmgp","mo_gas_concentrations.jl"))
# include(joinpath("rrtmgp","mo_gas_optics_rrtmgp.jl"))
include(joinpath("examples","rfmip-clear-sky","mo_rfmip_io.jl"))
end # module
