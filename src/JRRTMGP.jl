module JRRTMGP

include("misc.jl")
include("fortran_intrinsics.jl")
include(joinpath("rte","mo_util_array.jl"))
include(joinpath("rte","kernels","mo_fluxes_broadband_kernels.jl"))
include(joinpath("rte","kernels","mo_rte_solver_kernels.jl"))
include(joinpath("rte","mo_optical_props.jl"))
include(joinpath("rte","mo_fluxes.jl"))
include(joinpath("rte","mo_source_functions.jl"))
include(joinpath("extensions","mo_fluxes_byband_kernels.jl"))
# include(joinpath("extensions","mo_fluxes_byband.jl"))
include(joinpath("rte","mo_rte_sw.jl"))

include(joinpath("rrtmgp","mo_util_string.jl"))
include(joinpath("rrtmgp","mo_rrtmgp_constants.jl"))
include(joinpath("rrtmgp","kernels","mo_reorder_kernels.jl"))
include(joinpath("rrtmgp","mo_util_reorder.jl"))
include(joinpath("rrtmgp","mo_gas_optics.jl"))
include(joinpath("rrtmgp","kernels","mo_gas_optics_kernels.jl"))
include(joinpath("rrtmgp","mo_gas_concentrations.jl"))
include(joinpath("rrtmgp","mo_gas_optics_rrtmgp.jl"))
include(joinpath("examples","mo_load_coefficients.jl"))
include(joinpath("examples","rfmip-clear-sky","mo_rfmip_io.jl"))

# Drivers:
include(joinpath("examples","rfmip-clear-sky","rrtmgp_rfmip_sw.jl"))

end # module
