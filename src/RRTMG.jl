module RRTMG

# data_dir = joinpath("rrtmgp","data")

include("fortran_intrinsics.jl")
include(joinpath("rrtmgp","mo_util_string.jl"))
# include(joinpath("rrtmgp","mo_util_reorder.jl"))
# include(joinpath("rrtmgp","mo_gas_optics.jl"))
# include(joinpath("rrtmgp","mo_rrtmgp_constants.jl"))
# include(joinpath("rrtmgp","kernels","mo_reorder_kernels.jl"))
include(joinpath("rrtmgp","kernels","mo_gas_optics_kernels.jl"))
include(joinpath("rrtmgp","mo_gas_concentrations.jl"))
# include(joinpath("rrtmgp","mo_gas_optics_rrtmgp.jl"))

end # module
