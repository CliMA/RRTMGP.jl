module RTE

include(joinpath("rte","mo_util_array.jl"))
include(joinpath("rte","kernels","mo_optical_props_kernels.jl"))
include(joinpath("rte","mo_optical_props.jl"))

end # module
