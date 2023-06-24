include("gray_atm_utils.jl")

context = ClimaComms.context()
@time gray_atmos_lw_equil(context, OneScalar, Float64)
@time gray_atmos_lw_equil(context, TwoStream, Float64)

@time gray_atmos_sw_test(context, OneScalar, Float64, 1)
@time gray_atmos_sw_test(context, TwoStream, Float64, 1)
