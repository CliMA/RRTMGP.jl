FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64
include("gray_atm_utils.jl")

context = ClimaComms.context()
@time gray_atmos_lw_equil(context, OneScalar, FT)
@time gray_atmos_lw_equil(context, TwoStream, FT)

@time gray_atmos_sw_test(context, OneScalar, FT, 1)
@time gray_atmos_sw_test(context, TwoStream, FT, 1)
