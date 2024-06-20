FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64
include("gray_atm_utils.jl")

context = ClimaComms.context()
@time gray_atmos_lw_equil(context, NoScatLWRTE, FT)
@time gray_atmos_lw_equil(context, TwoStreamLWRTE, FT)

@time gray_atmos_sw_test(context, NoScatSWRTE, FT, 1)
@time gray_atmos_sw_test(context, TwoStreamSWRTE, FT, 1)
