FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64
include("all_sky_utils.jl")

context = ClimaComms.context()
@testset "Cloudy (all-sky, Two-stream calculations using lookup table method" begin
    @time all_sky(context, TwoStream, FT; ncol = 128, use_lut = true, cldfrac = FT(1))
end
@testset "Cloudy (all-sky), Two-stream calculations using Pade method" begin
    @time all_sky(context, TwoStream, FT; ncol = 128, use_lut = false, cldfrac = FT(1))
end
@testset "Cloudy (all-sky, Two-stream calculations using lookup table method with aerosols" begin
    @time all_sky(context, TwoStream, FT; ncol = 128, use_lut = true, cldfrac = FT(1), aero = true)
end
@testset "Cloudy (all-sky), Two-stream calculations using Pade method with aerosols" begin
    @time all_sky(context, TwoStream, FT; ncol = 128, use_lut = false, cldfrac = FT(1), aero = true)
end
