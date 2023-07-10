include("all_sky_utils.jl")

context = ClimaComms.context()
@testset "Cloudy (all-sky, Two-stream calculations using lookup table method" begin
    @time all_sky(context, TwoStream, Float64; ncol = 128, use_lut = true, cldfrac = Float64(1))
end
@testset "Cloudy (all-sky), Two-stream calculations using Pade method" begin
    @time all_sky(context, TwoStream, Float64; ncol = 128, use_lut = false, cldfrac = Float64(1))
end
