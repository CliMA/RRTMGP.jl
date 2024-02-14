FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64
FT = Float32

include("clear_sky_utils.jl")

context = ClimaComms.context()
@testset "testing clear sky 2-stream solver" begin
    @time clear_sky(context, TwoStream, SourceLW2Str, VmrGM, FT)
end
