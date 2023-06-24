include("clear_sky_utils.jl")

context = ClimaComms.context()
@testset "testing clear sky 2-stream solver" begin
    @time clear_sky(context, TwoStream, SourceLW2Str, VmrGM, Float64)
end
