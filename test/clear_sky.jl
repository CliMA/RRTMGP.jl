FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64

include("clear_sky_utils.jl")

context = ClimaComms.context()

@testset "clear-sky: NoScatLWRTE (OneScalar optics) & TwoStreamSWRTE (TwoStream optics)" begin
    @time clear_sky(context, OneScalar, TwoStream, NoScatLWRTE, TwoStreamSWRTE, VmrGM, FT)
end

@testset "clear-sky: TwoStreamLWRTE (TwoStream optics) & TwoStreamSWRTE (TwoStream optics)" begin
    @time clear_sky(context, TwoStream, TwoStream, TwoStreamLWRTE, TwoStreamSWRTE, VmrGM, FT)
end
