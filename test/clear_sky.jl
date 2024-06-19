FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64

include("clear_sky_utils.jl")

context = ClimaComms.context()

toler_lw_noscat = Dict(Float64 => Float64(1e-4), Float32 => Float32(0.05))
toler_lw_2stream = Dict(Float64 => Float64(4.5), Float32 => Float32(4.5))
toler_sw = Dict(Float64 => Float64(1e-3), Float32 => Float32(0.04))

@testset "clear-sky: NoScatLWRTE (OneScalar optics) & TwoStreamSWRTE (TwoStream optics)" begin
    @time clear_sky(context, NoScatLWRTE, TwoStreamSWRTE, VmrGM, FT, toler_lw_noscat, toler_sw)
end

@testset "clear-sky: TwoStreamLWRTE (TwoStream optics) & TwoStreamSWRTE (TwoStream optics)" begin
    @time clear_sky(context, TwoStreamLWRTE, TwoStreamSWRTE, VmrGM, FT, toler_lw_2stream, toler_sw)
end
