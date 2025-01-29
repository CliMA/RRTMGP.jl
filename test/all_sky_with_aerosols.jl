FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64
include("all_sky_with_aerosols_utils.jl")

context = ClimaComms.context()

toler_lw_noscat = Dict(Float64 => Float64(1e-5), Float32 => Float32(0.05))
toler_lw_2stream = Dict(Float64 => Float64(5), Float32 => Float32(5))
toler_sw = Dict(Float64 => Float64(1e-5), Float32 => Float32(0.06))

@testset "all-sky (gas + clouds + aerosols) calculations using lookup table method, with non-scattering LW and TwoStream SW solvers" begin
    @time all_sky_with_aerosols(
        context,
        NoScatLWRTE,
        TwoStreamSWRTE,
        FT,
        toler_lw_noscat,
        toler_sw;
        ncol = 128,
        cldfrac = FT(1),
    )
end

@testset "all-sky (gas + clouds + aerosols) Two-stream calculations using lookup table method" begin
    @time all_sky_with_aerosols(
        context,
        TwoStreamLWRTE,
        TwoStreamSWRTE,
        FT,
        toler_lw_2stream,
        toler_sw;
        ncol = 128,
        cldfrac = FT(1),
    )
end
