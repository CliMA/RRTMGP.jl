using Test

color1 = 130

printstyled("\n\nRRTMGP gray radiation tests\n", color = color1)
printstyled("===============================\n\n", color = color1)

@testset "RRTMGP gray radiation tests" begin
    include("gray_atm.jl")
end

printstyled("\n\nRRTMGP clear-sky tests \n", color = color1)
printstyled("==========================\n\n", color = color1)
@testset "RRTMGP clear-sky tests" begin
    context = ClimaComms.context()

    toler_lw_noscat = Dict(Float64 => Float64(1e-4), Float32 => Float32(0.05))
    toler_lw_2stream = Dict(Float64 => Float64(4.5), Float32 => Float32(4.5))
    toler_sw = Dict(Float64 => Float64(1e-3), Float32 => Float32(0.04))

    include("clear_sky_utils.jl")

    for FT in (Float32, Float64)
        clear_sky(context, OneScalar, TwoStream, NoScatLWRTE, TwoStreamSWRTE, VmrGM, FT, toler_lw_noscat, toler_sw)
        clear_sky(context, TwoStream, TwoStream, TwoStreamLWRTE, TwoStreamSWRTE, VmrGM, FT, toler_lw_2stream, toler_sw)
    end
end

printstyled("\n\nRRTMGP all-sky (cloudy) tests\n", color = color1)
printstyled("=================================\n\n", color = color1)
@testset "RRTMGP all-sky (cloudy) tests" begin
    context = ClimaComms.context()
    include("all_sky_utils.jl")
    for FT in (Float32, Float64)
        all_sky(
            context,
            TwoStream,
            TwoStream,
            TwoStreamLWRTE,
            TwoStreamSWRTE,
            FT;
            ncol = 128,
            use_lut = true,
            cldfrac = FT(1),
        )
    end
end

printstyled("\n\nAqua tests\n", color = color1)
printstyled("==============\n\n", color = color1)
@testset "Aqua" begin
    include("aqua.jl")
end
printstyled("****************************************************************\n", color = color1)
