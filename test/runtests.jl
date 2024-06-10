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
    include("clear_sky_utils.jl")
    for FT in (Float32, Float64)
        clear_sky(context, TwoStream, TwoStream, TwoStreamLWRTE, TwoStreamSWRTE, VmrGM, FT)
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
