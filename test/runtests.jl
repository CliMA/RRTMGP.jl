using Test

color1 = 130

@testset "Aerosol name map tests" begin
    include("aerosol_name_map.jl")
end

@testset "Grid adaptation tests" begin
    include("grid_adaptation.jl")
end

@testset "RRTMGP gray radiation tests" begin
    include("gray_atm.jl")
end

@testset "RRTMGP standalone tests" begin
    include("standalone.jl")
end

@testset "Flux-getter scalar-indexing guard (GPU)" begin
    include("scalar_indexing.jl")
end

printstyled("\n\nRRTMGP clear-sky tests \n", color = color1)
printstyled("==========================\n\n", color = color1)
@testset "RRTMGP clear-sky tests" begin
    context = ClimaComms.context()

    toler_lw_noscat =
        Dict(Float64 => Float64(1e-4), Float32 => Float32(5e-3))
    toler_lw_2stream = Dict(Float64 => Float64(4.5), Float32 => Float32(4.5))
    toler_sw = Dict(Float64 => Float64(1e-3), Float32 => Float32(0.04))

    include("clear_sky_utils.jl")

    for FT in (Float32, Float64)
        clear_sky(
            context,
            NoScatLWRTE,
            TwoStreamSWRTE,
            VmrGM,
            FT,
            toler_lw_noscat,
            toler_sw,
        )
        clear_sky(
            context,
            TwoStreamLWRTE,
            TwoStreamSWRTE,
            VmrGM,
            FT,
            toler_lw_2stream,
            toler_sw,
        )
    end
end

printstyled("\n\nRRTMGP cloudy-sky (gas + clouds) tests\n", color = color1)
printstyled("=================================\n\n", color = color1)
@testset "RRTMGP cloudy-sky (gas + clouds) tests" begin
    context = ClimaComms.context()

    toler_lw_noscat =
        Dict(Float64 => Float64(1e-5), Float32 => Float32(5e-3))
    toler_lw_2stream = Dict(Float64 => Float64(5), Float32 => Float32(5))
    toler_sw = Dict(Float64 => Float64(1e-5), Float32 => Float32(0.06))

    include("cloudy_sky_utils.jl")
    for FT in (Float32, Float64)
        cloudy_sky(
            context,
            NoScatLWRTE,
            TwoStreamSWRTE,
            FT,
            toler_lw_noscat,
            toler_sw;
            ncol = 128,
            cldfrac = FT(1),
        )
        cloudy_sky(
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
end

printstyled(
    "\n\nRRTMGP partial cloud fraction reproducibility tests\n",
    color = color1,
)
printstyled("=================================\n\n", color = color1)
@testset "RRTMGP partial cloud fraction reproducibility tests" begin
    include("partial_cloud_fraction.jl")
    context = ClimaComms.context()
    partial_cloud_fraction_test(context, Float64)
end

printstyled(
    "\n\nRRTMGP all-sky (gas + clouds + aerosols) tests\n",
    color = color1,
)
printstyled("=================================\n\n", color = color1)
@testset "RRTMGP all-sky (gas + clouds + aerosols) tests" begin
    context = ClimaComms.context()

    toler_lw_noscat =
        Dict(Float64 => Float64(1e-5), Float32 => Float32(5e-3))
    toler_lw_2stream = Dict(Float64 => Float64(5), Float32 => Float32(5))
    toler_sw = Dict(Float64 => Float64(1e-5), Float32 => Float32(0.06))

    include("all_sky_with_aerosols_utils.jl")
    for FT in (Float32, Float64)
        all_sky_with_aerosols(
            context,
            NoScatLWRTE,
            TwoStreamSWRTE,
            FT,
            toler_lw_noscat,
            toler_sw;
            ncol = 128,
            cldfrac = FT(1),
        )
        all_sky_with_aerosols(
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
end

printstyled("\n\ncos_zenith edge cases tests\n", color = color1)
printstyled("=================================\n\n", color = color1)
@testset "cos_zenith edge cases tests" begin
    include("cos_zenith_edge_cases.jl")
end

printstyled("\n\nFloat32 consistency tests\n", color = color1)
printstyled("=================================\n\n", color = color1)
include("float32_consistency.jl")

printstyled("\n\nLookup-table serialization tests\n", color = color1)
printstyled("=================================\n\n", color = color1)
include("lookup_serialization.jl")

printstyled("\n\nStandalone spectral solve(profile) tests\n", color = color1)
printstyled("=================================\n\n", color = color1)
include("standalone_spectral.jl")

printstyled("\n\nOptics utilities tests\n", color = color1)
printstyled("==============\n\n", color = color1)
include("optics_utils.jl")
printstyled(
    "****************************************************************\n",
    color = color1,
)

printstyled("\n\nAqua tests\n", color = color1)
printstyled("==============\n\n", color = color1)
@testset "Aqua" begin
    include("aqua.jl")
end
printstyled(
    "****************************************************************\n",
    color = color1,
)
