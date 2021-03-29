using Test

@testset "RRTMGP gray radiation tests" begin

    include("gray_atm.jl")

end

@testset "RRTMGP clear-sky tests" begin

    include("rfmip_clear_sky_lw.jl")
    include("rfmip_clear_sky_sw.jl")

end
