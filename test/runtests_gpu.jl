using Test

@testset "RRTMGP gray radiation tests" begin

    include("gray_atm.jl")
    println("*****************************")

end

@testset "RRTMGP clear-sky tests" begin

    include("rfmip_clear_sky_lw.jl")
    println("*****************************")
    include("rfmip_clear_sky_sw.jl")
    println("*****************************")

end
