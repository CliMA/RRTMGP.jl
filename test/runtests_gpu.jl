using Test

@testset "RRTMGP gray radiation tests" begin
    include("gray_atm.jl")
end
println("================================================================\n\n\n")

@testset "RRTMGP clear-sky tests" begin
    include("clear_sky.jl")
end
println("================================================================\n\n\n")

@testset "RRTMGP all-sky (cloudy) tests" begin
    include("all_sky.jl")
end
println("================================================================\n\n\n")
