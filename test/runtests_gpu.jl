using Test

@testset "RRTMGP gray radiation tests" begin

    include("gray_atm.jl")

end

@testset "RRTMGP clear-sky tests" begin
    include("clear_sky.jl")
end

@testset "RRTMGP clear-sky tests" begin
    include("all_sky.jl")
end
