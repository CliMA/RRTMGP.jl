using Test
import RRTMGP
using Aqua

@testset "Aqua tests - unbound args" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    Aqua.test_unbound_args(RRTMGP)

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    Aqua.test_ambiguities(RRTMGP; recursive = true)
end
