using Test
using RRTMGP
using Aqua

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    # Aqua.test_unbound_args(RRTMGP)
    ua = Aqua.detect_unbound_args_recursively(RRTMGP)
    @test length(ua) == 0

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(RRTMGP; recursive = true)
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("RRTMGP", pkgdir(last(x).module)), ambs)

    # Uncomment for debugging:
    # for method_ambiguity in ambs
    #     @show method_ambiguity
    # end
    @test length(ambs) == 0
end

@testset "Aqua tests (additional)" begin
    Aqua.test_undefined_exports(RRTMGP)
    Aqua.test_stale_deps(RRTMGP)
    Aqua.test_deps_compat(RRTMGP)
    Aqua.test_project_extras(RRTMGP)
    Aqua.test_piracies(RRTMGP)
end

nothing
