using Test
using RRTMGP
using RRTMGP.Optics: loc_lower, interp1d_equispaced, interp1d_loc_factor

@testset "RRTMGP optics utilities tests" begin
    FT = Float64
    # test loc_lower for equispaced data
    Δx = 0.05
    xeq = Vector(FT(0):0.05:1.5)
    neq = length(xeq)
    @test loc_lower(-0.3, Δx, neq, xeq) == 1 # below lower limit
    @test loc_lower(1.55, Δx, neq, xeq) == neq - 1 # above upper limit
    @test loc_lower(0.72, Δx, neq, xeq) == 15 # interior point
    @test loc_lower(1.1, Δx, neq, xeq) == 23 # coincident with an internal data point

    # test loc_lower for non-equispaced data
    x = vcat(Vector(FT(0):0.05:0.8), Vector(FT(0.825):0.025:1.5))
    n = length(x)

    @test loc_lower(-0.3, x) == 1 # below lower limit
    @test loc_lower(1.55, x) == n - 1 # above upper limit
    @test loc_lower(0.72, x) == 15 # interior point
    @test loc_lower(1.02, x) == 25 # interior point
    @test loc_lower(1.10, x) == 29 # coincident with an internal data point

    # test interp1d_equispaced
    yeq = 3 .* xeq .+ 4

    @test interp1d_equispaced(-0.3, xeq, yeq) == yeq[1] # below lower limit
    @test interp1d_equispaced(1.55, xeq, yeq) == yeq[end] # above upper limit
    @test interp1d_equispaced(0.72, xeq, yeq) == yeq[15] * (1 - 0.4) + yeq[16] * 0.4 # interior data point
    @test interp1d_equispaced(1.10, xeq, yeq) == yeq[23] # coincident with an internal data point

    # test interp1d_loc_factor for non-equispaced data
    @test interp1d_loc_factor(-0.3, x) == (1, FT(0)) # below lower limit
    @test interp1d_loc_factor(1.55, x) == (n - 1, FT(1)) # above upper limit
    @test interp1d_loc_factor(1.10, x) == (29, FT(0)) # coincident with an internal data point
    loc, fac = interp1d_loc_factor(FT(1.02), x) # internal point
    @test loc == 25 && fac ≈ FT(0.8)
end
