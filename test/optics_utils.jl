using Test
using RRTMGP
using RRTMGP.Optics:
    loc_lower,
    interp1d_equispaced,
    interp1d_loc_factor,
    compute_lookup_cld_log_props,
    add_cloud_optics_1scalar!,
    add_cloud_optics_2stream!

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
    @test interp1d_equispaced(0.72, xeq, yeq) ==
          yeq[15] * (1 - 0.4) + yeq[16] * 0.4 # interior data point
    @test interp1d_equispaced(1.10, xeq, yeq) == yeq[23] # coincident with an internal data point

    # test interp1d_loc_factor for non-equispaced data
    @test interp1d_loc_factor(-0.3, x) == (1, FT(0)) # below lower limit
    @test interp1d_loc_factor(1.55, x) == (n - 1, FT(1)) # above upper limit
    @test interp1d_loc_factor(1.10, x) == (29, FT(0)) # coincident with an internal data point
    loc, fac = interp1d_loc_factor(FT(1.02), x) # internal point
    @test loc == 25 && fac ≈ FT(0.8)

    # test compute_lookup_cld_log_props on a log-uniform grid
    nsize_log = 64
    r_lwr, r_upr = FT(5.0), FT(1100.0)
    r_grid = exp.(range(log(r_lwr), log(r_upr), length = nsize_log))
    # Table linear in log(r) is reconstructed to machine precision everywhere
    ext_tbl = FT(0.2) .- FT(0.02) .* log.(r_grid)
    ssa_tbl = fill(FT(0.9), nsize_log)
    asy_tbl = fill(FT(0.85), nsize_log)
    path = FT(10.0)

    # Zero water path returns zeros
    @test compute_lookup_cld_log_props(
        nsize_log,
        r_lwr,
        r_upr,
        ext_tbl,
        ssa_tbl,
        asy_tbl,
        FT(100.0),
        FT(0.0),
    ) == (FT(0), FT(0), FT(0))

    # Clamping below r_lwr and above r_upr
    τ_low, τs_low, τsg_low = compute_lookup_cld_log_props(
        nsize_log,
        r_lwr,
        r_upr,
        ext_tbl,
        ssa_tbl,
        asy_tbl,
        FT(2.0),
        path,
    )
    @test τ_low ≈ ext_tbl[1] * path
    @test τs_low ≈ ssa_tbl[1] * τ_low
    @test τsg_low ≈ asy_tbl[1] * τs_low

    τ_high, τs_high, τsg_high = compute_lookup_cld_log_props(
        nsize_log,
        r_lwr,
        r_upr,
        ext_tbl,
        ssa_tbl,
        asy_tbl,
        FT(2000.0),
        path,
    )
    @test τ_high ≈ ext_tbl[end] * path
    @test τs_high ≈ ssa_tbl[end] * τ_high
    @test τsg_high ≈ asy_tbl[end] * τs_high

    # Interior point (r = 300 µm > 90 µm ice ceiling)
    r_mid = FT(300.0)
    τ_mid, τs_mid, τsg_mid = compute_lookup_cld_log_props(
        nsize_log,
        r_lwr,
        r_upr,
        ext_tbl,
        ssa_tbl,
        asy_tbl,
        r_mid,
        path,
    )
    @test τ_mid ≈ (FT(0.2) - FT(0.02) * log(r_mid)) * path
    @test τs_mid ≈ FT(0.9) * τ_mid
    @test τsg_mid ≈ FT(0.85) * τs_mid
end
