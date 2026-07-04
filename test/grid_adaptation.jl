using Test
import RRTMGP
import RRTMGP.AtmosphericStates: AtmosphericState
import RRTMGP.Vmrs: VmrGM

# Build a minimal clear-sky AtmosphericState (no clouds/aerosols, VmrGM storage).
function _make_as(FT, nlay, ncol)
    layerdata = zeros(FT, 4, nlay, ncol)
    p_lev = zeros(FT, nlay + 1, ncol)
    t_lev = zeros(FT, nlay + 1, ncol)
    t_sfc = fill(FT(300), ncol)
    vmr = VmrGM(zeros(FT, nlay, ncol), zeros(FT, nlay, ncol), zeros(FT, 3))
    as = AtmosphericState(
        nothing,
        nothing,
        layerdata,
        p_lev,
        t_lev,
        t_sfc,
        vmr,
        nothing,
        nothing,
    )
    return as
end

@testset "interp!/extrap! arithmetic" begin
    p = zeros(2)
    T = zeros(2)
    pꜜ = [1.0, 2.0]
    Tꜜ = [10.0, 20.0]
    pꜛ = [3.0, 4.0]
    Tꜛ = [30.0, 40.0]

    RRTMGP.interp!(RRTMGP.ArithmeticMean(), p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
    @test p ≈ (pꜜ .+ pꜛ) ./ 2
    @test T ≈ (Tꜜ .+ Tꜛ) ./ 2

    RRTMGP.interp!(RRTMGP.GeometricMean(), p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
    @test p ≈ sqrt.(pꜜ .* pꜛ)
    @test T ≈ sqrt.(Tꜜ .* Tꜛ)

    # BestFit with equal temperatures -> pressure power law in z.
    pb = zeros(1)
    Tb = zeros(1)
    RRTMGP.interp!(
        RRTMGP.BestFit(),
        pb,
        Tb,
        [1.5],
        [1.0],
        [280.0],
        [1.0],
        [3.0],
        [280.0],
        [2.0],
    )
    @test Tb ≈ [280.0]
    @test pb[1] ≈ 1.0 * (3.0 / 1.0)^((1.5 - 1.0) / (2.0 - 1.0))

    # ArithmeticMean extrapolation (params unused for this method).
    pe = zeros(1)
    Te = zeros(1)
    RRTMGP.extrap!(
        RRTMGP.ArithmeticMean(),
        pe,
        Te,
        [400.0],
        [250.0],
        [600.0],
        [270.0],
        [300.0],
        nothing,
    )
    @test pe ≈ [(3 * 400.0 - 600.0) / 2]
    @test Te ≈ [(3 * 250.0 - 270.0) / 2]

    @test RRTMGP.requires_z(RRTMGP.BestFit())
    @test RRTMGP.requires_z(RRTMGP.HydrostaticBottom())
    @test !RRTMGP.requires_z(RRTMGP.ArithmeticMean())
    @test !RRTMGP.requires_z(RRTMGP.UniformZ())
    @test !RRTMGP.requires_z(RRTMGP.NoInterpolation())
end

@testset "interpolate_levels! (ArithmeticMean)" begin
    FT = Float64
    nlay, ncol = 4, 2
    as = _make_as(FT, nlay, ncol)
    p_lay = RRTMGP.AtmosphericStates.getview_p_lay(as)
    t_lay = RRTMGP.AtmosphericStates.getview_t_lay(as)
    p_lay .= FT[1000 1100; 800 850; 600 620; 400 410]
    t_lay .= FT[290 292; 270 272; 250 252; 230 232]

    RRTMGP.interpolate_levels!(
        as,
        RRTMGP.ArithmeticMean(),
        RRTMGP.SameAsInterpolation(),
        nothing,
    )
    p_lev = as.p_lev
    t_lev = as.t_lev
    # interior faces are arithmetic means of the adjacent layers
    for k in 2:nlay
        @test p_lev[k, :] ≈ (p_lay[k - 1, :] .+ p_lay[k, :]) ./ 2
        @test t_lev[k, :] ≈ (t_lay[k - 1, :] .+ t_lay[k, :]) ./ 2
    end
    # top face: extrapolation from the two layers below it
    @test t_lev[nlay + 1, :] ≈ (3 .* t_lay[nlay, :] .- t_lay[nlay - 1, :]) ./ 2
    # bottom face: SameAsInterpolation -> ArithmeticMean extrapolation
    @test t_lev[1, :] ≈ (3 .* t_lay[1, :] .- t_lay[2, :]) ./ 2

    # NoInterpolation is a no-op
    @test RRTMGP.interpolate_levels!(
        as,
        RRTMGP.NoInterpolation(),
        RRTMGP.SameAsInterpolation(),
        nothing,
    ) === as
end

@testset "add_isothermal_boundary_layer! + clip!" begin
    FT = Float64
    nlay, ncol = 4, 1   # nlay includes the extra (top) isothermal layer
    as = _make_as(FT, nlay, ncol)
    p_lay = RRTMGP.AtmosphericStates.getview_p_lay(as)
    t_lay = RRTMGP.AtmosphericStates.getview_t_lay(as)
    rh = RRTMGP.AtmosphericStates.getview_rel_hum(as)
    as.p_lev[:, 1] .= FT[1000, 800, 600, 400, 200]
    as.t_lev[:, 1] .= FT[300, 280, 260, 240, 220]
    rh[1:(nlay - 1), 1] .= FT[0.8, 0.6, 0.4]
    as.vmr.vmr_h2o[:, 1] .= FT[0.05, 0.04, 0.03, 0.0]

    p_min = FT(10)
    RRTMGP.add_isothermal_boundary_layer!(as, p_min)
    @test as.p_lev[end, 1] == p_min
    @test p_lay[end, 1] == (as.p_lev[end - 1, 1] + p_min) / 2
    @test t_lay[end, 1] == as.t_lev[end - 1, 1]
    @test as.t_lev[end, 1] == as.t_lev[end - 1, 1]
    @test rh[end, 1] == rh[end - 1, 1] == FT(0.4)
    @test as.vmr.vmr_h2o[end, 1] == as.vmr.vmr_h2o[end - 1, 1] == FT(0.03)

    # clip!
    p_lay[1, 1] = FT(-5)
    as.vmr.vmr_h2o[1, 1] = FT(-1)
    RRTMGP.clip!(as, p_min)
    @test p_lay[1, 1] == p_min
    @test as.vmr.vmr_h2o[1, 1] == 0
end
