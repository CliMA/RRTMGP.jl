using Test
import RRTMGP
import RRTMGP.AtmosphericStates: AtmosphericState
import RRTMGP.VolumeMixingRatios: VmrGM

# Closed-form unit tests for every level<->layer interpolation and bottom-
# extrapolation scheme in src/api/interpolation.jl. The laws being verified are
# the ones stated in that file's comment block: a constant lapse rate between
# the two known points, hydrostatic balance for isothermal columns, and the
# isentropic p ∝ T^(cₚ/R) law otherwise.

# Minimal clear-sky AtmosphericState (VmrGM storage), mirroring the helper in
# grid_adaptation.jl but under a different name so both files can be included
# into one test session.
function _make_interp_as(FT, nlay, ncol)
    layerdata = zeros(FT, 4, nlay, ncol)
    p_lev = zeros(FT, nlay + 1, ncol)
    t_lev = zeros(FT, nlay + 1, ncol)
    t_sfc = fill(FT(300), ncol)
    vmr = VmrGM(zeros(FT, nlay, ncol), zeros(FT, nlay, ncol), zeros(FT, 3))
    return AtmosphericState(
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
end

@testset "interp! closed forms (all schemes)" begin
    for FT in (Float32, Float64)
        pꜜ = FT[100000, 90000]
        Tꜜ = FT[300, 290]
        pꜛ = FT[80000, 70000]
        Tꜛ = FT[280, 270]
        p = zeros(FT, 2)
        T = zeros(FT, 2)

        # UniformZ: T is the arithmetic mean; p follows the isentropic
        # p ∝ T^B law, checked independently by the bounds here and the
        # isothermal limit below.
        RRTMGP.interp!(RRTMGP.UniformZ(), p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
        @test T ≈ (Tꜜ .+ Tꜛ) ./ 2
        # the interpolated face lies between the adjacent layers
        @test all(min.(pꜜ, pꜛ) .< p .< max.(pꜜ, pꜛ))
        @test all(min.(Tꜜ, Tꜛ) .< T .< max.(Tꜜ, Tꜛ))

        # isothermal limit: hydrostatic balance reduces UniformZ to the
        # geometric mean of the pressures.
        RRTMGP.interp!(RRTMGP.UniformZ(), p, T, pꜜ, Tꜜ, pꜛ, Tꜜ)
        @test T ≈ Tꜜ
        @test p ≈ sqrt.(pꜜ .* pꜛ)

        # UniformP: p is the arithmetic mean; T inverts the same power law,
        # checked independently by the bounds here and the isothermal limit
        # below.
        RRTMGP.interp!(RRTMGP.UniformP(), p, T, pꜜ, Tꜜ, pꜛ, Tꜛ)
        @test p ≈ (pꜜ .+ pꜛ) ./ 2
        @test all(min.(Tꜜ, Tꜛ) .< T .< max.(Tꜜ, Tꜛ))
        # isothermal limit (independent oracle): an isothermal layer stays
        # isothermal, whatever the pressure interpolation.
        RRTMGP.interp!(RRTMGP.UniformP(), p, T, pꜜ, Tꜜ, pꜛ, Tꜜ)
        @test T ≈ Tꜜ
        @test p ≈ (pꜜ .+ pꜛ) ./ 2

        # BestFit with distinct temperatures: T is linear in z and p follows
        # the power law fitted through the two layer points.
        z = FT[500]
        zꜜ = FT[0]
        zꜛ = FT[1000]
        pb = zeros(FT, 1)
        Tb = zeros(FT, 1)
        RRTMGP.interp!(
            RRTMGP.BestFit(),
            pb,
            Tb,
            z,
            pꜜ[1:1],
            Tꜜ[1:1],
            zꜜ,
            pꜛ[1:1],
            Tꜛ[1:1],
            zꜛ,
        )
        @test Tb ≈ [(Tꜜ[1] + Tꜛ[1]) / 2] # midpoint of a linear profile
        # the fitted pressure lies between the two layer pressures; its exact
        # value is checked by the dry-adiabat test below.
        @test pꜛ[1] < pb[1] < pꜜ[1]
    end
end

@testset "extrap! closed forms (all schemes)" begin
    for FT in (Float32, Float64)
        params = RRTMGP.default_parameters(FT)
        g = RRTMGP.Parameters.grav(params)
        cₚ = RRTMGP.Parameters.cp_d(params)
        R = RRTMGP.Parameters.R_d(params)

        p⁺ = FT[90000]
        T⁺ = FT[285]
        p⁺⁺ = FT[80000]
        T⁺⁺ = FT[275]
        Tₛ = FT[300]
        p = zeros(FT, 1)
        T = zeros(FT, 1)

        # GeometricMean continues the geometric progression of the two layers
        # above to the boundary face half a layer below, so the boundary/nearest
        # ratio is the square root of the layer-to-layer ratio.
        RRTMGP.extrap!(RRTMGP.GeometricMean(), p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
        @test T ./ T⁺ ≈ sqrt.(T⁺ ./ T⁺⁺)
        @test p ./ p⁺ ≈ sqrt.(p⁺ ./ p⁺⁺)
        @test all(T .> T⁺) && all(p .> p⁺) # continues downward

        # UniformZ: linear T continuation; the isothermal limit independently
        # pins the isentropic pressure law to the hydrostatic geometric mean.
        RRTMGP.extrap!(RRTMGP.UniformZ(), p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
        @test T ≈ (3 .* T⁺ .- T⁺⁺) ./ 2
        @test all(p .> p⁺) # continues downward below the lowest layer
        RRTMGP.extrap!(RRTMGP.UniformZ(), p, T, p⁺, T⁺, p⁺⁺, T⁺, Tₛ, params)
        @test T ≈ T⁺
        @test p ≈ sqrt.(p⁺ .* p⁺⁺)

        # UniformP: linear p continuation; an isothermal layer extrapolates
        # isothermally.
        RRTMGP.extrap!(RRTMGP.UniformP(), p, T, p⁺, T⁺, p⁺⁺, T⁺⁺, Tₛ, params)
        @test p ≈ (3 .* p⁺ .- p⁺⁺) ./ 2
        @test all(T .> T⁺) # continues downward below the lowest layer
        RRTMGP.extrap!(RRTMGP.UniformP(), p, T, p⁺, T⁺, p⁺⁺, T⁺, Tₛ, params)
        @test T ≈ T⁺
        @test p ≈ (3 .* p⁺ .- p⁺⁺) ./ 2

        # UseSurfaceTempAtBottom: T is the surface temperature; p follows the
        # dry-isentropic law p = p⁺ (T/T⁺)^(cₚ/R).
        RRTMGP.extrap!(
            RRTMGP.UseSurfaceTempAtBottom(),
            p,
            T,
            p⁺,
            T⁺,
            p⁺⁺,
            T⁺⁺,
            Tₛ,
            params,
        )
        @test T == Tₛ
        @test all(p .> p⁺) # Tₛ > T⁺, so the boundary pressure exceeds p⁺
        # identity limit: Tₛ = T⁺ recovers p = p⁺
        RRTMGP.extrap!(
            RRTMGP.UseSurfaceTempAtBottom(),
            p,
            T,
            p⁺,
            T⁺,
            p⁺⁺,
            T⁺⁺,
            T⁺,
            params,
        )
        @test p ≈ p⁺

        # HydrostaticBottom: dry-adiabatic lapse rate T = T⁺ + (g/cₚ)(z⁺ - z)
        # and the same isentropic pressure law.
        z = FT[0]
        z⁺ = FT[500]
        z⁺⁺ = FT[1500]
        RRTMGP.extrap!(
            RRTMGP.HydrostaticBottom(),
            p,
            T,
            z,
            p⁺,
            T⁺,
            z⁺,
            p⁺⁺,
            T⁺⁺,
            z⁺⁺,
            Tₛ,
            params,
        )
        @test T ≈ T⁺ .+ (g / cₚ) .* (z⁺ .- z)
        # pressure follows the isentropic law; its exact value is pinned by the
        # identity limit below and the dry-adiabat test.
        @test all(p .> p⁺) # continues downward below the lowest layer
        # identity limit: z = z⁺ recovers (p, T) = (p⁺, T⁺)
        RRTMGP.extrap!(
            RRTMGP.HydrostaticBottom(),
            p,
            T,
            z⁺,
            p⁺,
            T⁺,
            z⁺,
            p⁺⁺,
            T⁺⁺,
            z⁺⁺,
            Tₛ,
            params,
        )
        @test T ≈ T⁺
        @test p ≈ p⁺

        # BestFit extrapolation: linear T in z and the fitted power law.
        RRTMGP.extrap!(
            RRTMGP.BestFit(),
            p,
            T,
            z,
            p⁺,
            T⁺,
            z⁺,
            p⁺⁺,
            T⁺⁺,
            z⁺⁺,
            Tₛ,
            params,
        )
        @test T ≈ T⁺ .+ (T⁺⁺ .- T⁺) .* (z .- z⁺) ./ (z⁺⁺ .- z⁺)
        # the fitted pressure continues downward; its exact value is checked by
        # the dry-adiabat test below.
        @test all(p .> p⁺)
    end
end

# On a dry adiabat — T(z) linear with slope -g/cₚ and p(z) = p₀ (T/T₀)^(cₚ/R) —
# BestFit interpolation and HydrostaticBottom extrapolation are exact, so
# interpolate_levels! must reproduce the analytic level values from the layer
# values alone. This exercises the z-based branch of interpolate_levels!
# (center_z/face_z plumbing) end to end.
@testset "interpolate_levels! reproduces a dry adiabat (BestFit + HydrostaticBottom)" begin
    for FT in (Float32, Float64)
        params = RRTMGP.default_parameters(FT)
        g = RRTMGP.Parameters.grav(params)
        cₚ = RRTMGP.Parameters.cp_d(params)
        R = RRTMGP.Parameters.R_d(params)
        T₀ = FT(300)
        p₀ = FT(100000)
        T_ad(z) = T₀ - (g / cₚ) * z
        p_ad(z) = p₀ * (T_ad(z) / T₀)^(cₚ / R)

        nlay, ncol = 8, 2
        z_lev = collect(FT, range(FT(0), FT(8000); length = nlay + 1))
        z_lay = (z_lev[1:(end - 1)] .+ z_lev[2:end]) ./ 2
        face_z = repeat(z_lev, 1, ncol)
        center_z = repeat(z_lay, 1, ncol)

        as = _make_interp_as(FT, nlay, ncol)
        RRTMGP.AtmosphericStates.getview_p_lay(as) .= p_ad.(center_z)
        RRTMGP.AtmosphericStates.getview_t_lay(as) .= T_ad.(center_z)

        RRTMGP.interpolate_levels!(
            as,
            RRTMGP.BestFit(),
            RRTMGP.HydrostaticBottom(),
            params;
            center_z,
            face_z,
        )

        rtol = sqrt(eps(FT))
        @test as.t_lev ≈ T_ad.(face_z) rtol = rtol
        @test as.p_lev ≈ p_ad.(face_z) rtol = rtol
    end
end

# UniformP divides by log(pꜛ/pꜜ) and therefore assumes distinct layer
# pressures (see the comment in interpolation.jl). Pin the precondition: equal
# pressures produce NaN temperatures rather than silently wrong values.
@testset "UniformP precondition: distinct layer pressures" begin
    p = zeros(1)
    T = zeros(1)
    RRTMGP.interp!(RRTMGP.UniformP(), p, T, [500.0], [250.0], [500.0], [260.0])
    @test p == [500.0]
    @test isnan(T[1])
end
