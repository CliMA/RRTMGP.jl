using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
import RRTMGP.AngularDiscretizations: AngularDiscretization
import RRTMGP.BCs: LwBCs
import RRTMGP.RTE: NoScatLWRTE
import RRTMGP.RTESolver: solve_lw!

# The spectral testset at the bottom reuses the clear-sky state reader, which is
# written against the same names its own driver brings into scope.
using NCDatasets
using RRTMGP.LookUpTables
using RRTMGP.VolumeMixingRatios
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.Fluxes
using RRTMGP.ArtifactPaths
include("reference_files.jl")
include("read_clear_sky.jl")

# The non-scattering longwave integrates the Schwarzschild equation along a few
# discrete zenith angles and sums with quadrature weights. These tests check the
# quadrature itself against the exact hemispheric integral, and then check that
# the solver reproduces that behavior end to end.

# Exact hemispheric transmittance of an absorbing slab of optical depth τ,
# 2E₃(τ) = 2∫₀¹ exp(-τ/μ) μ dμ, by fine Simpson quadrature in μ (no dependency
# on a special-function package).
function two_E3(τ; n = 100_001)
    f(μ) = μ == 0 ? zero(τ) : exp(-τ / μ) * μ
    h = 1 / (n - 1)
    s = f(0.0) + f(1.0)
    for i in 1:(n - 2)
        s += (isodd(i) ? 4 : 2) * f(i * h)
    end
    return 2 * s * h / 3
end

const ΤS = (0.05, 0.2, 0.5, 1.0, 2.0, 5.0)

@testset "quadrature weights and secants" begin
    context = ClimaComms.context()
    gp = RRTMGP.RRTMGPGridParams(
        Float64;
        context,
        domain_nlay = 1,
        ncol = 1,
    )
    for n in 1:4
        ad = AngularDiscretization(gp, n)
        w = Array(ad.gauss_wts)
        D = Array(ad.gauss_Ds)
        @test ad.n_gauss_angles == n
        @test length(w) == n && length(D) == n
        # Weights normalized so that Σ wᵢ π Iᵢ recovers the hemispheric flux
        @test sum(w) ≈ 1
        @test all(>(0), w)
        # Secants exceed 1 (D = 1/μ with μ ∈ (0, 1]) and run steepest first
        @test all(>(1), D)
        @test issorted(D; rev = true)
    end
    # only 1-4 angles are tabulated
    @test_throws AssertionError AngularDiscretization(gp, 0)
    @test_throws AssertionError AngularDiscretization(gp, 5)
end

# The angular integral is exactly what more angles are meant to improve: the
# worst-case error over a range of optical depths must fall with each added
# angle. The single-angle case is the diffusivity approximation.
@testset "more angles integrate the hemisphere better" begin
    context = ClimaComms.context()
    gp = RRTMGP.RRTMGPGridParams(
        Float64;
        context,
        domain_nlay = 1,
        ncol = 1,
    )
    worst = map(1:4) do n
        ad = AngularDiscretization(gp, n)
        w, D = Array(ad.gauss_wts), Array(ad.gauss_Ds)
        maximum(τ -> abs(sum(w .* exp.(-τ .* D)) - two_E3(τ)), ΤS)
    end
    @test issorted(worst; rev = true)      # strictly improving
    @test worst[1] > 1e-2                  # one angle: a percent-level error
    @test worst[4] < 1e-3                  # four angles: better by >10x
end

# Exact check of the transport itself. For an isothermal layer the linear-in-τ
# correction vanishes (`lay_source == lev_source`), so one angle contributes
# exactly π·w·B·(1 - exp(-τD)) of downward flux at the surface, and the angles
# sum to π·B·(1 - Σ wᵢ exp(-τ Dᵢ)), which converges to the exact
# π·B·(1 - 2E₃(τ)). This pins the weights, the secant scaling, and the
# intensity/flux normalization without going through a full driver.
#
# `rte_lw_noscat_one_angle!` is a per-column kernel body: it writes its buffer
# element by element, which is legal on the GPU only from inside a launch. It
# runs here through `ClimaComms.@threaded` — a plain loop on the CPU, a kernel
# on the GPU — exactly as the real driver calls it, so the test exercises the
# device path it will run on.
@testset "one angle of transport is exact for an isothermal layer" begin
    FT = Float64
    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    nlay, ncol, nlev = 1, 1, 2
    params = RRTMGP.default_parameters(FT)
    gp = RRTMGP.RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
    B = FT(0.5)   # gray source function [W/m²/sr], any positive value will do
    τ_layer = FT(0.7)

    function surface_flux_dn(Ds, w_μ)
        op = Optics.OneScalar(gp)
        src = Sources.SourceLWNoScat(gp; params)
        flux = Fluxes.FluxLW(gp)
        bcs = LwBCs(fill!(DA{FT}(undef, 1, ncol), FT(1)), nothing)
        fill!(op.τ, τ_layer)
        fill!(src.lay_source, B)   # isothermal: layer and level sources agree,
        fill!(src.lev_source, B)   # so the linear-in-τ term drops out
        fill!(src.sfc_source, B)
        ClimaComms.@threaded device for gcol in 1:ncol
            RRTMGP.RTESolver.rte_lw_noscat_one_angle!(
                src,
                bcs,
                op,
                Ds,
                w_μ,
                gcol,
                flux,
                1,
                1,
                nlay,
                nlev,
            )
        end
        return Array(flux.flux_dn)[1, 1]
    end

    total = zero(FT)
    for n in 1:4
        ad = AngularDiscretization(gp, n)
        D, w = Array(ad.gauss_Ds), Array(ad.gauss_wts)
        n == 4 && (total = sum(i -> surface_flux_dn(D[i], w[i]), 1:n))
        for i in 1:n
            expected = FT(π) * w[i] * B * (1 - exp(-τ_layer * D[i]))
            @test surface_flux_dn(D[i], w[i]) ≈ expected rtol = 1e-14
        end
    end
    # the four-angle sum lands on the exact hemispheric flux
    exact = FT(π) * B * (1 - two_E3(τ_layer))
    @test isapprox(total, exact; rtol = 2e-3)
end

# The default must stay the single-angle diffusivity approximation, bit for bit:
# the reference-comparison tolerances and the Float32 ratchets are all measured
# against it.
@testset "one angle is the default and is unchanged" begin
    FT = Float64
    context = ClimaComms.context()
    DA = ClimaComms.array_type(ClimaComms.device(context))
    params = RRTMGP.default_parameters(FT)
    nlay, ncol = 30, 4

    function gray_noscat_flux(; kwargs...)
        gp = RRTMGP.RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
        as = RRTMGP.AtmosphericStates.setup_gray_as_pr_grid(
            context,
            nlay,
            DA{FT}([-60, -20, 20, 60]),
            FT(1.0e5),
            FT(9.0e3),
            RRTMGP.AtmosphericStates.GrayOpticalThicknessOGorman2008(FT),
            params,
            DA,
        )
        sfc_emis = fill!(DA{FT}(undef, 1, ncol), FT(1))
        slv = NoScatLWRTE(gp; params, sfc_emis, inc_flux = nothing, kwargs...)
        solve_lw!(slv, as)
        return copy(Array(slv.flux.flux_up)), slv.angle_disc.n_gauss_angles
    end

    F_default, n_default = gray_noscat_flux()
    F_explicit, n_explicit = gray_noscat_flux(; n_gauss_angles = 1)
    @test n_default == 1
    @test n_explicit == 1
    @test F_default == F_explicit # bit for bit

    # the same default holds through the Layer-2 solver
    solver = RRTMGP.solve_gray(FT; nlay, ncol).solver
    @test !hasproperty(solver.lws, :angle_disc) ||
          solver.lws.angle_disc.n_gauss_angles == 1
end

# The spectral path is the only one that accumulates over g-points as well as
# angles, so the flattened (igpt, imu) counter is exercised here alone. This
# also guards against the angle loop collapsing back to a single angle. The
# stored rte-rrtmgp reference fluxes are themselves a one-angle calculation, so
# adding angles moves the result away from them; what must hold is that
# successive refinements converge.
@testset "multi-angle clear-sky (spectral) solve" begin
    FT = Float64
    context = ClimaComms.context()
    DA = ClimaComms.array_type(ClimaComms.device(context))
    ncol = 8
    param_set = RRTMGP.default_parameters(FT)
    lw_file = get_lookup_filename(:gas, :lw)
    input_file = get_input_filename(:gas, :lw)
    lookup_lw, idx_gases = Dataset(lw_file, "r") do ds
        LookUpLW(ds, FT, DA)
    end
    ds_lw_in = Dataset(input_file, "r")
    (as, sfc_emis, _, _, _, _) = setup_clear_sky_as(
        context,
        ds_lw_in,
        idx_gases,
        1,
        lookup_lw,
        ncol,
        FT,
        VmrGM,
        param_set,
    )
    close(ds_lw_in)
    nlay, _ = AtmosphericStates.get_dims(as)
    gp = RRTMGP.RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)

    olr = map((1, 2, 4)) do n
        slv = NoScatLWRTE(
            gp;
            params = param_set,
            sfc_emis,
            inc_flux = nothing,
            n_gauss_angles = n,
        )
        solve_lw!(slv, as, lookup_lw)
        copy(Array(slv.flux.flux_up))
    end
    F1, F2, F4 = olr
    @test all(isfinite, F4)
    @test all(>(0), F4)
    @test F1 != F2 && F2 != F4
    @test maximum(abs, F2 .- F4) < maximum(abs, F1 .- F4)
    # the angular integral is a few W/m² correction, not a wholesale change
    @test maximum(abs, F1 .- F4) < 10
end

# Gray radiation is single-band and single-angle; asking for more must fail
# loudly rather than be silently ignored. The Layer-2 constructor knows the
# radiation method, so it catches this before a solve.
@testset "gray radiation rejects multiple angles" begin
    FT = Float64
    context = ClimaComms.context()
    DA = ClimaComms.array_type(ClimaComms.device(context))
    nlay, ncol = 20, 2
    params = RRTMGP.default_parameters(FT)
    gp = RRTMGP.RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
    as = AtmosphericStates.setup_gray_as_pr_grid(
        context,
        nlay,
        DA{FT}([-30, 30]),
        FT(1.0e5),
        FT(9.0e3),
        AtmosphericStates.GrayOpticalThicknessOGorman2008(FT),
        params,
        DA,
    )
    sfc_emis = fill!(DA{FT}(undef, 1, ncol), FT(1))
    cos_zenith = fill!(DA{FT}(undef, ncol), FT(0.5))
    toa = fill!(DA{FT}(undef, ncol), FT(1361))
    alb = fill!(DA{FT}(undef, 1, ncol), FT(0.2))
    bcs_lw = LwBCs(sfc_emis, nothing)
    bcs_sw = RRTMGP.BCs.SwBCs(cos_zenith, toa, alb, nothing, alb)

    # Layer 1: the state only arrives at solve time, so the check is there
    slv = NoScatLWRTE(
        gp;
        params,
        sfc_emis,
        inc_flux = nothing,
        n_gauss_angles = 2,
    )
    @test_throws ErrorException solve_lw!(slv, as)

    # Layer 2: the radiation method is known at construction
    @test_throws ErrorException RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.GrayRadiation(),
        params,
        bcs_lw,
        bcs_sw,
        as;
        op_lw = Optics.OneScalar(gp),
        op_sw = Optics.OneScalar(gp),
        n_gauss_angles = 2,
    )

    # The two-stream longwave solver (the default `op_lw`) has a fixed angular
    # treatment; asking for more angles must fail rather than be ignored. The
    # check precedes every use of the state and the lookups, so the mismatched
    # gray state never matters here.
    @test_throws ErrorException RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.ClearSkyRadiation(false),
        params,
        bcs_lw,
        bcs_sw,
        as;
        n_gauss_angles = 2,
    )

    # ... and the single-angle gray solver still runs through Layer 2, which
    # allocates no longwave band buffer for gray radiation
    solver = RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.GrayRadiation(),
        params,
        bcs_lw,
        bcs_sw,
        as;
        op_lw = Optics.OneScalar(gp),
        op_sw = Optics.OneScalar(gp),
    )
    @test solver.lws.fluxb === nothing
    RRTMGP.update_fluxes!(solver)
    @test all(isfinite, Array(RRTMGP.lw_flux_up(solver)))
    @test all(>(0), Array(RRTMGP.lw_flux_up(solver)))
end

nothing
