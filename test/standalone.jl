using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
import Adapt
import JET

# A CPU adaptor that forces a fresh copy of each array. It lets us verify — without a
# GPU — that RRTMGP's custom `adapt_structure` overloads rebuild the optics/sources
# scratch as views into ONE newly-adapted buffer instead of de-aliasing them into
# separate arrays (which a plain `@adapt_structure` would do on a real device transfer).
# This is a memory-topology / zero-allocation contract; results are correct either way
# because the scratch is recomputed on every solve.
struct CopyToArray end
Adapt.adapt_storage(::CopyToArray, x::AbstractArray) = copy(x)

@testset "gray standalone solve_gray" begin
    for FT in (Float32, Float64)
        nlay, ncol = 60, 9
        out = RRTMGP.solve_gray(FT; nlay, ncol)
        @test out isa RRTMGP.RadiationOutput
        @test out.solver isa RRTMGP.RRTMGPSolver
        @test eltype(Array(out.net)) == FT
        @test size(Array(out.net)) == (nlay + 1, ncol)
        @test size(Array(out.heating_rate)) == (nlay, ncol)
        for f in (out.lw_net, out.sw_net, out.net, out.heating_rate)
            @test all(isfinite, Array(f))
        end
        # net flux is the sum of the longwave and shortwave net fluxes
        @test Array(out.net) ≈ Array(out.lw_net) .+ Array(out.sw_net)
        # the sun is up, so there is downwelling shortwave flux
        @test maximum(Array(out.sw_dn)) > 0
    end
end

# prepare_atmosphere! is the preparation half of update_fluxes!: running it
# separately and then the full update must give the same fluxes as the full
# update alone (the cascade is idempotent for an already-consistent state).
@testset "prepare_atmosphere! (prepare/solve split)" begin
    solver = RRTMGP.solve_gray(Float64; nlay = 20, ncol = 4).solver
    RRTMGP.prepare_atmosphere!(solver) # callable on its own
    RRTMGP.update_fluxes!(solver)
    F1 = copy(Array(RRTMGP.net_flux(solver)))
    RRTMGP.prepare_atmosphere!(solver) # prepare again, then full update
    RRTMGP.update_fluxes!(solver)
    @test Array(RRTMGP.net_flux(solver)) == F1
end

# Renamed-module compat aliases (deprecated names, kept for one release).
@testset "renamed-module aliases" begin
    @test RRTMGP.Vmrs === RRTMGP.VolumeMixingRatios
    @test RRTMGP.GrayUtils === RRTMGP.GrayAtmosphere
end

@testset "default_parameters" begin
    p = RRTMGP.default_parameters(Float64)
    @test RRTMGP.Parameters.grav(p) == 9.81
    @test RRTMGP.Parameters.cp_d(p) ≈ RRTMGP.Parameters.R_d(p) / (2 / 7)
end

@testset "non-gray needs NCDatasets" begin
    @test RRTMGP._check_lookup_support(RRTMGP.GrayRadiation()) === nothing
    # Only assert the actionable error when NCDatasets is genuinely not loaded
    # (the full test suite loads it, which makes the spectral methods available).
    if isnothing(Base.get_extension(RRTMGP, :RRTMGPNCDatasetsExt))
        @test_throws ErrorException RRTMGP._check_lookup_support(
            RRTMGP.ClearSkyRadiation(false),
        )
    end
end

@testset "standard_atmosphere profiles" begin
    for FT in (Float32, Float64)
        nlay, ncol = 40, 3
        prof = RRTMGP.standard_atmosphere(FT; kind = :tropical, nlay, ncol)
        @test prof isa RRTMGP.AtmosphereProfile
        @test eltype(prof.p_lay) == FT
        @test size(prof.p_lay) == (nlay, ncol)
        @test size(prof.p_lev) == (nlay + 1, ncol)
        @test size(prof.t_sfc) == (ncol,)
        # all columns are identical
        @test prof.t_lay[:, 1] == prof.t_lay[:, end]
        # pressure decreases monotonically upward, from the surface value
        @test prof.p_lev[1, 1] ≈ 101325.0
        @test all(diff(prof.p_lev[:, 1]) .< 0)
        @test all(diff(prof.p_lay[:, 1]) .< 0)
        # layers sit between their bounding levels
        @test all(prof.p_lev[2:end, 1] .< prof.p_lay[:, 1] .<= prof.p_lev[1:(end - 1), 1])
        # surface temperature and the tropical tropopause
        @test prof.t_lev[1, 1] == 300
        @test 180 < minimum(prof.t_lev) < 210
        # water vapor decays to the stratospheric floor; ozone peaks aloft
        @test prof.vmr_h2o[1, 1] > 1e-2
        @test minimum(prof.vmr_h2o) ≈ 4.0e-6
        @test 20e3 <
              prof.z_lev[argmax(prof.vmr_o3[:, 1]), 1] <
              40e3
        @test prof.well_mixed_vmr["co2"] ≈ 420e-6
    end
    # the three kinds are ordered warm → cold at the surface
    t_sfc(kind) = RRTMGP.standard_atmosphere(Float64; kind).t_sfc[1]
    @test t_sfc(:tropical) > t_sfc(:midlatitude_summer) > t_sfc(:subarctic_winter)
    @test_throws ErrorException RRTMGP.standard_atmosphere(Float64; kind = :venus)
end

@testset "gray solve(profile)" begin
    for FT in (Float32, Float64)
        nlay, ncol = 30, 2
        prof = RRTMGP.standard_atmosphere(FT; nlay, ncol)
        out = RRTMGP.solve(prof; method = RRTMGP.GrayRadiation())
        @test out isa RRTMGP.RadiationOutput
        @test eltype(Array(out.net)) == FT
        @test size(Array(out.net)) == (nlay + 1, ncol)
        @test size(Array(out.heating_rate)) == (nlay, ncol)
        for f in (out.lw_up, out.lw_dn, out.sw_up, out.sw_dn, out.net)
            @test all(isfinite, Array(f))
        end
        @test Array(out.net) ≈ Array(out.lw_net) .+ Array(out.sw_net)
        # OLR is positive and bounded by σT⁴ of the warmest layer
        @test 0 < Array(out.lw_up)[end, 1] < 600
        # the profile is not aliased by the solve (arrays stay host-side, unclipped)
        @test prof.p_lay isa Array
    end
    # cloud/aerosol methods are rejected with actionable errors
    prof = RRTMGP.standard_atmosphere(Float64; nlay = 10)
    @test_throws ErrorException RRTMGP.solve(
        prof;
        method = RRTMGP.AllSkyRadiation(false, false),
    )
    @test_throws ErrorException RRTMGP.solve(
        prof;
        method = RRTMGP.ClearSkyRadiation(true),
    )
end

# Passing `op_sw = OneScalar` must route the constructor to a non-scattering
# `NoScatSWRTE` solver (rather than the default `TwoStreamSWRTE`), and the solve
# must run end to end. Uses gray radiation so no lookup tables are needed.
@testset "gray NoScat shortwave solver (op_sw = OneScalar)" begin
    context = ClimaComms.context()
    DA = ClimaComms.array_type(ClimaComms.device(context))
    for FT in (Float32, Float64)
        nlay, ncol = 40, 4
        gp = RRTMGP.RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
        params = RRTMGP.default_parameters(FT)
        lat = DA{FT}(range(FT(-80), FT(80); length = ncol))
        as = RRTMGP.AtmosphericStates.setup_gray_as_pr_grid(
            context,
            nlay,
            lat,
            FT(1.0e5),
            FT(9.0e3),
            RRTMGP.AtmosphericStates.GrayOpticalThicknessOGorman2008(FT),
            params,
            DA,
        )
        sfc_emis = fill!(DA{FT}(undef, 1, ncol), FT(1))
        bcs_lw = RRTMGP.BCs.LwBCs(sfc_emis, nothing)
        cos_zen = fill!(DA{FT}(undef, ncol), FT(0.5))
        toa = fill!(DA{FT}(undef, ncol), FT(1361))
        alb_dir = fill!(DA{FT}(undef, 1, ncol), FT(0.2))
        alb_dif = fill!(DA{FT}(undef, 1, ncol), FT(0.2))
        bcs_sw = RRTMGP.BCs.SwBCs(cos_zen, toa, alb_dir, nothing, alb_dif)

        op_sw = RRTMGP.Optics.OneScalar(gp)
        solver = RRTMGP.RRTMGPSolver(
            gp,
            RRTMGP.GrayRadiation(),
            params,
            bcs_lw,
            bcs_sw,
            as;
            op_sw,
        )
        @test solver.sws isa RRTMGP.RTE.NoScatSWRTE

        RRTMGP.update_fluxes!(solver)
        sw_net = Array(RRTMGP.sw_flux_net(solver))
        sw_dn = Array(RRTMGP.sw_flux_dn(solver))
        @test size(sw_net) == (nlay + 1, ncol)
        @test all(isfinite, sw_net)
        @test all(isfinite, sw_dn)
        @test maximum(sw_dn) > 0 # the sun is up, so there is downwelling shortwave

        # Non-scattering shortwave optics are rejected for spectral (non-gray)
        # radiation. The guard fires before `as` is inspected, so reusing the gray
        # state here is fine.
        @test_throws ErrorException RRTMGP.RRTMGPSolver(
            gp,
            RRTMGP.ClearSkyRadiation(false),
            params,
            bcs_lw,
            bcs_sw,
            as;
            op_sw,
        )
    end
end

# With an isothermal boundary layer, RRTMGP's internal arrays carry an extra top
# layer/level, but the getters must return domain-sized (masked) views so a
# caller sees arrays matching the physical domain.
@testset "gray solver with isothermal boundary layer" begin
    context = ClimaComms.context()
    DA = ClimaComms.array_type(ClimaComms.device(context))
    for FT in (Float32, Float64)
        domain_nlay, ncol = 40, 4
        gp = RRTMGP.RRTMGPGridParams(
            FT;
            context,
            domain_nlay,
            ncol,
            isothermal_boundary_layer = true,
        )
        nlay = gp.nlay # total (boundary-extended) layer count
        params = RRTMGP.default_parameters(FT)
        lat = DA{FT}(range(FT(-80), FT(80); length = ncol))
        as = RRTMGP.AtmosphericStates.setup_gray_as_pr_grid(
            context,
            nlay,
            lat,
            FT(1.0e5),
            FT(9.0e3),
            RRTMGP.AtmosphericStates.GrayOpticalThicknessOGorman2008(FT),
            params,
            DA,
        )
        sfc_emis = fill!(DA{FT}(undef, 1, ncol), FT(1))
        bcs_lw = RRTMGP.BCs.LwBCs(sfc_emis, nothing)
        cos_zen = fill!(DA{FT}(undef, ncol), FT(0.5))
        toa = fill!(DA{FT}(undef, ncol), FT(1361))
        alb = fill!(DA{FT}(undef, 1, ncol), FT(0.2))
        bcs_sw = RRTMGP.BCs.SwBCs(cos_zen, toa, alb, nothing, copy(alb))
        solver = RRTMGP.RRTMGPSolver(
            gp,
            RRTMGP.GrayRadiation(),
            params,
            bcs_lw,
            bcs_sw,
            as,
        )

        RRTMGP.update_fluxes!(solver)

        # the net-flux buffer is boundary-extended, host-facing (nlev, ncol); getters are domain-masked
        @test size(solver.net_flux_buffer, 1) == nlay + 1             # full levels
        @test size(RRTMGP.net_flux(solver), 1) == domain_nlay + 1     # domain levels
        @test size(RRTMGP.level_pressure(solver), 1) == domain_nlay + 1
        @test size(RRTMGP.lw_flux_net(solver), 1) == domain_nlay + 1
        @test size(RRTMGP.layer_temperature(solver), 1) == domain_nlay # domain layers
        @test size(RRTMGP.layer_pressure(solver), 1) == domain_nlay
        @test size(RRTMGP.heating_rate(solver)) == (domain_nlay, ncol)

        net = Array(RRTMGP.net_flux(solver))
        @test net ≈
              Array(RRTMGP.lw_flux_net(solver)) .+
              Array(RRTMGP.sw_flux_net(solver))
        @test all(isfinite, net)
        @test all(isfinite, Array(RRTMGP.heating_rate(solver)))
    end
end

# The optics and source scratch buffers store `τ`/`ssa`/`g`/`lev_source`/`albedo`/`src`
# as views into a single packed `layerdata`/`leveldata` buffer. Custom `adapt_structure`
# overloads must preserve that topology across a device transfer (adapt the buffer once,
# rebuild the views into it) so a checkpoint restore stays zero-allocation.
@testset "Adapt preserves optics/sources view topology" begin
    for FT in (Float32, Float64)
        nlay, ncol = 40, 4
        solver = RRTMGP.solve_gray(FT; nlay, ncol).solver # two-stream lw + sw
        F_before = copy(RRTMGP.net_flux(solver))

        solver2 = Adapt.adapt(CopyToArray(), solver)

        # the adapted buffers are genuinely fresh copies, not the originals
        @test solver2.lws.op.layerdata !== solver.lws.op.layerdata
        # every scratch view aliases its own solver's single adapted buffer
        @test parent(solver2.lws.op.τ) === solver2.lws.op.layerdata
        @test parent(solver2.lws.op.ssa) === solver2.lws.op.layerdata
        @test parent(solver2.lws.op.g) === solver2.lws.op.layerdata
        @test parent(solver2.sws.op.τ) === solver2.sws.op.layerdata
        @test parent(solver2.lws.src.lev_source) === solver2.lws.src.leveldata
        @test parent(solver2.lws.src.albedo) === solver2.lws.src.leveldata
        @test parent(solver2.lws.src.src) === solver2.lws.src.leveldata
        @test parent(solver2.sws.src.albedo) === solver2.sws.src.leveldata
        @test parent(solver2.sws.src.src) === solver2.sws.src.leveldata

        # functional: the adapted solver still solves and reproduces the fluxes
        RRTMGP.update_fluxes!(solver2)
        @test RRTMGP.net_flux(solver2) ≈ F_before
    end

    # non-scattering optics (single `τ` view) used by the gray/no-scatter solvers
    gp = RRTMGP.RRTMGPGridParams(
        Float64;
        context = ClimaComms.context(),
        domain_nlay = 10,
        ncol = 3,
    )
    op1 = Adapt.adapt(CopyToArray(), RRTMGP.Optics.OneScalar(gp))
    @test parent(op1.τ) === op1.layerdata
end

# Opt-in input validation: `RRTMGP.check_values[] = true` makes `update_fluxes!`
# validate solver inputs; off (the default) it costs a single branch — the
# zero-alloc/JET testset below runs with the toggle off and still asserts 0.
@testset "input validation toggle" begin
    solver = RRTMGP.solve_gray(Float64; nlay = 12, ncol = 3).solver
    RRTMGP.validate_inputs(solver) # valid inputs pass
    RRTMGP.cos_zenith(solver) .= 2 # corrupt an input (broadcast: device-safe)
    @test_throws ErrorException RRTMGP.validate_inputs(solver)
    try
        RRTMGP.check_values[] = true
        @test_throws ErrorException RRTMGP.update_fluxes!(solver)
    finally
        RRTMGP.check_values[] = false
    end
    RRTMGP.cos_zenith(solver) .= 1 / 2
    RRTMGP.update_fluxes!(solver) # runs again once repaired, toggle off
    @test all(isfinite, Array(RRTMGP.net_flux(solver)))
end

# Layer-2 `update_fluxes!` runs every radiation step, so it must be allocation-free
# and type-stable. The getters are type-stable because `_domain_view` always returns
# one concrete view type (no `Union{Array, SubArray}` from a runtime boundary-layer
# branch), and `update_fluxes!`/`update_net_fluxes!` mutate in place and return
# `nothing` rather than an escaping view. (`@allocated`/host-side JET are CPU-only.)
@testset "update_fluxes! zero-alloc + type-stable (Layer 2)" begin
    context = ClimaComms.context()
    solver = RRTMGP.solve_gray(Float64; nlay = 60, ncol = 10).solver
    # getters infer to a single concrete view type (device-independent)
    @test (@inferred RRTMGP.net_flux(solver)) isa SubArray
    @test (@inferred RRTMGP.lw_flux_net(solver)) isa SubArray
    RRTMGP.update_fluxes!(solver)        # warm up / compile
    RRTMGP.update_fluxes!(solver, 1234)
    # Single-threaded only: multi-threaded `@threaded` loops allocate task state, so the
    # zero-alloc contract (like the kernel tests) is asserted on `CPUSingleThreaded`.
    if ClimaComms.device(context) isa ClimaComms.CPUSingleThreaded
        JET.@test_opt RRTMGP.update_fluxes!(solver)
        JET.@test_opt RRTMGP.update_net_fluxes!(solver)
        # `@allocated`'s first measurement of a call includes one-time setup; discard it,
        # then assert the steady-state (what a host pays per radiation step) is zero.
        @allocated RRTMGP.update_fluxes!(solver)
        @test (@allocated RRTMGP.update_fluxes!(solver)) == 0
        @allocated RRTMGP.update_fluxes!(solver, 1234)
        @test (@allocated RRTMGP.update_fluxes!(solver, 1234)) == 0
        @allocated RRTMGP.update_net_fluxes!(solver)
        @test (@allocated RRTMGP.update_net_fluxes!(solver)) == 0
    end
end
