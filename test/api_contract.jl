using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
import RRTMGP.VolumeMixingRatios: Vmr, VmrGM
import Serialization

# Unit tests for the Layer-2 API contract: constructor validation, the
# informative getter-guard errors, input validation, boundary conditions, and
# the physical laws behind heating_rate and the deep-atmosphere flux scaling.

# Small gray problem shared by these tests (the gray path needs no lookup
# tables, so every case here runs after a bare `using RRTMGP`).
function _gray_pieces(FT; nlay = 20, ncol = 4)
    context = ClimaComms.context()
    DA = ClimaComms.array_type(ClimaComms.device(context))
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
    cos_zen = fill!(DA{FT}(undef, ncol), FT(0.5))
    toa = fill!(DA{FT}(undef, ncol), FT(1361))
    alb_dir = fill!(DA{FT}(undef, 1, ncol), FT(0.2))
    alb_dif = fill!(DA{FT}(undef, 1, ncol), FT(0.2))
    bcs_sw = RRTMGP.BCs.SwBCs(cos_zen, toa, alb_dir, nothing, alb_dif)
    return (; context, DA, gp, params, as, sfc_emis, bcs_sw)
end

function _gray_solver(FT; nlay = 20, ncol = 4, inc_flux = nothing, kwargs...)
    (; gp, params, as, sfc_emis, bcs_sw) = _gray_pieces(FT; nlay, ncol)
    bcs_lw = RRTMGP.BCs.LwBCs(sfc_emis, inc_flux)
    return RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.GrayRadiation(),
        params,
        bcs_lw,
        bcs_sw,
        as;
        kwargs...,
    )
end

@testset "RRTMGPSolver constructor validation" begin
    (; gp, params, as, sfc_emis, bcs_sw) = _gray_pieces(Float64)
    bcs_lw = RRTMGP.BCs.LwBCs(sfc_emis, nothing)

    # z-based interpolation/extrapolation without altitudes fails at
    # construction (not with a MethodError mid-solve)
    @test_throws "requires altitudes" RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.GrayRadiation(),
        params,
        bcs_lw,
        bcs_sw,
        as;
        interpolation = RRTMGP.BestFit(),
    )
    @test_throws "requires altitudes" RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.GrayRadiation(),
        params,
        bcs_lw,
        bcs_sw,
        as;
        bottom_extrapolation = RRTMGP.HydrostaticBottom(),
    )

    # per-band fluxes are meaningless for a single-band (gray) atmosphere
    @test_throws "not supported for GrayRadiation" RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.GrayRadiation(),
        params,
        bcs_lw,
        bcs_sw,
        as;
        spectral_fluxes = true,
    )

    # ... and require two-stream optics for both bands. The guard fires before
    # the lookup tables or atmospheric state are inspected, so a placeholder
    # bundle and the gray state suffice (no NetCDF read).
    placeholder = RRTMGP.LookupBundle(; nbnd_lw = 16, nbnd_sw = 14)
    @test_throws "requires two-stream optics" RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.ClearSkyRadiation(false),
        params,
        bcs_lw,
        bcs_sw,
        as;
        op_lw = RRTMGP.Optics.OneScalar(gp),
        lookups = placeholder,
        spectral_fluxes = true,
    )
end

@testset "getter guards raise informative errors (gray solver)" begin
    solver = _gray_solver(Float64)

    # optional workspaces that this solver was not built with
    @test_throws "clear-sky diagnostic fluxes" RRTMGP.clear_net_flux(solver)
    @test_throws "clear-sky diagnostic fluxes" RRTMGP.clear_lw_flux_up(solver)
    @test_throws "clear-sky diagnostic fluxes" RRTMGP.clear_sw_flux_dn(solver)
    @test_throws "spectral lookup tables are not present" RRTMGP.lw_band_bounds(
        solver,
    )
    @test_throws "spectral lookup tables are not present" RRTMGP.sw_band_bounds(
        solver,
    )
    @test_throws "aerosol lookup tables are not present" RRTMGP.aerosol_radius(
        solver,
        "dust1",
    )
    @test_throws "aerosol lookup tables are not present" RRTMGP.aerosol_column_mass_density(
        solver,
        "dust1",
    )
    @test_throws "spectral fluxes were not retained" RRTMGP.spectral_lw_flux_up(
        solver,
    )
    @test_throws "spectral fluxes were not retained" RRTMGP.spectral_sw_flux_dn(
        solver,
    )

    # configuration getters on the same solver
    @test RRTMGP.radiation_method(solver) isa RRTMGP.GrayRadiation
    @test RRTMGP.optical_thickness_parameter(solver) isa
          RRTMGP.AtmosphericStates.GrayOpticalThicknessOGorman2008
    @test RRTMGP.isothermal_boundary_layer(solver) == false
    @test RRTMGP.top_of_atmosphere_lw_flux_dn(solver) === nothing
    @test RRTMGP.top_of_atmosphere_diffuse_sw_flux_dn(solver) === nothing
end

@testset "validate_inputs names the offending field" begin
    # each case corrupts one input on a fresh solver and asserts the error
    # message identifies it
    cases = [
        ("level_pressure", s -> RRTMGP.level_pressure(s) .= -1),
        ("level_temperature", s -> RRTMGP.level_temperature(s) .= NaN),
        ("layer_pressure", s -> RRTMGP.layer_pressure(s) .= 0),
        ("layer_temperature", s -> RRTMGP.layer_temperature(s) .= -5),
        ("surface_temperature", s -> RRTMGP.surface_temperature(s) .= 0),
        ("toa_flux", s -> RRTMGP.toa_flux(s) .= -1),
        ("surface_emissivity", s -> RRTMGP.surface_emissivity(s) .= 2),
        (
            "direct_sw_surface_albedo",
            s -> RRTMGP.direct_sw_surface_albedo(s) .= 2,
        ),
        (
            "diffuse_sw_surface_albedo",
            s -> RRTMGP.diffuse_sw_surface_albedo(s) .= -1,
        ),
    ]
    for (name, corrupt!) in cases
        solver = _gray_solver(Float64; nlay = 10, ncol = 2)
        RRTMGP.validate_inputs(solver) # sane inputs pass
        corrupt!(solver)
        @test_throws name RRTMGP.validate_inputs(solver)
    end
end

@testset "_check_vmr rejects negative / non-finite mixing ratios" begin
    vmr_gm = VmrGM(zeros(2, 2), zeros(2, 2), zeros(3))
    @test RRTMGP._check_vmr(vmr_gm) == true
    vmr_gm.vmr_h2o[1] = -1
    @test_throws "vmr_h2o" RRTMGP._check_vmr(vmr_gm)
    vmr_gm.vmr_h2o[1] = 0
    vmr_gm.vmr_o3[2] = NaN
    @test_throws "vmr_o3" RRTMGP._check_vmr(vmr_gm)

    vmr = Vmr(zeros(3, 2, 2))
    @test RRTMGP._check_vmr(vmr) == true
    vmr.vmr[1] = -1
    @test_throws "vmr" RRTMGP._check_vmr(vmr)

    @test RRTMGP._check_vmr(nothing) == true
end

# Physical law: the downwelling longwave flux at the top level equals the
# prescribed incident flux (the LW top boundary condition); with no incident
# flux it is zero. Also pins the (ncol, ngpt) → physical-dim-first transpose
# contract of the TOA getters.
@testset "incident longwave flux boundary condition" begin
    for FT in (Float32, Float64)
        pieces = _gray_pieces(FT)
        ncol = size(pieces.sfc_emis, 2)

        solver0 = _gray_solver(FT)
        RRTMGP.update_fluxes!(solver0)
        @test all(iszero, Array(RRTMGP.lw_flux_dn(solver0))[end, :])

        inc = fill!(pieces.DA{FT}(undef, ncol, 1), FT(25))
        solver = _gray_solver(FT; inc_flux = inc)
        toa_dn = RRTMGP.top_of_atmosphere_lw_flux_dn(solver)
        @test size(toa_dn) == (1, ncol)
        @test all(==(FT(25)), Array(toa_dn))

        RRTMGP.update_fluxes!(solver)
        @test all(Array(RRTMGP.lw_flux_dn(solver))[end, :] .≈ FT(25))
        # the incident flux warms the surface's incoming longwave radiation
        @test all(
            Array(RRTMGP.lw_flux_dn(solver))[1, :] .>
            Array(RRTMGP.lw_flux_dn(solver0))[1, :],
        )
    end
end

# Physical law: `deep_atmosphere_inverse_scaling` multiplies the fluxes
# pointwise (geometric expansion of grid columns aloft), so a uniform factor c
# must scale every flux by exactly c relative to the unscaled solve.
@testset "deep-atmosphere metric scaling multiplies fluxes" begin
    for FT in (Float32, Float64)
        nlay, ncol = 20, 4
        base = _gray_solver(FT; nlay, ncol)
        RRTMGP.update_fluxes!(base)

        pieces = _gray_pieces(FT; nlay, ncol)
        scaling = fill!(pieces.DA{FT}(undef, nlay + 1, ncol), FT(2))
        bcs_lw = RRTMGP.BCs.LwBCs(pieces.sfc_emis, nothing)
        scaled = RRTMGP.RRTMGPSolver(
            pieces.gp,
            RRTMGP.GrayRadiation(),
            pieces.params,
            bcs_lw,
            pieces.bcs_sw,
            pieces.as;
            deep_atmosphere_inverse_scaling = scaling,
        )
        @test RRTMGP.deep_atmosphere_inverse_scaling(scaled) !== nothing
        RRTMGP.update_fluxes!(scaled)

        for getter in (
            RRTMGP.lw_flux_up,
            RRTMGP.lw_flux_dn,
            RRTMGP.lw_flux_net,
            RRTMGP.sw_flux_up,
            RRTMGP.sw_flux_dn,
            RRTMGP.sw_direct_flux_dn,
            RRTMGP.sw_flux_net,
            RRTMGP.net_flux,
        )
            @test Array(getter(scaled)) ≈ 2 .* Array(getter(base))
        end
    end
end

# Physical law: the radiative heating rate is the net-flux divergence in
# pressure coordinates, hr = (g/cₚ) ∂F_net/∂p. Recompute the finite difference
# from the getters and compare against heating_rate(s).
@testset "heating_rate matches the net-flux divergence" begin
    for FT in (Float32, Float64)
        out = RRTMGP.solve_gray(FT; nlay = 20, ncol = 3)
        params = out.solver.params
        g = RRTMGP.Parameters.grav(params)
        cₚ = RRTMGP.Parameters.cp_d(params)
        F = Array(out.net)
        p = Array(RRTMGP.level_pressure(out.solver))
        expected =
            (g / cₚ) .* (F[2:end, :] .- F[1:(end - 1), :]) ./
            (p[2:end, :] .- p[1:(end - 1), :])
        @test Array(out.heating_rate) ≈ expected
    end
end

@testset "gray clip! floors pressures and leaves temperatures alone" begin
    FT = Float64
    nlay, ncol = 3, 2
    p_lay = FT[50 60; 500 510; 5000 5010]
    p_lev = FT[10 20; 100 110; 1000 1010; 10000 10010]
    t_lay = fill(FT(150), nlay, ncol)
    t_lev = fill(FT(150), nlay + 1, ncol)
    z_lev = zeros(FT, nlay + 1, ncol)
    as = RRTMGP.AtmosphericStates.GrayAtmosphericState(
        zeros(FT, ncol),
        p_lay,
        p_lev,
        t_lay,
        t_lev,
        z_lev,
        fill(FT(300), ncol),
        RRTMGP.AtmosphericStates.GrayOpticalThicknessOGorman2008(FT),
    )
    p_min = FT(100)
    # the temperature clamp bounds are ignored for the gray state: its analytic
    # optics remain valid outside the lookup-table range
    @test RRTMGP.clip!(as, p_min; t_min = FT(160), t_max = FT(355)) === as
    @test all(≥(p_min), as.p_lay)
    @test all(≥(p_min), as.p_lev)
    @test as.p_lay[3, 1] == 5000 # values above the floor are untouched
    @test all(==(FT(150)), as.t_lay)
    @test all(==(FT(150)), as.t_lev)

    # update_concentrations! is a no-op for gray radiation
    device = ClimaComms.device(ClimaComms.context())
    params = RRTMGP.default_parameters(FT)
    @test RRTMGP.update_concentrations!(as, params, device) === as
end
