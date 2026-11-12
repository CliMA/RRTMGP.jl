using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import JET
using NCDatasets
using RRTMGP
using RRTMGP: RRTMGPGridParams, RRTMGPSolver
using RRTMGP.AtmosphericStates: AtmosphericState, CloudState, MaxRandomOverlap
using RRTMGP.BCs: LwBCs, SwBCs
using RRTMGP.LookUpTables: LookUpCld, LookUpPrecip
using RRTMGP.Optics:
    loc_lower,
    interp1d_equispaced,
    interp1d_loc_factor,
    compute_lookup_precip_props,
    add_cloud_optics_1scalar!,
    add_cloud_optics_2stream!
using RRTMGP.RTESolver: solve_lw!, solve_sw!
using RRTMGP.ArtifactPaths: get_lookup_filename

function _test_snow_solver_zero_alloc(solver, key)
    JET.@test_opt RRTMGP.update_fluxes!(solver, key)
    @allocated RRTMGP.update_fluxes!(solver, key)
    @test (@allocated RRTMGP.update_fluxes!(solver, key)) == 0
    return nothing
end

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

    # test compute_lookup_precip_props on a log-uniform grid
    nsize_log = 64
    r_lwr, r_upr = FT(5.0), FT(1100.0)
    r_grid = exp.(range(log(r_lwr), log(r_upr), length = nsize_log))
    # Table linear in log(r) is reconstructed to machine precision everywhere
    ext_tbl = FT(0.2) .- FT(0.02) .* log.(r_grid)
    ssa_tbl = fill(FT(0.9), nsize_log)
    asy_tbl = fill(FT(0.85), nsize_log)
    path = FT(10.0)

    # Zero water path returns zeros
    @test compute_lookup_precip_props(
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
    τ_low, τs_low, τsg_low = compute_lookup_precip_props(
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

    τ_high, τs_high, τsg_high = compute_lookup_precip_props(
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
    τ_mid, τs_mid, τsg_mid = compute_lookup_precip_props(
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

@testset "LookUpPrecip cloud optics kernels and solver integration" begin
    FT = Float64
    context = ClimaComms.context()
    DA = ClimaComms.array_type(ClimaComms.device(context))

    lkp_cld_lw = NCDataset(get_lookup_filename(:cloud, :lw), "r") do ds
        LookUpCld(ds, FT, DA)
    end
    lkp_cld_sw = NCDataset(get_lookup_filename(:cloud, :sw), "r") do ds
        LookUpCld(ds, FT, DA)
    end
    lkp_precip_lw = NCDataset(get_lookup_filename(:precip, :lw), "r") do ds
        LookUpPrecip(ds, FT, DA)
    end
    lkp_precip_sw = NCDataset(get_lookup_filename(:precip, :sw), "r") do ds
        LookUpPrecip(ds, FT, DA)
    end

    # 1. Direct kernel checks: add_cloud_optics_1scalar! & add_cloud_optics_2stream!
    nlay = 4
    cld_mask = trues(nlay)
    r_liq = fill(FT(10), nlay)
    r_ice = fill(FT(30), nlay)
    p_liq = zeros(FT, nlay)
    p_ice = zeros(FT, nlay)
    ice_rgh = 3
    r_sno = fill(FT(150), nlay) # > 90 µm cloud-ice ceiling
    p_sno = fill(FT(50), nlay)  # [g/m²]

    τ1 = fill(FT(0.1), nlay)
    @test_throws "LookUpPrecip" add_cloud_optics_1scalar!(
        τ1,
        cld_mask,
        r_liq,
        r_ice,
        p_liq,
        p_ice,
        ice_rgh,
        lkp_cld_lw,
        6;
        lkp_precip = nothing,
        cld_r_eff_sno = r_sno,
        cld_path_sno = p_sno,
    )
    add_cloud_optics_1scalar!(
        τ1,
        cld_mask,
        r_liq,
        r_ice,
        p_liq,
        p_ice,
        ice_rgh,
        lkp_cld_lw,
        6;
        lkp_precip = lkp_precip_lw,
        cld_r_eff_sno = r_sno,
        cld_path_sno = p_sno,
    )
    τ_150 = τ1[1]
    @test τ_150 > FT(0.1)

    # Increasing r_eff from 150 µm to 600 µm reduces optical depth (~1/r_eff)
    fill!(r_sno, FT(600))
    fill!(τ1, FT(0.1))
    add_cloud_optics_1scalar!(
        τ1,
        cld_mask,
        r_liq,
        r_ice,
        p_liq,
        p_ice,
        ice_rgh,
        lkp_cld_lw,
        6;
        lkp_precip = lkp_precip_lw,
        cld_r_eff_sno = r_sno,
        cld_path_sno = p_sno,
    )
    τ_600 = τ1[1]
    @test FT(0.1) < τ_600 < τ_150

    # Cloud mask gating: cld_mask == false suppresses snow contribution
    fill!(cld_mask, false)
    fill!(τ1, FT(0.1))
    add_cloud_optics_1scalar!(
        τ1,
        cld_mask,
        r_liq,
        r_ice,
        p_liq,
        p_ice,
        ice_rgh,
        lkp_cld_lw,
        6;
        lkp_precip = lkp_precip_lw,
        cld_r_eff_sno = r_sno,
        cld_path_sno = p_sno,
    )
    @test τ1[1] == FT(0.1)
    fill!(cld_mask, true)

    # TwoStream kernel (with and without delta scaling)
    τ2 = fill(FT(0.05), nlay)
    ssa2 = fill(FT(0.2), nlay)
    g2 = fill(FT(0.0), nlay)
    @test_throws "LookUpPrecip" add_cloud_optics_2stream!(
        τ2,
        ssa2,
        g2,
        cld_mask,
        r_liq,
        r_ice,
        p_liq,
        p_ice,
        ice_rgh,
        lkp_cld_sw,
        10;
        lkp_precip = nothing,
        cld_r_eff_sno = r_sno,
        cld_path_sno = p_sno,
        delta_scaling = true,
    )
    add_cloud_optics_2stream!(
        τ2,
        ssa2,
        g2,
        cld_mask,
        r_liq,
        r_ice,
        p_liq,
        p_ice,
        ice_rgh,
        lkp_cld_sw,
        10;
        lkp_precip = lkp_precip_sw,
        cld_r_eff_sno = r_sno,
        cld_path_sno = p_sno,
        delta_scaling = false,
    )
    @test τ2[1] > FT(0.05)
    @test FT(0.2) < ssa2[1] <= FT(1.0)
    @test FT(0.0) < g2[1] < FT(1.0)

    # 2. End-to-end solver integration (AllSkyRadiation & AllSkyRadiationWithClearSkyDiagnostics)
    nlay_atm = 30
    ncol = 2
    prof =
        RRTMGP.standard_atmosphere(FT; kind = :tropical, nlay = nlay_atm, ncol)
    gp = RRTMGPGridParams(FT; context, domain_nlay = nlay_atm, ncol)
    params = RRTMGP.default_parameters(FT)
    method_diag = RRTMGP.AllSkyRadiationWithClearSkyDiagnostics(false, true)
    method_allsky = RRTMGP.AllSkyRadiation(false, true)
    bundle = RRTMGP.lookup_tables(gp, method_diag)

    vmr_wm = zeros(FT, bundle.ngas_lw)
    for (name, val) in prof.well_mixed_vmr
        idx = get(bundle.idx_gases_lw, name, 0)
        idx > 0 && (vmr_wm[idx] = FT(val))
    end
    make_vmr() = RRTMGP.VolumeMixingRatios.VmrGM(
        DA{FT}(copy(prof.vmr_h2o)),
        DA{FT}(copy(prof.vmr_o3)),
        DA{FT}(copy(vmr_wm)),
    )

    make_cs() = CloudState(
        DA(fill(FT(10), nlay_atm, ncol)),
        DA(fill(FT(25), nlay_atm, ncol)),
        DA(zeros(FT, nlay_atm, ncol)),
        DA(zeros(FT, nlay_atm, ncol)),
        DA(zeros(FT, nlay_atm, ncol)),
        DA(zeros(FT, ncol)),
        DA(zeros(FT, ncol)),
        DA(falses(nlay_atm, ncol)),
        DA(falses(nlay_atm, ncol)),
        MaxRandomOverlap(),
        3;
        cld_r_eff_sno = DA(fill(FT(250), nlay_atm, ncol)),
        cld_path_sno = DA(zeros(FT, nlay_atm, ncol)),
    )
    function make_as()
        layerdata_h = zeros(FT, 4, nlay_atm, ncol)
        layerdata_h[2, :, :] .= prof.p_lay
        layerdata_h[3, :, :] .= prof.t_lay
        return AtmosphericState(
            DA(zeros(FT, ncol)),
            DA(copy(prof.lat)),
            DA(layerdata_h),
            DA(copy(prof.p_lev)),
            DA(copy(prof.t_lev)),
            DA(copy(prof.t_sfc)),
            make_vmr(),
            make_cs(),
            nothing,
        )
    end

    bcs_lw = LwBCs(DA(fill(FT(0.98), bundle.nbnd_lw, ncol)), nothing)
    bcs_sw = SwBCs(
        DA(fill(FT(0.6), ncol)),
        DA(fill(FT(1361), ncol)),
        DA(fill(FT(0.15), bundle.nbnd_sw, ncol)),
        nothing,
        DA(fill(FT(0.15), bundle.nbnd_sw, ncol)),
    )

    solver_diag = RRTMGPSolver(
        gp,
        method_diag,
        params,
        bcs_lw,
        bcs_sw,
        make_as();
        lookups = bundle,
    )
    solver_allsky = RRTMGPSolver(
        gp,
        method_allsky,
        params,
        bcs_lw,
        bcs_sw,
        make_as();
        lookups = bundle,
    )

    key = UInt32(20260417)
    RRTMGP.update_fluxes!(solver_diag, key)
    olr_nosnow = copy(Array(RRTMGP.lw_flux_up(solver_diag))[end, :])
    sw_up_nosnow = copy(Array(RRTMGP.sw_flux_up(solver_diag))[end, :])
    clr_lw_ref = copy(Array(RRTMGP.clear_lw_flux_up(solver_diag)))
    clr_sw_ref = copy(Array(RRTMGP.clear_sw_flux_up(solver_diag)))

    # Add upper-tropospheric snow layer (layers 18:24) in overcast conditions
    for slv in (solver_diag, solver_allsky)
        cf = RRTMGP.cloud_fraction(slv)
        swp = RRTMGP.cloud_snow_water_path(slv)
        sre = RRTMGP.cloud_snow_effective_radius(slv)
        fill!(cf, FT(0))
        fill!(swp, FT(0))
        fill!(sre, FT(250))
        Array_cf = Array(cf)
        Array_swp = Array(swp)
        Array_cf[18:24, :] .= FT(1)
        Array_swp[18:24, :] .= FT(120)
        copyto!(cf, Array_cf)
        copyto!(swp, Array_swp)
    end

    RRTMGP.update_fluxes!(solver_diag, key)
    RRTMGP.update_fluxes!(solver_allsky, key)

    # Upper-tropospheric snow reduces OLR and increases reflected TOA SW
    @test all(Array(RRTMGP.lw_flux_up(solver_diag))[end, :] .< olr_nosnow)
    @test all(Array(RRTMGP.sw_flux_up(solver_diag))[end, :] .> sw_up_nosnow)
    # Clear-sky diagnostics are unaffected by snow
    @test Array(RRTMGP.clear_lw_flux_up(solver_diag)) ≈ clr_lw_ref
    @test Array(RRTMGP.clear_sw_flux_up(solver_diag)) ≈ clr_sw_ref
    # Fused both-skies path matches standard AllSkyRadiation path
    @test Array(RRTMGP.lw_flux_net(solver_diag)) ≈
          Array(RRTMGP.lw_flux_net(solver_allsky))
    @test Array(RRTMGP.sw_flux_net(solver_diag)) ≈
          Array(RRTMGP.sw_flux_net(solver_allsky))

    # Zero-allocation and type-stability checks with active snow
    _test_snow_solver_zero_alloc(solver_diag, key)
    _test_snow_solver_zero_alloc(solver_allsky, key)

    # Functional-core error check when cld_path_sno !== nothing and lookup_*_precip === nothing
    @test_throws "LookUpPrecip" solve_lw!(
        solver_allsky.lws,
        solver_allsky.as,
        bundle.lookup_lw,
        bundle.lookup_lw_cld,
        nothing,
        nothing,
        nothing,
    )
    @test_throws "LookUpPrecip" solve_sw!(
        solver_allsky.sws,
        solver_allsky.as,
        bundle.lookup_sw,
        bundle.lookup_sw_cld,
        nothing,
        nothing,
        nothing,
    )
end
