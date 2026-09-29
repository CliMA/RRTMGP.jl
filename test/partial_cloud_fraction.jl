using Random
@isdefined(setup_cloudy_sky_as) || include("cloudy_sky_utils.jl")

"""
    partial_cloud_fraction_test(context, ::Type{FT}; ncol = 128)

Test partial cloud fraction with McICA stochastic cloud masking.

When `cld_frac < 1`, `build_cloud_mask!` samples the cloud mask stochastically
(maximum-random overlap). The sample is a deterministic function of
`(cloud_state.seed, column, g-point, layer)`, so it is selected by the key rather
than drawn from an RNG's state. This test verifies:
1. Fully overcast (`cld_frac = 1`) is deterministic — the mask is always
   `true` regardless of the key.
2. One key reproduces the same fluxes exactly, on CPU and on GPU alike.
3. A different key resamples, so successive radiation steps do not repeat one
   cloud field.
4. The same properties hold for spatially varying (per-layer, per-column)
   cloud fractions.
5. Diagnosed cloud cover lies in [0, 1].

Points 2 and 3 used to be asserted the other way around, and point 2 only on the
CPU: reproducibility came from `Random.seed!`, which cannot reach the device RNG
the GPU sampler drew from (RRTMGP.jl#316).
"""
function partial_cloud_fraction_test(
    context,
    ::Type{FT};
    ncol = 128,
) where {FT <: AbstractFloat}
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    overrides =
        (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)

    # ----- Read lookup tables (done once) -----
    lookup_lw, idx_gases = Dataset(get_lookup_filename(:gas, :lw), "r") do ds
        LookUpLW(ds, FT, DA)
    end
    lookup_sw, _ = Dataset(get_lookup_filename(:gas, :sw), "r") do ds
        LookUpSW(ds, FT, DA)
    end
    lookup_lw_cld = Dataset(get_lookup_filename(:cloud, :lw), "r") do ds
        LookUpCld(ds, FT, DA)
    end
    lookup_sw_cld = Dataset(get_lookup_filename(:cloud, :sw), "r") do ds
        LookUpCld(ds, FT, DA)
    end

    # ----- Set up atmospheric state with fully overcast clouds -----
    ds_in = Dataset(get_input_filename(:gas_clouds, :lw), "r")
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, cos_zenith, toa_flux, _ =
        setup_cloudy_sky_as(
            context,
            ds_in,
            idx_gases,
            lookup_lw,
            lookup_sw,
            lookup_lw_cld,
            lookup_sw_cld,
            FT(1),
            ncol,
            FT,
            param_set,
        )
    close(ds_in)

    # Base cloud fraction pattern: 1 where cloudy, 0 elsewhere
    base_cld_frac = copy(as.cloud_state.cld_frac)

    nlay, ncol = AtmosphericStates.get_dims(as)
    grid_params = RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)

    # ----- Create solvers -----
    slv_lw = TwoStreamLWRTE(
        grid_params;
        params = param_set,
        sfc_emis,
        inc_flux = nothing,
    )
    swbcs = (;
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse = nothing,
        sfc_alb_diffuse,
    )
    slv_sw = TwoStreamSWRTE(grid_params; swbcs...)

    color = :cyan
    printstyled(
        "Partial cloud fraction test with ncol = $ncol, FT = $FT\n",
        color = color,
    )
    printstyled("device = $device\n\n", color = color)

    # Helper: set cloud fraction, optionally seed RNG, run solvers,
    # and return copies of the result arrays.
    function run_solve!(cld_frac_values; key = nothing)
        as.cloud_state.cld_frac .= cld_frac_values
        # The McICA sample is a function of (key, column, g-point, layer), so
        # the key selects it. `Random.seed!` no longer has any bearing on it,
        # which is the point of the change -- it never did on the GPU.
        isnothing(key) || fill!(as.cloud_state.seed, UInt32(key))
        solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, nothing, nothing)
        solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld, nothing, nothing)
        return (
            lw_net = copy(Array(PermutedDimsArray(slv_lw.flux.flux_net, (2, 1)))),
            sw_net = copy(Array(PermutedDimsArray(slv_sw.flux.flux_net, (2, 1)))),
            cld_cover_lw = copy(Array(as.cloud_state.cld_cover_lw)),
            cld_cover_sw = copy(Array(as.cloud_state.cld_cover_sw)),
        )
    end

    key = 1234

    # ------------------------------------------------------------------
    # 1. Fully overcast (cldfrac = 1): deterministic
    #    The mask is always `true` regardless of random draws.
    # ------------------------------------------------------------------
    @testset "cldfrac = 1 (deterministic)" begin
        r1 = run_solve!(base_cld_frac)
        r2 = run_solve!(base_cld_frac)
        @test r1.lw_net == r2.lw_net
        @test r1.sw_net == r2.sw_net
    end

    # ------------------------------------------------------------------
    # 2. Uniform partial cloud (cldfrac = 0.5): one key reproduces, a new key
    #    resamples. Both halves matter: without the first the sampling is not
    #    reproducible, and without the second the key is stuck and every
    #    radiation step would draw the same clouds.
    # ------------------------------------------------------------------
    @testset "cldfrac = 0.5 (keyed sampling)" begin
        cld_frac_half = base_cld_frac .* FT(0.5)
        r1 = run_solve!(cld_frac_half; key = 7)
        r2 = run_solve!(cld_frac_half; key = 7)
        @test r1.lw_net == r2.lw_net
        @test r1.sw_net == r2.sw_net
        r3 = run_solve!(cld_frac_half; key = 8)
        diff_lw = maximum(abs.(r1.lw_net .- r3.lw_net))
        diff_sw = maximum(abs.(r1.sw_net .- r3.sw_net))
        println("cldfrac=0.5 key 7 vs 8: L∞ diff LW = $diff_lw, SW = $diff_sw")
        @test diff_lw > 0
        @test diff_sw > 0
    end

    # ------------------------------------------------------------------
    # 3. Uniform partial cloud (cldfrac = 0.5): reproducible on EVERY device.
    #    This used to be asserted on CPU only, because reproducibility came from
    #    `Random.seed!` and the GPU draws from a device RNG a host reseed cannot
    #    control. The sample is now a function of the key and the coordinates,
    #    so there is nothing device-specific left to exclude.
    # ------------------------------------------------------------------
    @testset "cldfrac = 0.5 (keyed, reproducible on any device)" begin
        cld_frac_half = base_cld_frac .* FT(0.5)
        r1 = run_solve!(cld_frac_half; key)
        r2 = run_solve!(cld_frac_half; key)
        @test r1.lw_net == r2.lw_net
        @test r1.sw_net == r2.sw_net
    end

    # ------------------------------------------------------------------
    # 4. Spatially varying cloud fractions: the same contract as 2
    # ------------------------------------------------------------------
    @testset "spatially varying cldfrac (keyed sampling)" begin
        Random.seed!(42)
        cld_frac_vary = base_cld_frac .* DA(rand(FT, nlay, ncol))
        r1 = run_solve!(cld_frac_vary; key = 11)
        r2 = run_solve!(cld_frac_vary; key = 11)
        @test r1.lw_net == r2.lw_net
        @test r1.sw_net == r2.sw_net
        r3 = run_solve!(cld_frac_vary; key = 12)
        diff_lw = maximum(abs.(r1.lw_net .- r3.lw_net))
        diff_sw = maximum(abs.(r1.sw_net .- r3.sw_net))
        println(
            "varying cldfrac key 11 vs 12: L∞ diff LW = $diff_lw, SW = $diff_sw",
        )
        @test diff_lw > 0
        @test diff_sw > 0
    end

    # ------------------------------------------------------------------
    # 5. Spatially varying cloud fractions: reproducible on every device, for
    #    the reason given in section 3.
    # ------------------------------------------------------------------
    @testset "spatially varying cldfrac (keyed, reproducible)" begin
        Random.seed!(42)
        cld_frac_vary = base_cld_frac .* DA(rand(FT, nlay, ncol))
        r1 = run_solve!(cld_frac_vary; key)
        r2 = run_solve!(cld_frac_vary; key)
        @test r1.lw_net == r2.lw_net
        @test r1.sw_net == r2.sw_net
    end

    # ------------------------------------------------------------------
    # 6. Cloud cover bounds for partial cloud
    # ------------------------------------------------------------------
    @testset "cloud cover bounds (partial cloud)" begin
        cld_frac_half = base_cld_frac .* FT(0.5)
        r = run_solve!(cld_frac_half; key)
        @test all(r.cld_cover_lw .>= 0)
        @test all(r.cld_cover_lw .<= 1)
        @test all(r.cld_cover_sw .>= 0)
        @test all(r.cld_cover_sw .<= 1)
    end
end
