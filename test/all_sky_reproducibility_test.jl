using Test
using Random
include("all_sky_utils.jl")

function run_all_sky(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{FT},
    ncol,
    use_lut;
    cldfrac = FT(1),
    seed = nothing,
) where {FT, SLVLW, SLVSW}
    device, as, lookup_lw, lookup_lw_cld, lookup_sw, lookup_sw_cld, slv_lw, slv_sw, _ =
        setup_all_sky_test(context, SLVLW, SLVSW, FT, ncol, use_lut, FT(1))
    as.cloud_state.cld_frac .*= cldfrac
    if !isnothing(seed)
        Random.seed!(seed)
    end
    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld)
    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld)
    return slv_lw, slv_sw
end

function all_sky_reproducibility_test(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{FT},
    ncol,
    use_lut,
) where {FT, SLVLW, SLVSW}
    DA = ClimaComms.array_type(context.device)
    nlev, _ = size_ds_all_sky(use_lut)
    nlay = nlev - 1
    seed = 1234

    color2 = :cyan
    printstyled(
        "Running all-sky reproducibility test on $(context.device) device with \n$FT precision, $SLVLW longwave solver and $SLVSW shortwave solver\n",
        color = color2,
    )
    printstyled("======================================================================\n", color = color2)

    printstyled("Testing with cloud fraction of 1.0\n", color = color2)
    cldfrac1_slv_lw_1, cldfrac1_slv_sw_1 = run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = FT(1))
    cldfrac1_slv_lw_2, cldfrac1_slv_sw_2 = run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = FT(1))

    diff1_lw = maximum(Array(abs.(cldfrac1_slv_lw_1.flux.flux_net .- cldfrac1_slv_lw_2.flux.flux_net)))
    diff1_sw = maximum(Array(abs.(cldfrac1_slv_sw_1.flux.flux_net .- cldfrac1_slv_sw_2.flux.flux_net)))

    println("L∞ difference in longwave flux_net = $diff1_lw")
    println("L∞ difference in shortwave flux_net = $diff1_sw")
    @test diff1_lw ≈ FT(0)
    @test diff1_sw ≈ FT(0)

    printstyled("Testing with cloud fraction of 0.5 without manually seeding random number generator\n", color = color2)
    cldfrac0p5_slv_lw_1, cldfrac0p5_slv_sw_1 = run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = FT(0.5))
    cldfrac0p5_slv_lw_2, cldfrac0p5_slv_sw_2 = run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = FT(0.5))

    diff0p5_lw = maximum(Array(abs.(cldfrac0p5_slv_lw_1.flux.flux_net .- cldfrac0p5_slv_lw_2.flux.flux_net)))
    diff0p5_sw = maximum(Array(abs.(cldfrac0p5_slv_sw_1.flux.flux_net .- cldfrac0p5_slv_sw_2.flux.flux_net)))

    println("L∞ difference in longwave flux_net = $diff0p5_lw")
    println("L∞ difference in shortwave flux_net = $diff0p5_sw")
    @test diff0p5_lw ≠ FT(0)
    @test diff0p5_sw ≠ FT(0)

    printstyled("Testing with cloud fraction of 0.5 and manually seed random number generator\n", color = color2)
    cldfrac0p5_slv_lw_3, cldfrac0p5_slv_sw_3 =
        run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = FT(0.5), seed = seed)
    cldfrac0p5_slv_lw_4, cldfrac0p5_slv_sw_4 =
        run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = FT(0.5), seed = seed)

    diff0p5_seed_lw = maximum(Array(abs.(cldfrac0p5_slv_lw_3.flux.flux_net .- cldfrac0p5_slv_lw_4.flux.flux_net)))
    diff0p5_seed_sw = maximum(Array(abs.(cldfrac0p5_slv_sw_3.flux.flux_net .- cldfrac0p5_slv_sw_4.flux.flux_net)))

    println("L∞ difference in longwave flux_net = $diff0p5_seed_lw")
    println("L∞ difference in shortwave flux_net = $diff0p5_seed_sw")
    @test diff0p5_seed_lw ≈ FT(0)
    @test diff0p5_seed_sw ≈ FT(0)

    printstyled(
        "Testing with spatially variable cloud fraction without manually seeding random number generator\n",
        color = color2,
    )

    cldfrac = DA(Random.rand(FT, nlay, ncol))
    cldfrac_random_slv_lw_1, cldfrac_random_slv_sw_1 =
        run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = cldfrac)
    cldfrac_random_slv_lw_2, cldfrac_random_slv_sw_2 =
        run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = cldfrac)

    diff_random_lw =
        maximum(Array(abs.(cldfrac_random_slv_lw_1.flux.flux_net .- cldfrac_random_slv_lw_2.flux.flux_net)))
    diff_random_sw =
        maximum(Array(abs.(cldfrac_random_slv_sw_1.flux.flux_net .- cldfrac_random_slv_sw_2.flux.flux_net)))

    println("L∞ difference in longwave flux_net = $diff_random_lw")
    println("L∞ difference in shortwave flux_net = $diff_random_sw")
    @test diff_random_lw ≠ FT(0)
    @test diff_random_sw ≠ FT(0)


    printstyled(
        "Testing with spatially variable cloud fraction with manually seeding random number generator\n",
        color = color2,
    )
    cldfrac_random_slv_lw_3, cldfrac_random_slv_sw_3 =
        run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = cldfrac, seed = seed)
    cldfrac_random_slv_lw_4, cldfrac_random_slv_sw_4 =
        run_all_sky(context, SLVLW, SLVSW, FT, ncol, use_lut, cldfrac = cldfrac, seed = seed)

    diff_random_seed_lw =
        maximum(Array(abs.(cldfrac_random_slv_lw_3.flux.flux_net .- cldfrac_random_slv_lw_4.flux.flux_net)))
    diff_random_seed_sw =
        maximum(Array(abs.(cldfrac_random_slv_sw_3.flux.flux_net .- cldfrac_random_slv_sw_4.flux.flux_net)))

    println("L∞ difference in longwave flux_net = $diff_random_seed_lw")
    println("L∞ difference in shortwave flux_net = $diff_random_seed_sw")
    @test diff_random_seed_lw ≈ FT(0)
    @test diff_random_seed_sw ≈ FT(0)
end
