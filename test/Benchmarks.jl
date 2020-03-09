using Test
using RRTMGP
using RRTMGP.OpticalProps
using Profile
using Pkg
Pkg.add("BenchmarkTools")
using BenchmarkTools

include("rfmip_clear_sky_sw.jl")
include("rfmip_clear_sky_lw.jl")
include("allsky.jl")
include("DataSetFiles.jl")

@testset "Benchmark suite" begin
    datafolder = RRTMGP.data_folder_rrtmgp()

    ds_clearsky_sw = dataset_dict(data_files_dict(datafolder, "sw"))
    ds_clearsky_lw = dataset_dict(data_files_dict(datafolder, "lw"))
    ds_allsky_sw = dataset_dict(data_files_dict(
        datafolder,
        "sw";
        allsky = true,
    ))
    ds_allsky_lw = dataset_dict(data_files_dict(
        datafolder,
        "lw";
        allsky = true,
    ))

    refresh_pr = true
    refresh_params = true
    refresh_master = true
    params_file = joinpath("BenchmarkResults", "params.json")
    master_file = joinpath("BenchmarkResults", "master.json")
    pr_file = joinpath("BenchmarkResults", "pr.json")
    refresh_params && rm(params_file, force = true)
    refresh_master && rm(master_file, force = true)
    refresh_pr && rm(pr_file, force = true)

    suite = BenchmarkGroup()
    !refresh_params && loadparams!(
        suite,
        BenchmarkTools.load(params_file)[1],
        :evals,
        :samples,
    )
    !refresh_master && (master = BenchmarkTools.load(master_file)[1])

    suite["allsky"] = BenchmarkGroup(["clouds"])
    for use_luts in (true, false)
        suite["allsky"]["sw", string(use_luts)] = @benchmarkable all_sky(
            $(ds_allsky_sw);
            use_luts = $(use_luts),
            λ_string = $("sw"),
        )
    end

    for use_luts in (true, false)
        suite["allsky"]["lw", string(use_luts)] = @benchmarkable all_sky(
            $(ds_allsky_lw);
            use_luts = $(use_luts),
            λ_string = $("lw"),
        )
    end
    suite["clearsky"] = BenchmarkGroup(["clear"])
    suite["clearsky"]["lw_1scl"] = @benchmarkable rfmip_clear_sky_lw(
        $(ds_clearsky_lw),
        OneScalar,
    )
    suite["clearsky"]["lw_2str"] = @benchmarkable rfmip_clear_sky_lw(
        $(ds_clearsky_lw),
        TwoStream,
    )
    suite["clearsky"]["sw_1scl"] = @benchmarkable rfmip_clear_sky_sw(
        $(ds_clearsky_sw),
        OneScalar,
    )
    suite["clearsky"]["sw_2str"] = @benchmarkable rfmip_clear_sky_sw(
        $(ds_clearsky_sw),
        TwoStream,
    )

    refresh_params && tune!(suite)

    pr = refresh_pr ? run(suite, verbose = true) :
         BenchmarkTools.load(pr_file)[1]
    refresh_master && (master = deepcopy(pr))

    refresh_params && BenchmarkTools.save(params_file, params(suite))
    refresh_master && BenchmarkTools.save(master_file, master)
    refresh_pr && BenchmarkTools.save(pr_file, pr)

    regs = regressions(judge(minimum(pr), minimum(master))) # a BenchmarkGroup containing the regressions
    @show leaves(regs) # an array of (ID, `TrialJudgement`) pairs

    close.(values(ds_clearsky_sw))
    close.(values(ds_clearsky_lw))
    close.(values(ds_allsky_sw))
    close.(values(ds_allsky_lw))

  # Profile.clear()
  # Profile.@profile run_driver(datafolder, OneScalar)
  # Profile.@profile run_driver(datafolder, TwoStream)
  # Profile.print()
end
