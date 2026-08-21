#=
GPU benchmark ratchet for the high-resolution solver kernels.

Measures the minimum wall time of `solve_lw!`/`solve_sw!` at high resolution
for the three production configurations (clear sky, all sky, all sky with
aerosols; two-stream longwave and shortwave) in both precisions, and compares
each against a recorded per-GPU baseline:

    perf/benchmark_baselines/<device_key>.txt

A measurement slower than baseline × (1 + tolerance) fails its `@test`. The
tolerance (default 20%, chosen for shared-cluster timing noise on the
Buildkite P100s) can be overridden with the `RRTMGP_BENCH_TOL_PCT` environment
variable.

Baselines:
- If no baseline file exists for this GPU, the script records one and passes;
  commit the file (CI uploads it as a Buildkite artifact) to arm the ratchet.
- To re-baseline intentionally (e.g. after an accepted performance change),
  run with `RRTMGP_BENCH_WRITE_BASELINE=1` and commit the updated file.

Usage (needs a CUDA device):

    CLIMACOMMS_DEVICE=CUDA julia --project=test perf/benchmark_ratchet.jl

This lane is advisory (soft-fail in Buildkite): it flags regressions loudly
without blocking merges, per the CliMA convention that `perf/` is not a gate.
=#

using Test
using Statistics
using Printf
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

const context = ClimaComms.context()
const device = ClimaComms.device(context)

if !(device isa ClimaComms.CUDADevice)
    @info "benchmark_ratchet.jl requires a CUDA device; nothing to do." device
    exit(0)
end

using CUDA

# The target column count used by the benchmark scripts: a cubed sphere of
# 30² horizontal elements × 6 panels × 4² quadrature points, i.e. 86,400
# columns at roughly 80-km spacing, across 64 levels / 63 layers.
highres_ncol(; helems = 30, nq = 4) = Int(helems * helems * 6 * nq * nq)
const highres_nlay = 63

# Pull in benchmark_clear_sky / benchmark_all_sky /
# benchmark_all_sky_with_aerosols; their sweep loops are guarded behind
# `PROGRAM_FILE == @__FILE__`, so the includes only define functions.
const test_dir = joinpath(dirname(@__DIR__), "test")
include(joinpath(test_dir, "clear_sky_highres_gpu_benchmark.jl"))
include(joinpath(test_dir, "cloudy_sky_highres_gpu_benchmark.jl"))
include(joinpath(test_dir, "all_sky_with_aerosols_highres_gpu_benchmark.jl"))

# One measured case: production two-stream solvers across identical uniform
# high-resolution grids. Each `run` returns (trial_lw, trial_sw) from
# BenchmarkTools.
cases = [
    (
        key = "clear_sky",
        run = FT -> benchmark_clear_sky(
            context,
            TwoStreamLWRTE,
            TwoStreamSWRTE,
            VmrGM,
            FT;
            ncol = highres_ncol(),
            nlay = highres_nlay,
        ),
    ),
    (
        key = "all_sky",
        run = FT -> benchmark_all_sky(
            context,
            TwoStreamLWRTE,
            TwoStreamSWRTE,
            FT;
            ncol = highres_ncol(),
            nlay = highres_nlay,
            cldfrac = FT(1),
        ),
    ),
    (
        key = "all_sky_with_aerosols",
        run = FT -> benchmark_all_sky_with_aerosols(
            context,
            TwoStreamLWRTE,
            TwoStreamSWRTE,
            FT;
            ncol = highres_ncol(),
            nlay = highres_nlay,
        ),
    ),
]

device_key() =
    replace(lowercase(CUDA.name(CUDA.device())), r"[^a-z0-9]+" => "_")

baseline_path(key) =
    joinpath(@__DIR__, "benchmark_baselines", string(key, ".txt"))

# Baseline file format: comment lines start with '#'; data lines are two
# whitespace-separated tokens, "<case>_<lw|sw>_<FT> <minimum_ns>". The key is a
# single token (no embedded spaces) so the round-trip is unambiguous.
function read_baseline(path)
    baseline = Dict{String, Float64}()
    isfile(path) || return nothing
    for line in eachline(path)
        line = strip(line)
        (isempty(line) || startswith(line, "#")) && continue
        name, ns = split(line)
        baseline[name] = parse(Float64, ns)
    end
    return baseline
end

function write_baseline(path, results; kind = "baseline")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# RRTMGP high-resolution GPU benchmark $kind (minimum ns)")
        println(io, "# device: $(CUDA.name(CUDA.device()))")
        println(io, "# julia: $(VERSION)")
        println(io, "# recorded: ", get(ENV, "BUILDKITE_BUILD_URL", "local run"))
        println(io, "# regenerate with RRTMGP_BENCH_WRITE_BASELINE=1")
        for (name, ns) in sort(collect(results); by = first)
            println(io, "$name $ns")
        end
    end
    return path
end

const tol = parse(Float64, get(ENV, "RRTMGP_BENCH_TOL_PCT", "20")) / 100
const write_mode = get(ENV, "RRTMGP_BENCH_WRITE_BASELINE", "0") == "1"
const path = baseline_path(device_key())
const baseline = read_baseline(path)

results = Dict{String, Float64}()
for FT in (Float32, Float64), case in cases
    @info "benchmark_ratchet: running $(case.key) ($FT) at high resolution"
    trial_lw, trial_sw = case.run(FT)
    # `minimum` is the recommended timing statistic: it is the most resilient to
    # GC/scheduler noise on a shared cluster (CliMA DeveloperGuides,
    # performance/allocation_debugging.md §5), which also makes the ratchet
    # threshold reliable.
    results[string(case.key, "_lw_", FT)] = minimum(trial_lw).time # ns
    results[string(case.key, "_sw_", FT)] = minimum(trial_sw).time # ns
end

# Always persist the measured minima as an artifact, whatever the mode.
const measured_path =
    joinpath(@__DIR__, "benchmark_baselines", "measured_$(device_key()).txt")
write_baseline(measured_path, results; kind = "measurement")

println("\ndevice: $(CUDA.name(CUDA.device())) (key: $(device_key()))")
println("tolerance: +$(round(Int, 100 * tol))% over baseline\n")
@printf(
    "%-40s %14s %14s %9s\n",
    "case",
    "minimum [ms]",
    "baseline [ms]",
    "ratio"
)
for (name, ns) in sort(collect(results); by = first)
    base = isnothing(baseline) ? missing : get(baseline, name, missing)
    @printf(
        "%-40s %14.3f %14s %9s\n",
        name,
        ns / 1e6,
        ismissing(base) ? "—" : @sprintf("%.3f", base / 1e6),
        ismissing(base) ? "—" : @sprintf("%.2fx", ns / base),
    )
end

if write_mode || isnothing(baseline)
    write_baseline(path, results)
    action = write_mode ? "rewrote" : "no baseline found for this GPU; recorded"
    @info "benchmark_ratchet: $action $path — commit it to arm the ratchet."
    println("\n--- baseline file contents ---")
    println(read(path, String))
else
    @testset "GPU benchmark ratchet ($(device_key()))" begin
        for (name, ns) in sort(collect(results); by = first)
            base = get(baseline, name, missing)
            if ismissing(base)
                @warn "no baseline entry for $name; rerun with " *
                      "RRTMGP_BENCH_WRITE_BASELINE=1 to record it"
                continue
            end
            @test ns <= base * (1 + tol)
        end
    end
end
