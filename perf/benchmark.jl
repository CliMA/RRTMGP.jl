#=
julia --project=examples
julia --project=examples perf/benchmark.jl
```
include(joinpath("perf", "benchmark.jl"))
```
=#
# ARGS[1]

using BenchmarkTools
using Suppressor
root_dir = joinpath(dirname(@__DIR__))
import ClimaComms
import Logging

@info "------------------------------------------------- Benchmark: gray_atm"
@suppress_out begin
    include(joinpath(root_dir, "test", "gray_atm_utils.jl"))
    gray_atmos_lw_equil(ClimaComms.context(), OneScalar, Float64; exfiltrate = true)
end
(; slv, max_threads) = Infiltrator.exfiltrated
@info "gray_atm lw"
solve_lw!(slv, max_threads) # compile first
trial = @benchmark solve_lw!(slv, max_threads)
show(stdout, MIME("text/plain"), trial)
println()

gray_atmos_sw_test(ClimaComms.context(), OneScalar, Float64, 1; exfiltrate = true)
(; slv, max_threads) = Infiltrator.exfiltrated
solve_sw!(slv, max_threads) # compile first
@info "gray_atm sw"
trial = @benchmark solve_sw!(slv, max_threads)
show(stdout, MIME("text/plain"), trial)
println()
@info "------------------------------------------------- Benchmark: clear_sky"
@suppress_out begin
    include(joinpath(root_dir, "test", "clear_sky_utils.jl"))
    context = ClimaComms.context()
    clear_sky(ClimaComms.context(), TwoStream, SourceLW2Str, VmrGM, Float64; exfiltrate = true)
end
(; slv, max_threads, lookup_sw, lookup_lw) = Infiltrator.exfiltrated

@info "clear_sky lw"
solve_lw!(slv, max_threads, lookup_lw) # compile first
trial = @benchmark solve_lw!(slv, max_threads, lookup_lw)
show(stdout, MIME("text/plain"), trial)
println()

@info "clear_sky sw"
solve_sw!(slv, max_threads, lookup_sw) # compile first
trial = @benchmark solve_sw!(slv, max_threads, lookup_sw)
show(stdout, MIME("text/plain"), trial)
println()

@info "------------------------------------------------- Benchmark: all_sky"
@suppress_out begin
    include(joinpath(root_dir, "test", "all_sky_utils.jl"))
    all_sky(ClimaComms.context(), TwoStream, Float64; use_lut = true, cldfrac = Float64(1), exfiltrate = true)
end

(; slv, max_threads, lookup_sw, lookup_sw_cld, lookup_lw, lookup_lw_cld) = Infiltrator.exfiltrated

solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld) # compile first
solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld) # compile first

@info "all_sky, lw, use_lut=true"
trial = @benchmark solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld)
show(stdout, MIME("text/plain"), trial)
println()
@info "all_sky, sw, use_lut=true"
trial = @benchmark solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld)
show(stdout, MIME("text/plain"), trial)
println()

@suppress_out begin
    all_sky(ClimaComms.context(), TwoStream, Float64; use_lut = false, cldfrac = Float64(1), exfiltrate = true)
end

(; slv, max_threads, lookup_sw, lookup_sw_cld, lookup_lw, lookup_lw_cld) = Infiltrator.exfiltrated

solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld) # compile first
solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld) # compile first

@info "all_sky, lw, use_lut=false"
trial = @benchmark solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld)
show(stdout, MIME("text/plain"), trial)
println()
@info "all_sky, sw, use_lut=false"
trial = @benchmark solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld)
show(stdout, MIME("text/plain"), trial)
println()