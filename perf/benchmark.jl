#=
julia --project=examples
julia --project=examples perf/benchmark.jl
```
include(joinpath("perf", "benchmark.jl"))
```
=#
# ARGS[1]
const FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64
using BenchmarkTools
using Suppressor
root_dir = joinpath(dirname(@__DIR__))
import ClimaComms
import Logging

@info "------------------------------------------------- Benchmark: gray_atm"
@suppress_out begin
    include(joinpath(root_dir, "test", "gray_atm_utils.jl"))
    gray_atmos_lw_equil(ClimaComms.context(), OneScalar, FT; exfiltrate = true)
end
(; slv, max_threads) = Infiltrator.exfiltrated
@info "gray_atm lw"
solve_lw!(slv, max_threads) # compile first
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv, $max_threads)
else
    @benchmark solve_lw!($slv, $max_threads)
end
show(stdout, MIME("text/plain"), trial)
println()

gray_atmos_sw_test(ClimaComms.context(), OneScalar, FT, 1; exfiltrate = true)
(; slv, max_threads) = Infiltrator.exfiltrated
solve_sw!(slv, max_threads) # compile first
@info "gray_atm sw"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_sw!($slv, $max_threads)
else
    @benchmark solve_sw!($slv, $max_threads)
end
show(stdout, MIME("text/plain"), trial)
println()
@info "------------------------------------------------- Benchmark: clear_sky"
# @suppress_out begin
include(joinpath(root_dir, "test", "clear_sky_utils.jl"))
context = ClimaComms.context()
clear_sky(ClimaComms.context(), TwoStream, SourceLW2Str, VmrGM, FT; exfiltrate = true)
# end
(; slv, max_threads, lookup_sw, lookup_lw) = Infiltrator.exfiltrated

@info "clear_sky lw"
solve_lw!(slv, max_threads, lookup_lw) # compile first
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv, $max_threads, $lookup_lw)
else
    @benchmark solve_lw!($slv, $max_threads, $lookup_lw)
end
show(stdout, MIME("text/plain"), trial)
println()

@info "clear_sky sw"
solve_sw!(slv, max_threads, lookup_sw) # compile first
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_sw!($slv, $max_threads, $lookup_sw)
else
    @benchmark solve_sw!($slv, $max_threads, $lookup_sw)
end
show(stdout, MIME("text/plain"), trial)
println()

@info "------------------------------------------------- Benchmark: all_sky"
# @suppress_out begin
ncols_all_sky_benchmark = 60000
include(joinpath(root_dir, "test", "all_sky_utils.jl"))
all_sky(
    ClimaComms.context(),
    TwoStream,
    FT;
    ncol = ncols_all_sky_benchmark,
    use_lut = true,
    cldfrac = FT(1),
    exfiltrate = true,
)
# end

(; slv, max_threads, lookup_sw, lookup_sw_cld, lookup_lw, lookup_lw_cld) = Infiltrator.exfiltrated

solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld) # compile first
solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld) # compile first

@info "all_sky, lw, use_lut=true"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv, $max_threads, $lookup_lw, $lookup_lw_cld)
else
    @benchmark solve_lw!($slv, $max_threads, $lookup_lw, $lookup_lw_cld)
end
show(stdout, MIME("text/plain"), trial)
println()
@info "all_sky, sw, use_lut=true"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_sw!($slv, $max_threads, $lookup_sw, $lookup_sw_cld)
else
    @benchmark solve_sw!($slv, $max_threads, $lookup_sw, $lookup_sw_cld)
end
show(stdout, MIME("text/plain"), trial)
println()

# @suppress_out begin
all_sky(
    ClimaComms.context(),
    TwoStream,
    FT;
    ncol = ncols_all_sky_benchmark,
    use_lut = false,
    cldfrac = FT(1),
    exfiltrate = true,
)
# end

(; slv, max_threads, lookup_sw, lookup_sw_cld, lookup_lw, lookup_lw_cld) = Infiltrator.exfiltrated

solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld) # compile first
solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld) # compile first

@info "all_sky, lw, use_lut=false"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv, $max_threads, $lookup_lw, $lookup_lw_cld)
else
    @benchmark solve_lw!($slv, $max_threads, $lookup_lw, $lookup_lw_cld)
end
show(stdout, MIME("text/plain"), trial)
println()
@info "all_sky, sw, use_lut=false"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_sw!($slv, $max_threads, $lookup_sw, $lookup_sw_cld)
else
    @benchmark solve_sw!($slv, $max_threads, $lookup_sw, $lookup_sw_cld)
end
show(stdout, MIME("text/plain"), trial)
println()
