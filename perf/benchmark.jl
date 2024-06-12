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
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import Logging

@info "------------------------------------------------- Benchmark: gray_atm"
@suppress_out begin
    include(joinpath(root_dir, "test", "gray_atm_utils.jl"))
    gray_atmos_lw_equil(ClimaComms.context(), OneScalar, NoScatLWRTE, FT; exfiltrate = true)
end
(; slv_lw, gray_as) = Infiltrator.exfiltrated
@info "gray_atm lw"
solve_lw!(slv_lw, gray_as) # compile first
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv_lw, $gray_as)
else
    @benchmark solve_lw!($slv_lw, $gray_as)
end
show(stdout, MIME("text/plain"), trial)
println()

gray_atmos_sw_test(ClimaComms.context(), OneScalar, NoScatSWRTE, FT, 1; exfiltrate = true)
(; slv_sw, as) = Infiltrator.exfiltrated
solve_sw!(slv_sw, as) # compile first
@info "gray_atm sw"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_sw!($slv_sw, $as)
else
    @benchmark solve_sw!($slv_sw, $as)
end
show(stdout, MIME("text/plain"), trial)
println()
@info "------------------------------------------------- Benchmark: clear_sky"
# @suppress_out begin
include(joinpath(root_dir, "test", "clear_sky_utils.jl"))
context = ClimaComms.context()
clear_sky(ClimaComms.context(), TwoStream, TwoStream, TwoStreamLWRTE, TwoStreamSWRTE, VmrGM, FT; exfiltrate = true)
# end
(; slv_lw, slv_sw, as, lookup_sw, lookup_lw) = Infiltrator.exfiltrated

@info "clear_sky lw"
solve_lw!(slv_lw, as, lookup_lw) # compile first
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv_lw, $as, $lookup_lw)
else
    @benchmark solve_lw!($slv_lw, $as, $lookup_lw)
end
show(stdout, MIME("text/plain"), trial)
println()

@info "clear_sky sw"
solve_sw!(slv_sw, as, lookup_sw) # compile first
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_sw!($slv_sw, $as, $lookup_sw)
else
    @benchmark solve_sw!($slv_sw, $as, $lookup_sw)
end
show(stdout, MIME("text/plain"), trial)
println()

@info "------------------------------------------------- Benchmark: all_sky"
# @suppress_out begin
include(joinpath(root_dir, "test", "all_sky_utils.jl"))
all_sky(
    ClimaComms.context(),
    TwoStream,
    TwoStream,
    TwoStreamLWRTE,
    TwoStreamSWRTE,
    FT;
    use_lut = true,
    cldfrac = FT(1),
    exfiltrate = true,
)
# end

(; slv_lw, slv_sw, as, lookup_sw, lookup_sw_cld, lookup_lw, lookup_lw_cld) = Infiltrator.exfiltrated

solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld) # compile first
solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld) # compile first

@info "all_sky, lw, use_lut=true"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv_lw, $as, $lookup_lw, $lookup_lw_cld)
else
    @benchmark solve_lw!($slv_lw, $as, $lookup_lw, $lookup_lw_cld)
end
show(stdout, MIME("text/plain"), trial)
println()
@info "all_sky, sw, use_lut=true"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_sw!($slv_sw, $as, $lookup_sw, $lookup_sw_cld)
else
    @benchmark solve_sw!($slv_sw, $as, $lookup_sw, $lookup_sw_cld)
end
show(stdout, MIME("text/plain"), trial)
println()

#=
# @suppress_out begin
all_sky(
    ClimaComms.context(),
    TwoStream,
    TwoStream,
    TwoStreamLWRTE,
    TwoStreamSWRTE,
    FT;
    use_lut = false,
    cldfrac = FT(1),
    exfiltrate = true,
)
# end

(; slv_lw, slv_sw, as, lookup_sw, lookup_sw_cld, lookup_lw, lookup_lw_cld) = Infiltrator.exfiltrated

solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld) # compile first
solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld) # compile first

@info "all_sky, lw, use_lut=false"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv_lw, $as, $lookup_lw, $lookup_lw_cld)
else
    @benchmark solve_lw!($slv_lw, $as, $lookup_lw, $lookup_lw_cld)
end
show(stdout, MIME("text/plain"), trial)
println()
@info "all_sky, sw, use_lut=false"
device = ClimaComms.device(ClimaComms.context())
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_sw!($slv_sw, $as, $lookup_sw, $lookup_sw_cld)
else
    @benchmark solve_sw!($slv_sw, $as, $lookup_sw, $lookup_sw_cld)
end
show(stdout, MIME("text/plain"), trial)
println()
=#
