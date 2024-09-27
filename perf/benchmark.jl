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

context = ClimaComms.context()
@info "------------------------------------------------- Benchmark: gray_atm"
include(joinpath(root_dir, "test", "gray_atm_utils.jl"))
gray_as, slv_lw, _ = setup_gray_atmos_lw_equil_test(context, NoScatLWRTE, FT)

@info "gray_atm lw"
solve_lw!(slv_lw, gray_as) # compile first
device = ClimaComms.device(context)
trial = if device isa ClimaComms.CUDADevice
    using CUDA
    @benchmark CUDA.@sync solve_lw!($slv_lw, $gray_as)
else
    @benchmark solve_lw!($slv_lw, $gray_as)
end
show(stdout, MIME("text/plain"), trial)
println()

as, slv_sw, _ = setup_gray_atmos_sw_test(context, NoScatSWRTE, FT, 1)

solve_sw!(slv_sw, as) # compile first
@info "gray_atm sw"
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

toler_lw_noscat = Dict(Float64 => Float64(1e-4), Float32 => Float32(0.04))
toler_lw_2stream = Dict(Float64 => Float64(4.5), Float32 => Float32(4.5))
toler_sw = Dict(Float64 => Float64(1e-3), Float32 => Float32(0.04))

device, as, lookup_lw, lookup_sw, slv_lw, slv_sw =
    setup_clear_sky_test(context, TwoStreamLWRTE, TwoStreamSWRTE, VmrGM, FT, toler_lw_2stream, toler_sw)
# end

@info "clear_sky lw"
solve_lw!(slv_lw, as, lookup_lw) # compile first
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

toler_lw_noscat = Dict(Float64 => Float64(1e-5), Float32 => Float32(0.05))
toler_lw_2stream = Dict(Float64 => Float64(5), Float32 => Float32(5))
toler_sw = Dict(Float64 => Float64(1e-5), Float32 => Float32(0.06))
use_lut = true
cldfrac = FT(1)
ncol = 128
device, as, lookup_lw, lookup_lw_cld, lookup_sw, lookup_sw_cld, slv_lw, slv_sw, _ =
    setup_all_sky_test(ClimaComms.context(), TwoStreamLWRTE, TwoStreamSWRTE, FT, ncol, use_lut, cldfrac)

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
