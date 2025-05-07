#=
julia --project=examples
julia --project=examples perf/flame.jl gray_atm.jl
julia --project=examples perf/flame.jl clear_sky.jl
julia --project=examples perf/flame.jl all_sky.jl
```
include(joinpath("perf", "flame.jl"))
```
=#
# ARGS[1]

filename = isempty(ARGS) ? "all_sky.jl" : ARGS[1]
root_dir = joinpath(dirname(@__DIR__))
output_dir = joinpath(root_dir, "flame_graphs", splitext(filename)[1])
file = joinpath(root_dir, "test", filename)
@info "Collecting flame for $file"
@info "Output directory $output_dir"

@info "collect profile"
import Profile
import ProfileCanvas
Profile.clear()
prof = Profile.@profile begin
    include(file)
end
mkpath(output_dir)

results = Profile.fetch()
Profile.clear()
ProfileCanvas.html_file(joinpath(output_dir, "flame.html"), results)

# Nearly all allocations are for inference, which
# is difficult to map to src code. Need to split
# test directory into utils so that we target
# specific function calls.

# @info "collecting allocations"
# Profile.Allocs.clear()
# Profile.Allocs.@profile sample_rate = 0.01 include(file)
# results = Profile.Allocs.fetch()
# Profile.Allocs.clear()
# profile = ProfileCanvas.view_allocs(results)
# ProfileCanvas.html_file(joinpath(output_dir, "allocs.html"), profile)
