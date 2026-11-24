using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
import RRTMGP.RTE: LWGptGroups, SWGptGroups, gpt_group_bytes

# `default_ngpt_groups` bounds the group count by what the widened workspaces
# cost, and `gpt_group_bytes` counts that cost from the planes the constructors
# allocate. A field added to any of the structs those constructors build would
# make the count silently low, so it is checked against the arrays themselves
# rather than against a figure measured once.

# Bytes of every distinct array a workspace holds. The views share parents, so
# each is resolved to the array that owns the memory and counted once.
function workspace_bytes(x, seen = Base.IdSet{Any}())
    if x isa AbstractArray
        root = x
        while root !== parent(root)
            root = parent(root)
        end
        root in seen && return 0
        push!(seen, root)
        return sizeof(root)
    end
    isbits(x) && return 0
    total = 0
    for name in fieldnames(typeof(x))
        isdefined(x, name) &&
            (total += workspace_bytes(getfield(x, name), seen))
    end
    return total
end

# The workspaces are built on the CPU whatever the test device is: `sizeof` of a
# device array is its handle, not the memory the count is about.
@testset "grouped workspace bytes ($FT)" for FT in (Float32, Float64)
    context = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    params = RRTMGP.default_parameters(FT)
    for ngroups in (1, 2, 5), (domain_nlay, ncol) in ((10, 16), (31, 48))
        grid_params = RRTMGP.RRTMGPGridParams(
            FT;
            context,
            domain_nlay,
            ncol,
            isothermal_boundary_layer = false,
        )
        (; nlay) = grid_params
        built =
            workspace_bytes(LWGptGroups(grid_params, ngroups; params)) +
            workspace_bytes(SWGptGroups(grid_params, ngroups))
        counted = gpt_group_bytes(FT, ncol, nlay, ngroups)
        # The optics and the mask are over layers and are charged a level each,
        # which is the whole of the slack: six optics planes and two mask bytes
        # per wide column. Nothing else may differ, in either direction.
        slack = ncol * ngroups * (6 * sizeof(FT) + 2 * sizeof(Bool))
        @test counted - built == slack
    end
end
