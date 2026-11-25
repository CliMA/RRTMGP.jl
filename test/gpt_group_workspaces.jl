using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
import RRTMGP.RTE:
    LWGptGroups, SWGptGroups, default_ngpt_groups, gpt_group_bytes

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

# The byte bound is a device property, and CUDA reports memory unsigned, so the
# count can come back a `UInt64` on a GPU and nowhere else -- where it no longer
# matches the `ngpt_groups::Int` keyword it is the default for. This runs on
# whatever device the tests are on, which is the only way to see that.
@testset "the default group count is an Int" begin
    context = ClimaComms.context()
    grid_params = RRTMGP.RRTMGPGridParams(
        Float32;
        context,
        domain_nlay = 80,
        ncol = 2048,
        isothermal_boundary_layer = false,
    )
    @test default_ngpt_groups(grid_params) isa Int
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
        # With and without the per-band aerosol stores, at band counts that
        # differ between the two bands as the real ones do.
        for (nbnd_lw, nbnd_sw) in ((0, 0), (16, 14))
            built =
                workspace_bytes(
                    LWGptGroups(
                        grid_params,
                        ngroups;
                        params,
                        naero_bands = nbnd_lw,
                    ),
                ) + workspace_bytes(
                    SWGptGroups(grid_params, ngroups; naero_bands = nbnd_sw),
                )
            counted =
                gpt_group_bytes(FT, ncol, nlay, ngroups, nbnd_lw + nbnd_sw)
            # The optics and the mask are over layers and are charged a level
            # each, which is the whole of the slack: six optics planes and two
            # mask bytes per wide column. Nothing else may differ, in either
            # direction -- the aerosol stores are over layers and counted so.
            slack = ncol * ngroups * (6 * sizeof(FT) + 2 * sizeof(Bool))
            @test counted - built == slack
        end
    end
end

# The aerosol store is paid once whatever the group count, so it comes off the
# byte bound before the groups are priced: more bands can only leave fewer
# groups affordable, and never fewer than the one that needs no workspaces.
@testset "the aerosol store counts against the default group count" begin
    context = ClimaComms.context()
    grid_params = RRTMGP.RRTMGPGridParams(
        Float32;
        context,
        domain_nlay = 80,
        ncol = 2048,
        isothermal_boundary_layer = false,
    )
    counts = [
        default_ngpt_groups(grid_params; naero_bands) for
        naero_bands in (0, 30, 10^4, 10^7)
    ]
    @test issorted(counts; rev = true)
    @test first(counts) == default_ngpt_groups(grid_params)
    @test last(counts) == 1
    @test gpt_group_bytes(Float32, 2048, 80, 0, 30) == 3 * 30 * 2048 * 80 * 4
end
