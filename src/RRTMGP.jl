module RRTMGP


import ClimaComms

# TODO: should this be moved to ClimaComms?
#=
    @threaded device for ... end

A threading macro that uses Julia native threading if the
device is a `CPUMultiThreaded` type, otherwise return
the original expression without `Threads.@threads`. This is
done to avoid overhead from `Threads.@threads`, and the device
is used (instead of checking `Threads.nthreads() == 1`) so
that this is statically inferred.

## References
 - https://discourse.julialang.org/t/threads-threads-with-one-thread-how-to-remove-the-overhead/58435
 - https://discourse.julialang.org/t/overhead-of-threads-threads/53964

=#
macro threaded(device, expr)
    return esc(quote
        let
            if $device isa ClimaComms.CPUMultiThreaded
                Threads.@threads $(expr)
            else
                @assert $device isa ClimaComms.CPUSingleThreaded
                $(expr)
            end
        end
    end)
end

include("Parameters.jl")
import .Parameters as RP

include(joinpath("optics", "Vmrs.jl"))
include(joinpath("optics", "LookUpTables.jl"))
include(joinpath("optics", "AngularDiscretizations.jl"))
include(joinpath("optics", "AtmosphericStates.jl"))
include(joinpath("optics", "Sources.jl"))
include(joinpath("optics", "Optics.jl"))
include(joinpath("optics", "Fluxes.jl"))
include(joinpath("optics", "GrayUtils.jl"))
include(joinpath("optics", "BCs.jl"))
include(joinpath("optics", "RTE.jl"))

include(joinpath("rte", "RTESolver.jl"))

end # module
