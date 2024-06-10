module RTESolver

using StaticArrays
import ClimaComms
using Adapt

using ..AngularDiscretizations
using ..Vmrs
using ..AtmosphericStates
using ..Fluxes
using ..BCs
using ..Sources
using ..RTE
using ..Optics
using ..LookUpTables

export solve_lw!, solve_sw!

include("longwave1scalar.jl")
include("longwave2stream.jl")
include("shortwave1scalar.jl")
include("shortwave2stream.jl")

"""
    solve_lw!((; context, flux, src, bcs, op)::NoScatLWRTE, as::GrayAtmosphericState, max_threads)

Non-scattering RTE solver for the longwave problem, using gray optics.
"""
solve_lw!((; context, flux, src, bcs, op)::NoScatLWRTE, as::GrayAtmosphericState, max_threads) =
    rte_lw_noscat_solve!(context.device, flux, src, bcs, op, max_threads, as)

"""
    solve_lw!((; context, flux, src, bcs, op)::TwoStreamLWRTE, as::GrayAtmosphericState, max_threads)

`Two Stream` RTE solver for the longwave problem, using gray optics.
"""
solve_lw!((; context, flux, src, bcs, op)::TwoStreamLWRTE, as::GrayAtmosphericState, max_threads) =
    rte_lw_2stream_solve!(context.device, flux, src, bcs, op, max_threads, as)

"""
    solve_lw!(
        (; context, fluxb, flux, src, bcs, op)::NoScatLWRTE,
        as::AtmosphericState,
        max_threads,
        lookup_lw::LookUpLW,
        lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing},
    )

Non-scattering RTE solver for the longwave problem, using RRTMGP optics.
"""
solve_lw!(
    (; context, fluxb, flux, src, bcs, op)::NoScatLWRTE,
    as::AtmosphericState,
    max_threads,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing},
) = rte_lw_noscat_solve!(context.device, fluxb, flux, src, bcs, op, max_threads, as, lookup_lw, lookup_lw_cld)

"""
    solve_lw!(
        (; context, fluxb, flux, src, bcs, op)::TwoStreamLWRTE,
        as::AtmosphericState,
        max_threads,
        lookup_lw::LookUpLW,
        lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing},
    )

`Two Stream` RTE solver for the longwave problem, using RRTMGP optics.
"""
solve_lw!(
    (; context, fluxb, flux, src, bcs, op)::TwoStreamLWRTE,
    as::AtmosphericState,
    max_threads,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing},
) = rte_lw_2stream_solve!(context.device, fluxb, flux, src, bcs, op, max_threads, as, lookup_lw, lookup_lw_cld)


"""
    solve_sw!((; context, flux, bcs, op)::NoScatSWRTE, as::GrayAtmosphericState, max_threads)

Non-scattering RTE solver for the shortwave problem, using gray optics.
"""
solve_sw!((; context, flux, bcs, op)::NoScatSWRTE, as::GrayAtmosphericState, max_threads) =
    rte_sw_noscat_solve!(context.device, flux, op, bcs, max_threads, as) # non-scattering solver, gray optics

"""
    solve_sw!((; context, flux, src, bcs, op)::TwoStreamSWRTE, as::GrayAtmosphericState, max_threads)

`Two Stream` RTE solver for the shortwave problem, using gray optics.
"""
solve_sw!((; context, flux, src, bcs, op)::TwoStreamSWRTE, as::GrayAtmosphericState, max_threads) =
    rte_sw_2stream_solve!(context.device, flux, op, bcs, src, max_threads, as)

"""
    solve_sw!(
        (; context, fluxb, flux, bcs, op)::NoScatSWRTE,
        as::AtmosphericState,
        max_threads,
        lookup_sw::LookUpSW,
        lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing},
    )

Non-scattering RTE solver for the shortwave problem, using RRTMGP optics.
"""
solve_sw!(
    (; context, fluxb, flux, bcs, op)::NoScatSWRTE,
    as::AtmosphericState,
    max_threads,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing},
) = rte_sw_noscat_solve!(context.device, fluxb, flux, op, bcs, max_threads, as, lookup_sw)

"""
    solve_sw!(
        (; context, fluxb, flux, src, bcs, op)::TwoStreamSWRTE,
        as::AtmosphericState,
        max_threads,
        lookup_sw::LookUpSW,
        lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing},
    )

`Two Stream` RTE solver for the shortwave problem, using RRTMGP optics.
"""
solve_sw!(
    (; context, fluxb, flux, src, bcs, op)::TwoStreamSWRTE,
    as::AtmosphericState,
    max_threads,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing},
) = rte_sw_2stream_solve!(context.device, fluxb, flux, op, bcs, src, max_threads, as, lookup_sw, lookup_sw_cld)

end
