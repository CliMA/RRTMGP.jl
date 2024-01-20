module RTESolver

using StaticArrays
import ClimaComms
using CUDA
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
    solve_lw!(
        slv::Solver,
        max_threads::Int,
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    )

Solver for the longwave radiation problem
"""
solve_lw!(slv::Solver, max_threads, lookup_lw = nothing, lookup_lw_cld = nothing) =
    solve_lw!(slv.context.device, slv, max_threads, lookup_lw, lookup_lw_cld)

# OneScalar solvers with gray optics
solve_lw!(
    device::ClimaComms.AbstractDevice,
    (; flux_lw, src_lw, bcs_lw, op, as)::Solver{<:Any, <:GrayAtmosphericState, <:OneScalar},
    max_threads::Int,
    ::Nothing,
    ::Nothing,
) = rte_lw_noscat_solve!(device, flux_lw, src_lw, bcs_lw, op, max_threads, as)

# 2-stream solvers with gray optics
solve_lw!(
    device::ClimaComms.AbstractDevice,
    (; flux_lw, src_lw, bcs_lw, op, as)::Solver{<:Any, <:GrayAtmosphericState, <:TwoStream},
    max_threads::Int,
    ::Nothing,
    ::Nothing,
) = rte_lw_2stream_solve!(device, flux_lw, src_lw, bcs_lw, op, max_threads, as)

# OneScalar solvers with RRTMGP optics
solve_lw!(
    device::ClimaComms.AbstractDevice,
    (; fluxb_lw, flux_lw, src_lw, bcs_lw, op, as)::Solver{<:Any, <:AtmosphericState, <:OneScalar},
    max_threads::Int,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing},
) = rte_lw_noscat_solve!(device, fluxb_lw, flux_lw, src_lw, bcs_lw, op, max_threads, as, lookup_lw, lookup_lw_cld)

# 2-stream solvers with RRTMGP optics
solve_lw!(
    device::ClimaComms.AbstractDevice,
    (; fluxb_lw, flux_lw, src_lw, bcs_lw, op, as)::Solver{<:Any, <:AtmosphericState, <:TwoStream},
    max_threads::Int,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing},
) = rte_lw_2stream_solve!(device, fluxb_lw, flux_lw, src_lw, bcs_lw, op, max_threads, as, lookup_lw, lookup_lw_cld)

"""
    solve_sw!(
        slv::Solver,
        max_threads::Int,
        lookup_sw::Union{LookUpSW, Nothing} = nothing,
        lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    )

Solver for the shortwave radiation problem
"""
solve_sw!(
    slv::Solver,
    max_threads::Int,
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) = solve_sw!(slv.context.device, slv, max_threads, lookup_sw, lookup_sw_cld)

# non-scattering solver, gray optics
solve_sw!(
    device::ClimaComms.AbstractDevice,
    (; flux_sw, op, bcs_sw, as)::Solver{<:Any, <:GrayAtmosphericState, <:OneScalar},
    max_threads::Int,
    ::Nothing,
    ::Nothing,
) = rte_sw_noscat_solve!(device, flux_sw, op, bcs_sw, max_threads, as) # non-scattering solver, gray optics

# 2-stream solver, gray optics
solve_sw!(
    device::ClimaComms.AbstractDevice,
    (; flux_sw, op, bcs_sw, src_sw, as)::Solver{<:Any, <:GrayAtmosphericState, <:TwoStream},
    max_threads::Int,
    ::Nothing,
    ::Nothing,
) = rte_sw_2stream_solve!(device, flux_sw, op, bcs_sw, src_sw, max_threads, as)

# non-scattering solver, RRTMGP optics
solve_sw!(
    device::ClimaComms.AbstractDevice,
    (; fluxb_sw, flux_sw, op, bcs_sw, as)::Solver{<:Any, <:AtmosphericState, <:OneScalar},
    max_threads::Int,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) = rte_sw_noscat_solve!(device, fluxb_sw, flux_sw, op, bcs_sw, max_threads, as, slv.lookup_sw, slv.lookup_sw_cld)

# 2-stream solver, RRTMGP optics
solve_sw!(
    device::ClimaComms.AbstractDevice,
    (; fluxb_sw, flux_sw, op, bcs_sw, src_sw, as)::Solver{<:Any, <:AtmosphericState, <:TwoStream},
    max_threads::Int,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
) = rte_sw_2stream_solve!(device, fluxb_sw, flux_sw, op, bcs_sw, src_sw, max_threads, as, lookup_sw, lookup_sw_cld)

end
