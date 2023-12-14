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
function solve_lw!(
    slv::Solver,
    max_threads::Int,
    lookup_lw::Union{LookUpLW, Nothing} = nothing,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    (; as, op, bcs_lw, src_lw, flux_lw, fluxb_lw) = slv
    context = RTE.context(slv)
    DA = RTE.array_type(slv)
    no_args = lookup_lw isa Nothing && lookup_lw_cld isa Nothing
    if no_args
        major_gpt2bnd = DA(UInt8[1])
        flux = flux_lw
    else
        (; major_gpt2bnd) = lookup_lw
        flux = fluxb_lw
    end

    rte_lw_solve!(context, flux, flux_lw, src_lw, bcs_lw, op, major_gpt2bnd, max_threads, as, lookup_lw, lookup_lw_cld)
    return nothing
end

"""
    solve_sw!(
        slv::Solver,
        max_threads::Int,
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
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
#--------------------------------------------------------------------------------------------------
#=
"""
    solve_sw!(
        slv::Solver,
        max_threads::Int,
        lookup_lw::Union{LookUpLW, Nothing} = nothing,
        lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
    )

Solver for the shortwave radiation problem
"""
function solve_sw!(
    slv::Solver,
    max_threads::Int,
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    (; as, op, bcs_sw, src_sw, flux_sw, fluxb_sw) = slv
    context = RTE.context(slv)
    device = context.device

    no_args = lookup_sw isa Nothing && lookup_sw_cld isa Nothing
    if no_args
        DA = RTE.array_type(slv)
        FT = RTE.float_type(slv)
        flux = flux_sw
    else
        flux = fluxb_sw
    end

    # solving radiative transfer equation
    if slv.op isa OneScalar
        rte_sw_noscat_solve!(device, flux, flux_sw, op, bcs_sw, max_threads, as, lookup_sw, lookup_sw_cld) # no-scattering solver
    else
        rte_sw_2stream_solve!(device, flux, flux_sw, op, bcs_sw, src_sw, max_threads, as, lookup_sw, lookup_sw_cld) # 2-stream solver
    end
    return nothing
end
=#
end
