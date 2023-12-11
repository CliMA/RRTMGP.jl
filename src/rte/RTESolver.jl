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
function solve_sw!(
    slv::Solver,
    max_threads::Int,
    lookup_sw::Union{LookUpSW, Nothing} = nothing,
    lookup_sw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    (; as, op, bcs_sw, src_sw, flux_sw, fluxb_sw) = slv
    context = RTE.context(slv)

    no_args = lookup_sw isa Nothing && lookup_sw_cld isa Nothing
    if no_args
        DA = RTE.array_type(slv)
        FT = RTE.float_type(slv)
        major_gpt2bnd = DA(UInt8[1])
        solar_src_scaled = DA(FT[1])
        flux = flux_sw
    else
        (; major_gpt2bnd, solar_src_scaled) = lookup_sw
        flux = fluxb_sw
    end

    # solving radiative transfer equation
    if slv.op isa OneScalar
        rte_sw_noscat_solve!(
            context,
            flux,
            flux_sw,
            op,
            bcs_sw,
            solar_src_scaled,
            max_threads,
            as,
            lookup_sw,
            lookup_sw_cld,
        ) # no-scattering solver
    else
        rte_sw_2stream_solve!(
            context,
            flux,
            flux_sw,
            op,
            bcs_sw,
            src_sw,
            major_gpt2bnd,
            solar_src_scaled,
            max_threads,
            as,
            lookup_sw,
            lookup_sw_cld,
        ) # 2-stream solver
    end
    return nothing
end

end
