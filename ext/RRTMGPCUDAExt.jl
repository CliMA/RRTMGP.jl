module RRTMGPCUDAExt

import ClimaComms
import RRTMGP.Parameters as RP
import RRTMGP.AngularDiscretizations.AngularDiscretization
import RRTMGP.Fluxes: FluxLW, FluxSW
import RRTMGP.Fluxes: compute_net_flux!, set_flux_to_zero!
import RRTMGP.Fluxes: set_band_flux_to_zero!, accumulate_band_flux!
import RRTMGP.Sources: SourceLWNoScat
import RRTMGP.Sources: SourceLW2Str
import RRTMGP.Sources: SourceSW2Str
import RRTMGP.BCs: LwBCs, SwBCs
import RRTMGP.Optics: OneScalar, TwoStream
import RRTMGP.Optics
import RRTMGP.Optics: compute_optical_props!
import RRTMGP.Optics: add_cloud_optics_lw!, add_cloud_optics_sw!
import RRTMGP.Optics: compute_col_gas!
import RRTMGP.Optics: compute_col_gas_kernel!
import RRTMGP.Optics: compute_relative_humidity!
import RRTMGP.Optics: compute_relative_humidity_kernel!
import RRTMGP.AtmosphericStates: setup_gray_as_pr_grid!
import RRTMGP.AtmosphericStates: setup_gray_as_pr_grid_kernel!
import RRTMGP.AtmosphericStates: GrayAtmosphericState
import RRTMGP.AtmosphericStates: AtmosphericState
import RRTMGP.AtmosphericStates: CloudState
import RRTMGP.AtmosphericStates: AerosolState
import RRTMGP.AtmosphericStates: TransposedStateCache
import RRTMGP.AtmosphericStates
import RRTMGP.GrayAtmosphere: update_profile_lw!
import RRTMGP.GrayAtmosphere: compute_gray_heating_rate!
import RRTMGP.GrayAtmosphere: compute_gray_heating_rate_kernel!
import RRTMGP.GrayAtmosphere: update_profile_lw_kernel!
import RRTMGP.LookUpTables:
    LookUpLW, LookUpCld, LookUpPrecip, LookUpSW, LookUpAerosolMerra
import RRTMGP.RTESolver: rte_lw_noscat_solve!
import RRTMGP.RTESolver: rte_lw_noscat_one_angle!
import RRTMGP.RTESolver: rte_lw_2stream_solve!
import RRTMGP.RTESolver: rte_lw_2stream!
import RRTMGP.RTESolver: rte_sw_2stream_solve!
import RRTMGP.RTESolver: rte_sw_2stream!
import RRTMGP.RTESolver: rte_sw_noscat!
import RRTMGP.RTESolver: rte_sw_noscat_solve!
# device-agnostic per-(g-point, column) bodies shared with the CPU drivers
import RRTMGP.RTESolver: lw_noscat_gpt_col!, lw_2stream_gpt_col!
import RRTMGP.RTESolver: sw_noscat_gpt_col!, sw_2stream_gpt_col!
import RRTMGP.RTESolver: _compute_aero_mask!, _build_cloud_mask!
import RRTMGP.RTESolver: rte_lw_2stream_solve_both_skies!
import RRTMGP.RTE: LWGptGroups, SWGptGroups, n_gpt_groups
import RRTMGP.RTE: GPT_GROUP_MEMORY_DIVISOR, _gpt_group_byte_bound
import RRTMGP.RTESolver: lw_2stream_gpt_col_both_skies!
import RRTMGP.RTESolver: rte_sw_2stream_solve_both_skies!
import RRTMGP.RTESolver: sw_2stream_gpt_col_both_skies!
import Adapt
import CUDA
import CUDA: threadIdx, blockIdx, blockDim, @cuda, i32

# `Int`, because CUDA reports memory as a `UInt64`.
_gpt_group_byte_bound(::ClimaComms.CUDADevice) =
    fld(Int(CUDA.totalmem(CUDA.device())), GPT_GROUP_MEMORY_DIVISOR)

_max_threads_cuda() = 256

function _configure_threadblock(max_threads, nitems)
    nthreads = min(max_threads, nitems)
    nblocks = cld(nitems, nthreads)
    return (nthreads, nblocks)
end

_configure_threadblock(nitems) =
    _configure_threadblock(_max_threads_cuda(), nitems)

"""
    Int32Array(array)

Wrapper of a device array whose linear indices are computed with 32-bit integer
arithmetic (64-bit integers take two registers and several instructions on
GPUs), which limits it to `typemax(Int32)` elements. Bounds are only checked for
the linear indices, by the device array.
"""
struct Int32Array{T, N, A <: AbstractArray{T, N}} <: AbstractArray{T, N}
    array::A
end
Adapt.adapt_structure(to, a::Int32Array) = Int32Array(Adapt.adapt(to, a.array))
Base.size(a::Int32Array) = size(a.array)
Base.@propagate_inbounds Base.getindex(a::Int32Array, I::Integer...) =
    a.array[linear_index(a, I)]
Base.@propagate_inbounds Base.setindex!(a::Int32Array, v, I::Integer...) =
    (a.array[linear_index(a, I)] = v)
# The offset from the first element is computed in 32 bits, but the 1 is added in
# 64 bits, so that it cancels the 1 subtracted from the index of the device array
# (otherwise, a pointer to the element before the first is kept for every array).
@inline function linear_index(a, I)
    offset =
        length(I) == 1 ? only(I) % Int32 - Int32(1) :
        offset_32(size(a.array), I)
    return Int(offset) + 1
end
@inline offset_32(_, ::Tuple{}) = Int32(0)
@inline offset_32(dims, I) =
    first(I) % Int32 - Int32(1) +
    first(dims) % Int32 * offset_32(Base.tail(dims), Base.tail(I))

struct Int32Arrays end
Adapt.adapt_storage(::Int32Arrays, a::CUDA.CuArray) =
    length(a) <= typemax(Int32) ? Int32Array(a) :
    throw(ArgumentError("GPU kernels support at most 2^31 - 1 array elements"))
int32_arrays(args) = Adapt.adapt(Int32Arrays(), args)

include(joinpath("cuda", "gray_atmospheric_states.jl"))
include(joinpath("cuda", "rte_longwave_2stream.jl"))
include(joinpath("cuda", "optics.jl"))
include(joinpath("cuda", "rte_shortwave_noscat.jl"))
include(joinpath("cuda", "gray_atmosphere.jl"))
include(joinpath("cuda", "rte_shortwave_2stream.jl"))
include(joinpath("cuda", "rte_longwave_noscat.jl"))
include(joinpath("cuda", "rte_gpt_groups.jl"))

end
