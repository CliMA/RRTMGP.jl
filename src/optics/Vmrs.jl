module Vmrs
using DocStringExtensions
using Adapt
import ..RRTMGPGridParams

export AbstractVmr, Vmr, VmrGM, get_vmr

"""
    AbstractVmr{FT}
"""
abstract type AbstractVmr{FT <: AbstractFloat} end
"""
    VmrGM{FT,FTA1D,FTA2D} <: AbstractVmr{FT}

Volume mixing ratios for various gases in the atmosphere. This struct can be used
when only H₂O and O₃ concentrations vary spatially and a global mean is used 
for all other gases.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct VmrGM{
    FT <: AbstractFloat,
    FTA1D <: AbstractArray{FT, 1},
    FTA2D <: AbstractArray{FT, 2},
} <: AbstractVmr{FT}
    "volume mixing ratio of H₂O"
    vmr_h2o::FTA2D
    "volume mixing ratio of Ozone"
    vmr_o3::FTA2D
    "volume mixing ratio of all other gases, which are independent of location and column"
    vmr::FTA1D
end
VmrGM(vmr_h2o, vmr_o3, vmr) =
    VmrGM{eltype(vmr_h2o), typeof(vmr), typeof(vmr_h2o)}(vmr_h2o, vmr_o3, vmr)
Adapt.@adapt_structure VmrGM

"""
    VolumeMixingRatioGlobalMean(grid_params::RRTMGPGridParams; vmr_h2o, vmr_o3, vmr)

Return a [`VmrGM`](@ref) given:
 - `grid_params::RRTMGPGridParams`: grid parameters
 - `vmr_h2o`: `(nlay, ncol)` volume mixing ratio of H₂O
 - `vmr_o3`: `(nlay, ncol)` volume mixing ratio of O₃
 - `vmr`: `(ngas,)` global-mean volume mixing ratios of the well-mixed gases
"""
function VolumeMixingRatioGlobalMean(
    grid_params::RRTMGPGridParams;
    vmr_h2o::AbstractArray{FT, 2},
    vmr_o3::AbstractArray{FT, 2},
    vmr::AbstractArray{FT, 1},
) where {FT}
    @assert FT == eltype(grid_params)
    return VmrGM{FT, typeof(vmr), typeof(vmr_h2o)}(vmr_h2o, vmr_o3, vmr)
end

"""
    Vmr{FT,FTA3D} <: AbstractVmr{FT}

Volume mixing ratios for various gases in the atmosphere. This struct can be used
when concentrations vary spatially for all gases.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Vmr{FT <: AbstractFloat, FTA3D <: AbstractArray{FT, 3}} <:
       AbstractVmr{FT}
    "volume mixing ratio of all gases as a function of location and column"
    vmr::FTA3D
end
Vmr(vmr) = Vmr{eltype(vmr), typeof(vmr)}(vmr)
Adapt.@adapt_structure Vmr
"""
    get_vmr(
        vmr::VmrGM{FT},
        ig::Int,
        ilay::Int,
        icol::Int,
    ) where {FT<:AbstractFloat}

Obtain volume mixing ratio of gas `ig` for layer `ilay` of column `icol`.
"""
@inline function get_vmr(
    vmr::VmrGM{FT},
    ig::Int,
    ilay::Int,
    icol::Int,
) where {FT <: AbstractFloat}
    if ig == 0
        return FT(1)
    elseif ig == 1 # h2o / h2o_foreign / h2o_self-continua
        return @inbounds vmr.vmr_h2o[ilay, icol]
    elseif ig == 3 # ozone
        return @inbounds vmr.vmr_o3[ilay, icol]
    else # other gases
        return @inbounds vmr.vmr[ig]
    end
end

"""
    get_vmr(
        vmr::Vmr{FT},
        ig::Int,
        ilay::Int,
        icol::Int,
    ) where {FT<:AbstractFloat}

Obtain volume mixing ratio of gas `ig` for layer `ilay` of column `icol`.
"""
@inline function get_vmr(
    vmr::Vmr{FT},
    ig::Int,
    ilay::Int,
    icol::Int,
) where {FT <: AbstractFloat}
    if ig == 0
        return FT(1)
    else
        return @inbounds vmr.vmr[ig, ilay, icol]
    end
end

end
