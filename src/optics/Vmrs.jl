module Vmrs
using CUDA
using DocStringExtensions
using Adapt

export AbstractVmr, Vmr, VmrGM, init_vmr, get_vmr, Vmrc, VmrGMc

"""
    AbstractVmr
"""
abstract type AbstractVmr end

abstract type AbstractVmrc <: AbstractVmr end
"""
    VmrGM{FT,FTA1D,FTA2D} <: AbstractVmr

Volume mixing ratios for various gases in the atmosphere. This struct can be used
when only H₂O and O₃ concentrations vary spatially and a global mean is used 
for all other gases.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct VmrGM{FT <: AbstractFloat, FTA1D <: AbstractArray{FT, 1}, FTA2D <: AbstractArray{FT, 2}} <: AbstractVmr
    "volume mixing ratio of H₂O"
    vmr_h2o::FTA2D
    "volume mixing ratio of Ozone"
    vmr_o3::FTA2D
    "volume mixing ratio of all other gases, which are independent of location and column"
    vmr::FTA1D
end
VmrGM(vmr_h2o, vmr_o3, vmr) = VmrGM{eltype(vmr_h2o), typeof(vmr), typeof(vmr_h2o)}(vmr_h2o, vmr_o3, vmr)
Adapt.@adapt_structure VmrGM

struct VmrGMc{FTA1DV, FTA1D} <: AbstractVmrc
    vmr_h2o::FTA1DV
    vmr_o3::FTA1DV
    vmr::FTA1D
end
Adapt.@adapt_structure VmrGMc

VmrGMc(vmr::VmrGM, gcol) = VmrGMc(view(vmr.vmr_h2o, :, gcol), view(vmr.vmr_o3, :, gcol), vmr.vmr)
get_col_view(vmr::VmrGM, gcol) = VmrGMc(vmr, gcol)

function VmrGM(
    ::Type{FT},
    ::Type{DA},
    vmr_h2o::FTA2D,
    vmr_o3::FTA2D,
    vmr::FTA1D,
) where {FT <: AbstractFloat, FTA1D <: AbstractArray{FT, 1}, FTA2D <: AbstractArray{FT, 2}, DA}
    return VmrGM{FT, typeof(vmr), typeof(vmr_h2o)}(vmr_h2o, vmr_o3, vmr)
end

"""
    Vmr{FT,FTA3D} <: AbstractVmr

Volume mixing ratios for various gases in the atmosphere. This struct can be used
when concentrations vary spatially for all gases.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Vmr{FT <: AbstractFloat, FTA3D <: AbstractArray{FT, 3}} <: AbstractVmr
    "volume mixing ratio of all gases as a function of location and column"
    vmr::FTA3D
end
Vmr(vmr) = Vmr{eltype(vmr), typeof(vmr)}(vmr)
Adapt.@adapt_structure Vmr

struct Vmrc{FTA2D} <: AbstractVmrc
    vmr::FTA2D
end
Adapt.@adapt_structure Vmrc

Vmrc(vmr::Vmr, gcol) = Vmrc(view(vmr.vmr, :, gcol, :))
get_col_view(vmr::Vmr, gcol) = Vmrc(vmr, gcol)

"""
    get_vmr(
        vmr::VmrGM{FT},
        ig::Int,
        ilay::Int,
        icol::Int,
    ) where {FT<:AbstractFloat}

Obtain volume mixing ratio of gas `ig` for layer `ilay` of column `icol`.
"""
@inline function get_vmr(vmr::VmrGM{FT}, ig::Int, ilay::Int, icol::Int) where {FT <: AbstractFloat}
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

@inline function get_vmr(vmr::VmrGMc, ig::Int, ilay::Int)
    if ig == 0
        FT = eltype(vmr.vmr_h2o)
        return FT(1)
    elseif ig == 1 # h2o / h2o_foreign / h2o_self-continua
        return @inbounds vmr.vmr_h2o[ilay]
    elseif ig == 3 # ozone
        return @inbounds vmr.vmr_o3[ilay]
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
@inline function get_vmr(vmr::Vmr{FT}, ig::Int, ilay::Int, icol::Int) where {FT}
    if ig == 0
        return FT(1)
    else
        return @inbounds vmr.vmr[ilay, icol, ig]
    end
end

@inline function get_vmr(vmr::Vmrc, ig::Int, ilay::Int)
    if ig == 0
        FT = eltype(vmr.vmr)
        return FT(1)
    else
        return @inbounds vmr.vmr[ilay, ig]
    end
end

"""
    init_vmr(
        ngas,
        nlay,
        ncol,
        ::Type{FT},
        ::Type{DA};
        gm::Bool = true,
    ) where {FT<:AbstractFloat,DA}

Initialize the Vmr struct.
"""
function init_vmr(
    ngas::Int,
    nlay::Int,
    ncol::Int,
    ::Type{FT},
    ::Type{DA};
    gm::Bool = true,
) where {FT <: AbstractFloat, DA}
    if gm
        FTA1D = DA{FT, 1}
        FTA2D = DA{FT, 2}
        return VmrGM{FT, FTA1D, FTA2D}(FTA2D(zeros(nlay, ncol)), FTA2D(zeros(nlay, ncol)), FTA1D(zeros(ngas)))
    else
        FTA3D = DA{FT, 3}
        return Vmr{FT, FTA3D}(FTA3D(zeros(nlay, ncol, ngas)))
    end
end

end
