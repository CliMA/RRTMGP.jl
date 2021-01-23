module Vmrs
using KernelAbstractions
using CUDA
using ..Device: array_type, array_device
using DocStringExtensions
using Adapt

export Vmr, get_vmr

struct Vmr{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}
    "volume mixing ratio of H₂O"
    vmr_h2o::FTA2D
    "volume mixing ratio of Ozone"
    vmr_o3::FTA2D
    "volume mixing ratio of all other gases, which are independent of location and column"
    vmr::FTA1D
end
Adapt.@adapt_structure Vmr

function get_vmr(
    vmr::Vmr{FT,FTA1D,FTA2D},
    ig::Int,
    ilay::Int,
    icol::Int,
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}
    if ig == 0
        return FT(1)
    elseif ig == 1 # h2o / h2o_foreign / h2o_self-continua
        return vmr.vmr_h2o[ilay, icol]
    elseif ig == 3 # ozone
        return vmr.vmr_o3[ilay, icol]
    else
        return vmr.vmr[ig]
    end
end

function Vmr(
    ngas::Int,
    nlay::Int,
    ncol::Int,
    ::Type{FT},
    ::Type{DA},
) where {FT<:AbstractFloat,DA}
    FTA1D = DA{FT,1}
    FTA2D = DA{FT,2}

    return Vmr{FT,FTA1D,FTA2D}(
        FTA2D(zeros(nlay, ncol)),
        FTA2D(zeros(nlay, ncol)),
        FTA1D(zeros(ngas)),
    )
end

end
