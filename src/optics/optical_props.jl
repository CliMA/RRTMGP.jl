"""
    OneScalar{FTA2D,AD} <: AbstractOpticalProps

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OneScalar{D, V} <: AbstractOpticalProps
    "storage for optical thickness"
    layerdata::D
    "view into optical depth"
    τ::V
end

function Adapt.adapt_structure(to, op::OneScalar)
    layerdata = Adapt.adapt(to, op.layerdata)
    τ = view(layerdata, 1, :, :)
    return OneScalar{typeof(layerdata), typeof(τ)}(layerdata, τ)
end

function OneScalar(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    layerdata = DA{FT, 3}(undef, 1, nlay, ncol)
    τ = view(layerdata, 1, :, :)
    V = typeof(τ)
    return OneScalar{typeof(layerdata), V}(layerdata, τ)
end


"""
    TwoStream{FTA2D} <: AbstractOpticalProps

Two stream approximation for optical properties, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStream{D, V} <: AbstractOpticalProps
    "storage for optical depth, single scattering albedo and asymmerty parameter"
    layerdata::D
    "view into optical depth"
    τ::V
    "view into single scattering albedo"
    ssa::V
    "view into asymmetry parameter"
    g::V
end

function Adapt.adapt_structure(to, op::TwoStream)
    layerdata = Adapt.adapt(to, op.layerdata)
    τ = view(layerdata, 1, :, :)
    ssa = view(layerdata, 2, :, :)
    g = view(layerdata, 3, :, :)
    return TwoStream{typeof(layerdata), typeof(τ)}(layerdata, τ, ssa, g)
end

function TwoStream(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    layerdata = DA{FT, 3}(zeros(3, nlay, ncol))
    V = typeof(view(layerdata, 1, :, :))
    return TwoStream{typeof(layerdata), V}(
        layerdata,
        view(layerdata, 1, :, :),
        view(layerdata, 2, :, :),
        view(layerdata, 3, :, :),
    )
end

