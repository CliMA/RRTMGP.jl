"""
    OneScalar{D, V} <: AbstractOpticalProps

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission.

# Fields
- `layerdata`: Storage for optical thickness.
- `τ`: View into optical depth.
"""
struct OneScalar{D, V} <: AbstractOpticalProps
    layerdata::D
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
    TwoStream{D, V} <: AbstractOpticalProps

Two stream approximation for optical properties, used in
calculations accounting for extinction and emission.

# Fields
- `layerdata`: Storage for optical depth, single scattering albedo and asymmetry
  parameter.
- `τ`: View into optical depth.
- `ssa`: View into single scattering albedo.
- `g`: View into asymmetry parameter.
"""
struct TwoStream{D, V} <: AbstractOpticalProps
    layerdata::D
    τ::V
    ssa::V
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

