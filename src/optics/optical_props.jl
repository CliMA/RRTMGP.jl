import ..Fluxes: _coalesced_3d, _coalesced_3d_zeros

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
    τ = view(layerdata, :, :, 1)
    return OneScalar{typeof(layerdata), typeof(τ)}(layerdata, τ)
end

function OneScalar(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    layerdata = _coalesced_3d(DA, FT, ncol, nlay, 1)
    τ = view(layerdata, :, :, 1)
    V = typeof(τ)
    return OneScalar{typeof(layerdata), V}(layerdata, τ)
end


"""
    TwoStream{D, V, AD, A} <: AbstractOpticalProps

Two stream approximation for optical properties, used in
calculations accounting for extinction and emission.

# Fields
- `layerdata`: Storage for optical depth, single scattering albedo and asymmetry
  parameter.
- `τ`: View into optical depth `(ncol, nlay)`.
- `ssa`: View into single scattering albedo `(ncol, nlay)`.
- `g`: View into asymmetry parameter `(ncol, nlay)`.
- `aerodata`: Storage for one band's aerosol optics, or `nothing`.
- `aero`: Views `(τ, ssa, g)` into `aerodata`, each `(ncol, nlay)`, or
  `nothing`. See [`with_aerosol_band_store`](@ref).
"""
struct TwoStream{D, V, AD, A} <: AbstractOpticalProps
    layerdata::D
    τ::V
    ssa::V
    g::V
    aerodata::AD
    aero::A
end

@inline _aerosol_band_views(::Nothing) = nothing
@inline _aerosol_band_views(aerodata) = (;
    τ = view(aerodata, :, :, 1),
    ssa = view(aerodata, :, :, 2),
    g = view(aerodata, :, :, 3),
)

function TwoStream(layerdata, aerodata = nothing)
    τ = view(layerdata, :, :, 1)
    ssa = view(layerdata, :, :, 2)
    g = view(layerdata, :, :, 3)
    aero = _aerosol_band_views(aerodata)
    return TwoStream{
        typeof(layerdata),
        typeof(τ),
        typeof(aerodata),
        typeof(aero),
    }(
        layerdata,
        τ,
        ssa,
        g,
        aerodata,
        aero,
    )
end

Adapt.adapt_structure(to, op::TwoStream) =
    TwoStream(Adapt.adapt(to, op.layerdata), Adapt.adapt(to, op.aerodata))

function TwoStream(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    return TwoStream(_coalesced_3d_zeros(DA, FT, ncol, nlay, 3))
end

"""
    with_aerosol_band_store(op::TwoStream)

`op` with storage for one band's aerosol optics, sharing `op`'s own arrays.

Aerosol optics depend on the band, the column and the layer but not on the
g-point. With this store, `compute_optical_props!` looks them up at the first
g-point of each band and adds the stored values to the band's other g-points,
instead of repeating the lookup for every one. The optics are the same bit for
bit either way. An `op` that already has a store is returned as it is.
"""
with_aerosol_band_store(op::TwoStream) =
    isnothing(op.aerodata) ? TwoStream(op.layerdata, zero(op.layerdata)) : op
