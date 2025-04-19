"""
    AbstractIndexOrder

An abstract type representing index orderings for multidimensional arrays used in radiation data layouts.
"""
abstract type AbstractIndexOrder end

"""
    NVCOrder <: AbstractIndexOrder

Index order with dimensions (N, vertical level, column).
"""
struct NVCOrder <: AbstractIndexOrder end

"""
    VCOrder <: AbstractIndexOrder

Index order with dimensions (vertical level, column).
"""
struct VCOrder <: AbstractIndexOrder end

"""
    NCOrder <: AbstractIndexOrder

Index order with dimensions (N, column).
"""
struct NCOrder <: AbstractIndexOrder end

"""
    NOrder <: AbstractIndexOrder

Index order with a single N dimension.
"""
struct NOrder <: AbstractIndexOrder end

Base.ndims(::Type{NVCOrder}) = 3
Base.ndims(::Type{VCOrder}) = 2
Base.ndims(::Type{NCOrder}) = 2
Base.ndims(::Type{NOrder}) = 1

"""
    RRTMGPData{O, A, T, N} <: AbstractArray{T, N}

A wrapper around a multidimensional array for storing radiation-related data,
with a parameterized index ordering type `O`.

# Fields

- `array::A`: The underlying array storage.

# Example

```julia
RRTMGPData{NCOrder}(rand(2,2))
```
"""
struct RRTMGPData{O <: AbstractIndexOrder, A, T, N} <: AbstractArray{T, N}
    array::A
end

function RRTMGPData{O}(array::A) where {O, A <: AbstractArray}
    @assert ndims(O) == ndims(array)
    RRTMGPData{O, A, eltype(array), ndims(array)}(array)
end

"""
    NVCData(gp::RRTMGPGridParams, ::Type{T}; N, nlay)
    NVCData(gp::RRTMGPGridParams; N, nlay)

Construct a `RRTMGPData` object with index order (N, vertical level, column).
"""
NVCData(gp::RRTMGPGridParams, ::Type{T}; N = 1, nlay) where {T} = RRTMGPData(NVCOrder, gp, T, N, nlay, gp.ncol)

NVCData(gp::RRTMGPGridParams; N = 1, nlay) = NVCData(gp, eltype(gp); N, nlay)

"""
    VCData(gp::RRTMGPGridParams, ::Type{T}; nlay)
    VCData(gp::RRTMGPGridParams; nlay)
    VCData(fn, gp::RRTMGPGridParams; nlay)

Construct a `RRTMGPData` object with index order (vertical level, column).
"""
VCData(gp::RRTMGPGridParams, ::Type{T}; nlay) where {T} = RRTMGPData(VCOrder, gp, T, nlay, gp.ncol)

VCData(gp::RRTMGPGridParams; nlay) = VCData(gp, eltype(gp); nlay)

"""
    NCData(gp::RRTMGPGridParams, ::Type{T}; N)
    NCData(gp::RRTMGPGridParams; N)

Construct a `RRTMGPData` object with index order (N, column).
"""
NCData(gp::RRTMGPGridParams, ::Type{T}; N = 1) where {T} = RRTMGPData(NCOrder, gp, T, N, gp.ncol)
NCData(gp::RRTMGPGridParams; N = 1) = NCData(gp, eltype(gp); N)

"""
    NData(gp::RRTMGPGridParams, ::Type{T}; N)
    NData(gp::RRTMGPGridParams; N)

Construct a `RRTMGPData` object with index order (N).
"""
NData(gp::RRTMGPGridParams, ::Type{T}; N = 1) where {T} = RRTMGPData(NOrder, gp, T, N)
NData(gp::RRTMGPGridParams; N = 1) = NData(gp, eltype(gp); N)

"""
    RRTMGPData(
    	::Type{<:AbstractIndexOrdering},
    	grid_params::RRTMGPGridParams,
    	::Type{Eltype},
    	s...
	)

Construct a `RRTMGPData` object of type `O` and element type `Eltype`, with
dimensions given by `s...`. The data is allocated using the
`ClimaComms.array_type` appropriate for the grid.
"""
function RRTMGPData(::Type{O}, grid_params::RRTMGPGridParams, ::Type{Eltype}, s...) where {O, Eltype}
    DA = ClimaComms.array_type(grid_params)
    return RRTMGPData{O}(DA{Eltype}(undef, s...))
end

Base.parent(data::RRTMGPData) = data.array
Base.size(data::RRTMGPData) = Base.size(parent(data))
Base.length(data::RRTMGPData) = Base.length(parent(data))
Base.eltype(data::RRTMGPData) = Base.eltype(parent(data))

Base.Array(data::RRTMGPData{O}) where {O} = RRTMGPData{O}(Base.Array(parent(data)))
Base.copyto!(data::RRTMGPData, array::AbstractArray) = Base.copyto!(parent(data), array)

Base.@propagate_inbounds Base.getindex(data::RRTMGPData, I::CartesianIndex) = Base.getindex(parent(data), I.I...)

Base.@propagate_inbounds Base.getindex(data::RRTMGPData, indices...) = Base.getindex(parent(data), indices...)

Base.@propagate_inbounds Base.setindex!(data::RRTMGPData, X, I::CartesianIndex) =
    Base.setindex!(parent(data), X, I.I...)

Base.@propagate_inbounds Base.setindex!(data::RRTMGPData, X, indices...) = Base.setindex!(parent(data), X, indices...)

Base.fill!(data::RRTMGPData, value::Number) = fill!(parent(data), value)

#####
##### Custom operators
#####

"""
	set_domain!(data::RRTMGPData, value, gp::RRTMGPGridParams)

Sets `data` in the domain to `value` (excludes extra layer).
"""
function set_domain!(data::RRTMGPData, value, gp::RRTMGPGridParams)
    set_cols!(domain_view(gp.isothermal_boundary_layer, data), value)
    return data
end

"""
	set_cols!(data, value::Union{Number, AbstractArray})

Sets the columns in data to values in value, based on the index order.
"""
function set_cols! end

set_cols!(data::RRTMGPData, value::Real) = fill!(data, value)

# TODO: specialize this more carefully for each ordering
function set_cols!(data::RRTMGPData, value::AbstractArray{<:Real})
    if ndims(data) == 2 && ndims(value) == 1 && size(data, 1) == 1
        copyto!(view(parent(data), 1, :), value)
    elseif ndims(data) == 2
        if size(value) == size(data)
            copyto!(data, value)
        elseif size(value) == (size(data, 1),)
            for col in eachcol(data)
                copyto!(col, value)
            end
        elseif size(value) == (1, size(data, 2))
            for (icol, col) in enumerate(eachcol(data))
                fill!(col, value[1, icol])
            end
        else
            error("expected RRTMGPData to have an array of size $(size(data)), \
                   ($(size(data, 1)),), or (1, $(size(data, 2))); received \
                   an array of size $(size(value))")
        end
    elseif ndims(data) == 3 && ndims(value) == 2 && size(data, 1) == 1
        copyto!(view(parent(data), 1, :, :), value)
    else
        if size(value) == size(data)
            copyto!(data, value)
        else
            error("expected lhs to be an array of size $(size(data)); \
                   received an array of size $(size(value))")
        end
    end
end


index_order_type(::RRTMGPData{O}) where {O} = O

import Adapt
Adapt.adapt_structure(to, data::RRTMGPData) = RRTMGPData{index_order_type(data)}(Adapt.adapt(to, parent(data)))
