"""
    GasConcentrations

Encapsulates a collection of volume mixing ratios (concentrations) of gases.
 Each concentration is associated with a name, normally the chemical formula.

Values may be provided as scalars, 1-dimensional profiles (nlay), or 2-D fields (ncol,nlay).
 (nlay and ncol are determined from the input arrays; self-consistency is enforced)
 example:
 set_vmr!(gas_concs, h2o(), values(:,:))
 set_vmr!(gas_concs, o3() , values(:)  )
 set_vmr!(gas_concs, co2(), value      )

Values can be requested as profiles (valid only if there are no 2D fields present in the object)
 or as 2D fields. Values for all columns are returned although the entire collection
 can be subsetted in the column dimension

Subsets can be extracted in the column dimension
"""
module GasConcentrations

using DocStringExtensions
using ..Utilities
using ..Gases
export GasConcs, GasConcsPGP
export set_vmr!

"""
    GasConcs{FT}

Encapsulates a collection of volume mixing ratios (concentrations) of gases.
 Each concentration is associated with a name, normally the chemical formula.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GasConcs{FT<:AbstractFloat,I<:Int}
  "gas names"
  gas_names::Vector{AbstractGas}
  "gas concentrations arrays"
  concs::Array{FT,3}
  "number of columns"
  ncol::I
  "number of layers"
  nlay::I
end
function GasConcs(::Type{FT}, ::Type{I}, gas_names, ncol::I, nlay::I, ngases::I=length(gas_names)) where {FT<:AbstractFloat,I<:Int}
  concs = zeros(FT, ngases, ncol, nlay)
  return GasConcs{FT,I}(gas_names, concs, ncol, nlay)
end

"""
    GasConcsPGP{FT}

Encapsulates a collection of volume mixing ratios (concentrations) of gases per grid point.
Each concentration is associated with a name, normally the chemical formula.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GasConcsPGP{FT<:AbstractFloat}
  "gas names"
  gas_names::Vector{AbstractGas}
  "gas concentrations arrays"
  concs::Array{FT,1}
end

function Base.convert(::Type{GasConcs}, data::Array{GasConcsPGP{FT}}) where {FT}
  s = size(data)
  ngases = length(first(data).concs)
  return GasConcs{FT,typeof(ngases)}(first(data).gas_names,
                      Array{FT}([data[i,j].concs[k] for k in 1:ngases, i in 1:s[1], j in 1:s[2]]), s...)
end

Base.convert(::Type{Array{GasConcsPGP}}, data::GasConcs{FT}) where {FT} =
  [GasConcsPGP{FT}(data.gas_names, data.concs[:,i,j]) for i in 1:size(data.concs,2), j in 1:size(data.concs,3)]


"""
    set_vmr!(this::GasConcs{FT}, gas::AbstractGas, w) where FT

Set volume mixing ratio (vmr)
"""
function set_vmr!(this::GasConcs{FT}, gas::AbstractGas, w::FT) where FT
  @assert !(w < FT(0) || w > FT(1))
  igas = loc_in_array(gas, this.gas_names)
  this.concs[igas,:,:] .= w
  this.gas_names[igas] = gas
end
function set_vmr!(this::GasConcs{FT}, gas::AbstractGas, w::Vector{FT}) where FT
  @assert !any(w .< FT(0)) || any(w .> FT(1))
  @assert !(this.nlay ≠ nothing && length(w) ≠ this.nlay)
  igas = loc_in_array(gas, this.gas_names)
  @assert igas ≠ -1 # assert gas is found
  this.concs[igas,:,:] .= reshape(w, 1, this.nlay)
  this.gas_names[igas] = gas
end
function set_vmr!(this::GasConcs{FT}, gas::AbstractGas, w::Array{FT, 2}) where FT
  @assert !any(w .< FT(0)) || any(w .> FT(1))
  @assert !(this.ncol ≠ nothing && size(w, 1) ≠ this.ncol)
  @assert !(this.nlay ≠ nothing && size(w, 2) ≠ this.nlay)
  igas = loc_in_array(gas, this.gas_names)
  @assert igas ≠ -1 # assert gas is found
  this.concs[igas,:,:] .= w
  this.gas_names[igas] = gas
end

end # module
