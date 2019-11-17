"""
    mo_gas_concentrations

Encapsulates a collection of volume mixing ratios (concentrations) of gases.
 Each concentration is associated with a name, normally the chemical formula.

Values may be provided as scalars, 1-dimensional profiles (nlay), or 2-D fields (ncol,nlay).
 (nlay and ncol are determined from the input arrays; self-consistency is enforced)
 example:
 set_vmr!(gas_concs, "h2o", values(:,:))
 set_vmr!(gas_concs, "o3" , values(:)  )
 set_vmr!(gas_concs, "co2", value      )

Values can be requested as profiles (valid only if there are no 2D fields present in the object)
 or as 2D fields. Values for all columns are returned although the entire collection
 can be subsetted in the column dimension

Subsets can be extracted in the column dimension
"""
module mo_gas_concentrations

using DocStringExtensions
using ..fortran_intrinsics
using ..Utilities
export ty_gas_concs
export set_vmr!, get_vmr!
export GasConcSize

const GAS_NOT_IN_LIST = -1

struct GasConcSize{N,I}
  ncol::I
  nlay::I
  s::NTuple{N,I}
  nconcs::I
end

struct conc_field{FT}
  conc::Array{FT,2}
end

"""
    ty_gas_concs{FT}

Encapsulates a collection of volume mixing ratios (concentrations) of gases.
 Each concentration is associated with a name, normally the chemical formula.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ty_gas_concs{FT}
  gas_name::Vector{String}
  concs
  ncol#::Int
  nlay#::Int
  gsc::GasConcSize
end

function ty_gas_concs(::Type{FT}, gas_names, ncol, nlay, gsc::GasConcSize) where FT
  concs = [conc_field(Array{FT}(undef, gsc.s...)) for i in 1:gsc.nconcs]
  return ty_gas_concs{FT}(gas_names, concs, gsc.s..., gsc)
end

# Set concentrations --- scalar, 1D, 2D
function set_vmr!(this::ty_gas_concs{FT}, gas::String, w::FT) where FT
  @assert !(w < FT(0) || w > FT(1))
  igas = loc_in_array(gas, this.gas_name)
  this.concs[igas].conc .= w
  this.gas_name[igas] = gas
end

function set_vmr!(this::ty_gas_concs{FT}, gas::String, w::Vector{FT}) where FT
  @assert !any(w .< FT(0)) || any(w .> FT(1))
  @assert !(this.nlay ≠ nothing && length(w) ≠ this.nlay)
  igas = loc_in_array(gas, this.gas_name)
  @assert igas ≠ GAS_NOT_IN_LIST
  this.concs[igas].conc .= reshape(w, 1, this.nlay)
  this.gas_name[igas] = gas
end

function set_vmr!(this::ty_gas_concs, gas::String, w::Array{FT, 2}) where FT
  @assert !any(w .< FT(0)) || any(w .> FT(1))
  @assert !(this.ncol ≠ nothing && size(w, 1) ≠ this.ncol)
  @assert !(this.nlay ≠ nothing && size(w, 2) ≠ this.nlay)
  igas = loc_in_array(gas, this.gas_name)
  @assert igas ≠ GAS_NOT_IN_LIST
  this.concs[igas].conc .= w
  this.gas_name[igas] = gas
end

"""
    get_vmr!(array::Array{FT,1}, this::ty_gas_concs{FT}, gas::String) where FT

Volume mixing ratio ( nlay dependence only)
"""
function get_vmr!(array::AbstractArray{FT,1}, this::ty_gas_concs{FT}, gas::String) where FT
  igas = loc_in_array(gas, this.gas_name)
  @assert igas ≠ GAS_NOT_IN_LIST
  @assert !(this.ncol ≠ nothing && this.ncol ≠ size(array,1))
  @assert !(this.nlay ≠ nothing && this.nlay ≠ size(array,2))
  conc = this.concs[igas].conc

  if size(conc, 1) > 1     # Concentration stored as 2D
    array .= conc
  elseif size(conc, 2) > 1 # Concentration stored as 1D
    array .= reshape(conc[1,:], this.ncol, this.nlay)
  else                     # Concentration stored as scalar
    array .= conc[1,1]
  end
  return nothing
end
function get_vmr!(array::AbstractArray{FT,2}, this::ty_gas_concs{FT}, gas::String) where FT
  igas = loc_in_array(gas, this.gas_name)
  @assert igas ≠ GAS_NOT_IN_LIST
  @assert !(this.ncol ≠ nothing && this.ncol ≠ size(array,1))
  @assert !(this.nlay ≠ nothing && this.nlay ≠ size(array,2))
  conc = this.concs[igas].conc

  if size(conc, 1) > 1                      # Concentration stored as 2D
    array .= conc
  elseif size(this.concs[igas].conc, 2) > 1 # Concentration stored as 1D
    array .= spread(conc[1,:], 1, max(this.ncol, size(array,1)))
  else                                      # Concentration stored as scalar
    array .= conc[1,1]
  end
  return nothing
end

end # module
