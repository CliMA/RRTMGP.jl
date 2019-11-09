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
export ty_gas_concs
export set_vmr!, get_vmr

const GAS_NOT_IN_LIST = -1

mutable struct conc_field{FT}
  conc::Array{FT,2}
end

"""
    ty_gas_concs{FT}

Encapsulates a collection of volume mixing ratios (concentrations) of gases.
 Each concentration is associated with a name, normally the chemical formula.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct ty_gas_concs{FT}
  gas_name::Vector{String}
  concs
  ncol#::Int
  nlay#::Int
end

function ty_gas_concs(::Type{FT}, gas_names, ncol, nlay) where FT
  gas_name = Vector([])
  concs = Vector([])
  ncol, nlay = nothing, nothing
  return ty_gas_concs{FT}(gas_name, concs, ncol, nlay)
end

# Set concentrations --- scalar, 1D, 2D
function set_vmr!(this::ty_gas_concs{FT}, gas::String, w, sw) where FT
  if length(sw)==1 && length(w)==1
    set_vmr!(this, gas, w[1])
  elseif length(sw)==1 && length(w)>1
    set_vmr!(this, gas, w)
  else
    set_vmr!(this, gas, w)
  end
end

function set_vmr!(this::ty_gas_concs{FT}, gas::String, w::FT) where FT # result(error_msg)
  # real(FT),            intent(in   ) :: w
  @assert !(w < FT(0) || w > FT(1))

  igas = find_gas(this, gas)

  conc = Array{FT}(undef, 1,1)
  fill!(conc, w)
  gas_name = trim(gas)

  if igas == GAS_NOT_IN_LIST
    push!(this.concs, conc_field(conc))
    push!(this.gas_name, gas_name)
    igas = length(this.concs)
  end
  this.concs[igas].conc = conc
  this.gas_name[igas] = gas_name
end

function set_vmr!(this::ty_gas_concs{FT}, gas::String, w::Vector{FT}) where FT
  # real(FT), dimension(:),intent(in   ) :: w

  @assert !any(w .< FT(0)) || any(w .> FT(1))
  @assert !(this.nlay ≠ nothing && length(w) ≠ this.nlay)

  this.nlay = length(w)

  igas = find_gas(this, gas)
  conc = reshape(w, 1, this.nlay)
  gas_name = trim(gas)
  if igas == GAS_NOT_IN_LIST || igas == length(this.concs)+1
    push!(this.concs, conc_field(conc))
    igas = length(this.concs)
  end
  this.concs[igas].conc = conc
  this.gas_name[igas] = gas_name
end

function set_vmr!(this::ty_gas_concs, gas::String, w::Array{FT}) where FT
  # real(FT), dimension(:,:),intent(in   ) :: w

  @assert !any(w .< FT(0)) || any(w .> FT(1))

  @assert !(this.ncol ≠ nothing && size(w, 1) ≠ this.ncol)
  this.ncol = size(w, 1)

  @assert !(this.nlay ≠ nothing && size(w, 2) ≠ this.nlay)

  this.nlay = size(w, 2)

  conc = w
  gas_name = trim(gas)
  igas = find_gas(this, gas)

  if igas == GAS_NOT_IN_LIST
    push!(this.concs, conc_field(conc))
    push!(this.gas_name, gas_name)
    igas = length(this.concs)
  end
  this.concs[igas].conc = conc
  this.gas_name[igas] = gas_name
end

# Return volume mixing ratio as 1D or 2D array
# 1D array ( lay depdendence only)
#
function get_vmr(this::ty_gas_concs{FT}, gas::String) where FT #result(error_msg)
  # real(FT), dimension(:),   intent(out) :: array

  igas = find_gas(this, gas)
  @assert igas ≠ GAS_NOT_IN_LIST
  conc = this.concs[igas].conc

  this.ncol == nothing && (this.ncol = 1)

  array = Array{FT}(undef, this.ncol, this.nlay)

  if size(conc, 1) > 1     # Concentration stored as 2D
    array .= conc[:,:]
  elseif size(conc, 2) > 1 # Concentration stored as 1D
    array .= reshape(conc[1,:], this.ncol, this.nlay)
  else                     # Concentration stored as scalar
    fill!(array, conc[1,1])
  end

  @assert !(this.ncol ≠ nothing && this.ncol ≠ size(array,1))
  @assert !(this.nlay ≠ nothing && this.nlay ≠ size(array,2))

  return array

end

# 2D array (col, lay)
function get_vmr(this::ty_gas_concs, gas::String, array::Array{FT,2}) where FT
  # real(FT), dimension(:,:), intent(out) :: array

  igas = find_gas(this, gas)
  @assert igas ≠ GAS_NOT_IN_LIST

  # Is the requested array the correct size?
  @assert !(this.ncol ≠ nothing && this.ncol ≠ size(array,1))
  @assert !(this.nlay ≠ nothing && this.nlay ≠ size(array,2))

  if size(this.concs[igas].conc, 1) > 1      # Concentration stored as 2D
    array = this.concs[igas].conc[:,:]
  elseif size(this.concs(igas).conc, 2) > 1 # Concentration stored as 1D
    array = spread(this.concs(igas).conc(1,:), dim=1, ncopies=max(this.ncol, size(array,1)))
  else                                                   # Concentration stored as scalar
    array = this.concs[igas].conc[1,1]
  end
  return array

end

# Extract a subset of n columns starting with column 'start'
function get_subset_range(this::ty_gas_concs, start::I, n::I, subset::ty_gas_concs) where I
  # class(ty_gas_concs),      intent(in   ) :: this
  # integer,                  intent(in   ) :: start, n
  # class(ty_gas_concs),      intent(inout) :: subset
  @assert !(n <= 0)
  @assert !(start < 1 )
  @assert !(this.ncol ≠ nothing && start > this.ncol || start+n-1 > this.ncol )

  this.nlay = 0
  this.ncol = 0
  # reset!(subset)
  # These two arrays should be the same length
  subset.gas_name = Array(undef, size(this.gas_name))
  subset.concs = Array(undef, size(this.concs))
  subset.nlay = this.nlay
  subset.ncol = fmerge(n, 0, this.ncol ≠ nothing)
  subset.gas_name[:] .= this.gas_name[:]

  for i = 1:size(this.gas_name)
    #
    # Preserve scalar/1D/2D representation in subset,
    #   but need to ensure at least extent 1 in col dimension (ncol = 0 means no gas exploits this dimension)
    #
    allocate(subset.concs[i].conc(min(max(subset.ncol,1), size(this.concs[i].conc, 1)),
                                        min(    subset.nlay,    size(this.concs[i].conc, 2))))
    if size(this.concs[i].conc, 1) > 1      # Concentration stored as 2D
      subset.concs[i].conc[:,:] .= this.concs[i].conc[start:(start+n-1),:]
    else
      subset.concs[i].conc[:,:] .= this.concs[i].conc[:,:]
    end
  end

end

# Inquiry functions
get_num_gases(this::ty_gas_concs) = size(this.gas_name)

get_gas_names(this) = this.gas_name[:]

# find gas in list; GAS_NOT_IN_LIST if not found
function find_gas(this::ty_gas_concs, gas)
  L = lowercase.(strip.(this.gas_name)) .== lowercase(trim(gas))
  !any(L) && return GAS_NOT_IN_LIST
  return argmax(L)
end

end # module
