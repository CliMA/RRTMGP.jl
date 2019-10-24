module mo_simple_netcdf

using ..fortran_intrinsics
using NCDatasets

export read_field, write_field, var_exists

function read_field(ds, varName, args...)
  @assert haskey(ds, varName)
  field = ds[varName][:]
  return field
end

function write_field(ds, varName, var)
  ndim = length( size(var) )
  if ndim == 1
    defVar(ds,varName, var, ("dim1"))
  elseif ndim == 2
    defVar(ds,varName, var, ("dim1","dim2"))
  elseif ndim == 3
    defVar(ds,varName, var, ("dim1","dim2","dim3"))
  elseif ndim == 4
    defVar(ds,varName, var, ("dim1","dim2","dim3","dim4"))
  else
    error("write_field: variables with more than 4 dimensions not supported at this point" * trim(varName))
  end
end

var_exists(ds, varName) = haskey(ds,varName)

end # module