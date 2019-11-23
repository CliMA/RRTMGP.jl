module mo_simple_netcdf

using NCDatasets

export read_field

function read_field(ds, varName, args...)
  @assert haskey(ds, varName)
  field = ds[varName][:]
  return field
end

end # module