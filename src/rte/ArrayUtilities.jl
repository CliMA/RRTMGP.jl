"""
    ArrayUtilities

Utility array functions
"""
module ArrayUtilities

export any_vals_less_than, any_vals_outside, reshape_for_comparison

# Values less than a floor
any_vals_less_than(array, check_value) = minimum(array) < check_value
function any_vals_less_than(array, mask, check_value)
  temp = array[mask]
  return isempty(temp) ? false : minimum(temp) < check_value
end

# Values outside a range
function any_vals_outside(array, checkMin, checkMax)
  minValue = minimum(array)
  maxValue = maximum(array)
  return minValue < checkMin || maxValue > checkMax
end
function any_vals_outside(array, mask, checkMin, checkMax)
  temp = array[mask]
  if isempty(temp)
    return false
  end
  minValue = minimum(temp)
  maxValue = maximum(temp)
  return minValue < checkMin || maxValue > checkMax
end

function reshape_for_comparison(flux::Array{FT}, nlay::I, ncol::I, nexp::I) where {FT, I}
  temp = Array{FT}(undef,size(flux,2),size(flux,1),size(flux,3))
  for i = 1:size(flux,3)
    temp[:,:,i] = transpose(flux[:,:,i])
  end
  return reshape(temp,nlay+1,ncol,nexp)
end

end
