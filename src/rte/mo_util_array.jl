# This code is part of Radiative Transfer for Energetics (RTE)
#
# Contacts: Robert Pincus and Eli Mlawer
# email:  rrtmgp@aer.com
#
# Copyright 2015-2019,  Atmospheric and Environmental Research and
# Regents of the University of Colorado.  All right reserved.
#
# Use and duplication is permitted under the terms of the
#    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
# -------------------------------------------------------------------------------------------------
module mo_util_array

export any_vals_less_than, any_vals_outside, zero_array!, reshape_for_comparison

# Values less than a floor
any_vals_less_than(array, check_value) = minimum(array) < check_value
function any_vals_less_than(array, mask, check_value)
  temp = array[mask]
  if isempty(temp)
    return false
  end
  return minimum(temp) < check_value
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

# Initializing arrays to 0
zero_array!(array::Array{FT}) where FT = (array .= FT(0))


function reshape_for_comparison(flux::Array{FT}, nlay::I, ncol::I, nexp::I) where {FT, I}
  temp = Array{FT}(undef,size(flux,2),size(flux,1),size(flux,3))
  for i = 1:size(flux,3)
    temp[:,:,i] = transpose(flux[:,:,i])
  end
  return reshape(temp,nlay+1,ncol,nexp)
end

end
