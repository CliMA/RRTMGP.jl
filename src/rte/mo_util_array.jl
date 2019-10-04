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

export any_vals_less_than, any_vals_outside, zero_array!

#-------------------------------------------------------------------------------------------------
# Values less than a floor
#-------------------------------------------------------------------------------------------------
any_vals_less_than(array, minVal) = any(array .< minVal)

#-------------------------------------------------------------------------------------------------
# Values outside a range
#-------------------------------------------------------------------------------------------------
any_vals_outside(array, minVal, maxVal) = any(array .< minVal || array .> maxVal)

#-------------------------------------------------------------------------------------------------
# Initializing arrays to 0
#-------------------------------------------------------------------------------------------------
zero_array!(array::Array{FT}) where FT = (array .= FT(0))

end
