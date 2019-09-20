# This code is part of
# RRTM for GCM Applications - Parallel (RRTMGP)
#
# Eli Mlawer and Robert Pincus
# Andre Wehe and Jennifer Delamere
# email:  rrtmgp@aer.com
#
# Copyright 2018,  Atmospheric and Environmental Research and
# Regents of the University of Colorado.  All right reserved.
#
# Use and duplication is permitted under the terms of the
#    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
#
# Description: Kernels to permute arrays

module mo_reorder_kernels

export reorder_123x312_kernel, reorder_123x321_kernel

# ----------------------------------------------------------------------------
function reorder_123x312_kernel(d1, d2, d3, array_in, array_out)

  for i2 = 1:d2
    for i1 = 1:d1
      for i3 = 1:d3
        array_out[i3,i1,i2] = array_in[i1,i2,i3]
      end
    end
  end

end
# ----------------------------------------------------------------------------
function reorder_123x321_kernel(d1, d2, d3, array_in, array_out) & 

  for i1 = 1:d1
    for i2 = 1:d2
      for i3 = 1:d3
        array_out[i3,i2,i1] = array_in[i1,i2,i3]
      end
    end
  end

end
# ----------------------------------------------------------------------------

end
