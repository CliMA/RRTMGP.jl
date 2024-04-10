
# Borrowed from:
#    https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/rrtmgp-frontend/mo_aerosol_optics_rrtmgp_merra.F90

# This code is part of Radiative Transfer for Energetics (RTE)
#
# Contacts: Robert Pincus and Eli Mlawer
# email:  rrtmgp@aer.com
#
# Copyright 2015-2018,  Atmospheric and Environmental Research and
# Regents of the University of Colorado.  All right reserved.
#
# Use and duplication is permitted under the terms of the
#    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
# -------------------------------------------------------------------------------------------------
# Provides aerosol optical properties as a function of aerosol size (radius), aerosol mass,
# and relative humidity for the RRTMGP spectral bands.
#   Based on climatoligical aerosol optical properties used in MERRA2 as derived from the
#     GOCART model for 15 aerosol types, including dust and sea salt each for five size bins,
#     one sulfate type, and both hydrophobic and hydrophilic black carbon and organic carbon.
#   Input aerosol optical data are stored in look-up tables.
#
#   References for the gocart interactive aerosols:
#     Chin et al., jgr, 2000 (https://doi.org/10.1029/2000jd900384)
#     Chin et al., jas, 2002 (https://doi.org/10.1175/1520-0469(2002)059<0461:TAOTFT>2.0.CO;2)
#     Colarco et al., jgr, 2010 (https://doi.org/10.1029/2009jd012820)
#
#   References for merra2 aerosol reanalysis:
#     Randles et al., j. clim., 2017 (https://doi.org/10.1175/jcli-d-16-0609.1)
#     Buchard et al., j. clim., 2017 (https://doi.org/10.1175/jcli-d-16-0613.1)
#
# The class can be used as-is but is also intended as an example of how to extend the RTE framework

module AerosolOptics

import ..Optics: AbstractOpticalProps
# intrinsics
# TODO: fix these
function freshape(x; shape, order)
  reshape(x, shape, order)
end
maxval(v) = max(v...)
minval(v) = min(v...)
fepsilon(v) = eps(v)
allocated(v) = !isnothing(v)

include("aerosol_optics.jl")

end # module
