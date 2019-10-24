# This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
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
# Physical and mathematical constants used in gas optics calculation
#   If the host model in which RRTGMP is embedded has defined these constants elsewhere
#   the model definitions can be used instead by renaming. For example,
# use  mo_model_constants, only k_boltz => boltzman_k, ...
#   where the syntax is local_name => original_name
#   and all the local names need to be defined
#
# "Constants" specific to the earth's atmosphere should also be made consistent with the
#   host model but may be changed in a call to init_constants(), normally at initialization
# -------------------------------------------------------------------------------------------------
module mo_rrtmgp_constants
  # use mo_rte_kind, only: FT

  export k_boltz,
         m_h2o,
         avogad,
         R_univ_gconst,
         m_dry,
         grav,
         cp_dry

  # -----------------------------------------
  # Physical constants, 2018 SI defintion of metric system
  #   doi:10.1088/1681-7575/aa950a (see also https://www.nist.gov/si-redefinition/meet-constants)
  # Boltzmann constant [J/K] = [(kg m^2)/(K s^2)]
  k_boltz(::Type{FT}) where FT = FT(1.380649e-23)

  #  molecular weight of water [kg/mol]
  m_h2o(::Type{FT}) where FT = FT(0.018016)

  # Avogadro's number [molec/mol]
  avogad(::Type{FT}) where FT = FT(6.02214076e23)

  # Universal gas constant [J/(mol K)]
  R_univ_gconst(::Type{FT}) where FT = avogad(FT) * k_boltz(FT)

  # -----------------------------------------
  #
  # Constants specific to the earth's atmosphere -- changeable in init() because they
  #   might be different on e.g. other planets

  # molecular weight of dry air [kg/mol]
  m_dry(::Type{FT}) where FT = FT(0.028964)

  # Gravity at Earth's surface [m/s2]
  grav(::Type{FT}) where FT = FT(9.80665)

  # Specific heat at constant pressure for dry air [J/(K kg)]
  cp_dry(::Type{FT}) where FT = FT(1004.64)

end
