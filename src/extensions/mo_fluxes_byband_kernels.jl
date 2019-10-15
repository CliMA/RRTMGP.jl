# This code is part of Radiative Transfer for Energetics (RTE)
#
# Eli Mlawer and Robert Pincus
# Andre Wehe and Jennifer Delamere
# email:  rrtmgp@aer.com
#
# Copyright 2015-2018,  Atmospheric and Environmental Research and
# Regents of the University of Colorado.  All right reserved.
#
# Use and duplication is permitted under the terms of the
#    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
#
# Kernels for computing byband fluxes by summing over all elements in the spectral dimension
#
module mo_fluxes_byband_kernels
  # use, intrinsic :: iso_c_binding
  # use mo_rte_kind, only: FT
  # implicit none
  # private
  # public :: sum_byband, net_byband
  export sum_byband, net_byband

  # interface net_byband
  #   module procedure net_byband_full, net_byband_precalc
  # end interface net_byband

# contains
# ----------------------------------------------------------------------------
  #
  # Spectral reduction over all points
  #
  function sum_byband(ncol, nlev, ngpt, nbnd, band_lims, spectral_flux)# bind (C)
    # integer,                               intent(in ) :: ncol, nlev, ngpt, nbnd
    # integer,  dimension(2,          nbnd), intent(in ) :: band_lims
    # real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux
    # real(FT), dimension(ncol, nlev, nbnd), intent(out) :: byband_flux

    # integer :: icol, ilev, igpt, ibnd
    #$acc parallel loop collapse(3) copyin(spectral_flux, band_lims) copyout(byband_flux)
    byband_flux = Array{FT}(undef, ncol, nlev, nbnd)
    for ibnd = 1:nbnd
      for ilev = 1:nlev
        for icol = 1:ncol
          byband_flux[icol, ilev, ibnd] =  spectral_flux[icol, ilev, band_lims[1, ibnd]]
          for igpt = band_lims[1,ibnd]+1: band_lims[2,ibnd]
            byband_flux[icol, ilev, ibnd] = byband_flux[icol, ilev, ibnd] + spectral_flux[icol, ilev, igpt]
          end
        end
      end
    end
    return byband_flux
  end
  # ----------------------------------------------------------------------------
  #
  # Net flux: Spectral reduction over all points
  #
  function net_byband(ncol, nlev, ngpt, nbnd, band_lims, spectral_flux_dn, spectral_flux_up) # bind (C)
    # integer,                               intent(in ) :: ncol, nlev, ngpt, nbnd
    # integer,  dimension(2,          nbnd), intent(in ) :: band_lims
    # real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux_dn, spectral_flux_up
    # real(FT), dimension(ncol, nlev, nbnd), intent(out) :: byband_flux_net

    # integer :: icol, ilev, igpt, ibnd
    byband_flux_net = Array{FT}(undef, ncol, nlev, nbnd)

    #$acc parallel loop collapse(3) copyin(spectral_flux_dn, spectral_flux_up, band_lims) copyout(byband_flux_net)
    for ibnd = 1:nbnd
      for ilev = 1:nlev
        for icol = 1:ncol
          igpt = band_lims[1,ibnd]
          byband_flux_net[icol, ilev, ibnd] = spectral_flux_dn[icol, ilev, igpt] -
                                              spectral_flux_up[icol, ilev, igpt]
          for igpt = band_lims[1,ibnd]+1:band_lims[2,ibnd]
            byband_flux_net[icol, ilev, ibnd] = byband_flux_net[icol, ilev, ibnd] +
                                                spectral_flux_dn[icol, ilev, igpt] -
                                                spectral_flux_up[icol, ilev, igpt]
          end
        end
      end
    end
    return byband_flux_net
  end
  # ----------------------------------------------------------------------------
  function net_byband(ncol, nlev, nbnd, byband_flux_dn, byband_flux_up) # bind (C)
    # integer,                               intent(in ) :: ncol, nlev, nbnd
    # real(FT), dimension(ncol, nlev, nbnd), intent(in ) :: byband_flux_dn, byband_flux_up
    # real(FT), dimension(ncol, nlev, nbnd), intent(out) :: byband_flux_net

    return byband_flux_dn .- byband_flux_up
  end

end # module
