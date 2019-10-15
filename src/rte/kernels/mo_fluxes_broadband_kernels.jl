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
#
# Kernels for computing broadband fluxes by summing over all elements in the spectral dimension
#
# -------------------------------------------------------------------------------------------------



module mo_fluxes_broadband_kernels

export sum_broadband, net_broadband

#  use, intrinsic :: iso_c_binding
#  use mo_rte_kind, only: FT
#  implicit none
#  private
#  public :: sum_broadband, net_broadband

#  interface net_broadband
#    module procedure net_broadband_full, net_broadband_precalc
#  end interface net_broadband
#contains



  # ----------------------------------------------------------------------------
    #
    # Spectral reduction over all points
    #
  function sum_broadband(ncol, nlev, ngpt, spectral_flux::Array{FT}) where FT #bind(C, name="sum_broadband")
#    integer,                               intent(in ) :: ncol, nlev, ngpt
#    real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux
#    real(FT), dimension(ncol, nlev),       intent(out) :: broadband_flux

#    integer :: icol, ilev, igpt

#    #$acc enter data copyin(spectral_flux) create(broadband_flux)
#    #$acc parallel loop collapse(2)

    broadband_flux = Array{FT}(undef, ncol, nlev)

    for ilev = 1:nlev
      for icol = 1:ncol
        broadband_flux[icol, ilev] =  spectral_flux[icol, ilev, 1]
      end
    end

#    #$acc parallel loop collapse(3)
    for igpt = 2:ngpt
      for ilev = 1:nlev
        for icol = 1:ncol
          #$acc atomic update
          broadband_flux[icol, ilev] = broadband_flux[icol, ilev] + spectral_flux[icol, ilev, igpt]
        end
      end
    end
    return broadband_flux

#    #$acc exit data delete(spectral_flux) copyout(broadband_flux)

  end




  # ----------------------------------------------------------------------------
  #
  # Net flux: Spectral reduction over all points
  #
  function net_broadband(ncol, nlev, ngpt, spectral_flux_dn::Array{FT}, spectral_flux_up) where FT # &
#    bind(C, name="net_broadband_full")
#    integer,                               intent(in ) :: ncol, nlev, ngpt
#    real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux_dn, spectral_flux_up
#    real(FT), dimension(ncol, nlev),       intent(out) :: broadband_flux_net

#    integer  :: icol, ilev, igpt
#    real(FT) :: diff_

#    #$acc enter data copyin(spectral_flux_dn, spectral_flux_up) create(broadband_flux_net)
#    #$acc parallel loop collapse(2)
    broadband_flux_net = Array{FT}(undef, ncol, nlev)
    for ilev = 1:nlev
      for icol = 1:ncol
        diff_ = spectral_flux_dn[icol, ilev, 1   ] - spectral_flux_up[icol, ilev,     1]
        broadband_flux_net[icol, ilev] = diff_
      end
    end
#    #$acc parallel loop collapse(3)
    for igpt = 2:ngpt
      for ilev = 1:nlev
        for icol = 1:ncol
          diff_ = spectral_flux_dn[icol, ilev, igpt] - spectral_flux_up[icol, ilev, igpt]
#          #$acc atomic update
          broadband_flux_net[icol, ilev] = broadband_flux_net[icol, ilev] + diff_
        end
      end
    end
    return broadband_flux_net
#    #$acc exit data delete(spectral_flux_dn, spectral_flux_up) copyout(broadband_flux_net)
  end
  # ----------------------------------------------------------------------------
  #
  # Net flux when bradband flux up and down are already available
  #
  function net_broadband(ncol, nlev, flux_dn, flux_up) # &
#    bind(C, name="net_broadband_precalc")
#    integer,                         intent(in ) :: ncol, nlev
#    real(FT), dimension(ncol, nlev), intent(in ) :: flux_dn, flux_up
#    real(FT), dimension(ncol, nlev), intent(out) :: broadband_flux_net

#    integer  :: icol, ilev
#    #$acc enter data copyin(flux_dn, flux_up) create(broadband_flux_net)
#    #$acc parallel loop collapse(2)
    return flux_dn .- flux_up
    # for ilev = 1:nlev
    #   for icol = 1:ncol
    #      broadband_flux_net[icol,ilev] = flux_dn[icol,ilev] - flux_up[icol,ilev]
    #    end
    # end
#    #$acc exit data delete(flux_dn, flux_up) copyout(broadband_flux_net)
  end
  # ----------------------------------------------------------------------------
end
