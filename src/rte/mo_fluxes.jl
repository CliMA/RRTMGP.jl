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
# Compute output quantities from RTE based on spectrally-resolved flux profiles
#    This module contains an abstract class and a broadband implmentation that sums over all spectral points
#    The abstract base class defines the routines that extenstions must implement: reduce() and are_desired()
#    The intent is for users to extend it as required, using mo_flxues_broadband as an example
#
# -------------------------------------------------------------------------------------------------
module mo_fluxes

using ..mo_optical_props
using ..mo_fluxes_broadband_kernels


  # -----------------------------------------------------------------------------------------------
  #
  # Abstract base class
  #   reduce() function accepts spectral flux profiles, computes desired outputs
  #   are_desired() returns a logical - does it makes sense to invoke reduce()?
  #
  # -----------------------------------------------------------------------------------------------

abstract type ty_fluxes{T} end

  # type, abstract, public :: ty_fluxes
  # contains
  #   procedure(reduce_abstract),      deferred, public :: reduce
  #   procedure(are_desired_abstract), deferred, public :: are_desired
  # end type ty_fluxes
  # -----------------------------------------------------------------------------------------------
  #
  # Class implementing broadband integration for the complete flux profile
  #   Data components are pointers so results can be written directly into memory
  #
  # -----------------------------------------------------------------------------------------------

struct ty_fluxes_broadband{T} <: ty_fluxes{T}
  flux_up::Array{T,2}
  flux_dn::Array{T,2}
  flux_net::Array{T,2}
  flux_dn_dir::Array{T,2}
end

  # type, extends(ty_fluxes), public :: ty_fluxes_broadband
  #   real(wp), dimension(:,:), pointer :: flux_up => NULL(), flux_dn => NULL()
  #   real(wp), dimension(:,:), pointer :: flux_net => NULL()    # Net (down - up)
  #   real(wp), dimension(:,:), pointer :: flux_dn_dir => NULL() # Direct flux down
  # contains
  #   procedure, public :: reduce      => reduce_broadband
  #   procedure, public :: are_desired => are_desired_broadband
  # end type ty_fluxes_broadband
  # -----------------------------------------------------------------------------------------------

  # --------------------------------------------------------------------------------------
  #
  # Broadband fluxes -- simply sum over the spectral dimension and report the whole profile
  #
  # --------------------------------------------------------------------------------------
"""
    class(ty_fluxes_broadband),        intent(inout) :: this
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up # Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn # Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_optical_props),           intent(in   ) :: spectral_disc  #< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional,
                                       intent(in   ) :: gpt_flux_dn_dir# Direct flux down
    character(len=128)                               :: error_msg
    # ------
    integer :: ncol, nlev, ngpt
    # ------
"""
  function reduce_broadband(this::ty_fluxes_broadband,
                            gpt_flux_up::A,
                            gpt_flux_dn::A,
                            spectral_disc::ty_optical_props,
                            top_at_1::Bool,
                            gpt_flux_dn_dir=nothing) where A

    ncol = size(gpt_flux_up, dim=1)
    nlev = size(gpt_flux_up, dim=2)
    ngpt = size(gpt_flux_up, dim=3)
    error_msg = ""

    #
    # Check array sizes
    #  Input arrays
    #
    if any([size(gpt_flux_dn, 1) ≠ ncol, size(gpt_flux_dn, 2) ≠ nlev, size(gpt_flux_dn, 3) ≠ ngpt])
      error_msg = "reduce: gpt_flux_dn array incorrectly sized"
      return error_msg
    end
    if gpt_flux_dn_dir ≠ nothing
      if any([size(gpt_flux_dn_dir, 1) ≠ ncol, size(gpt_flux_dn_dir, 2) ≠ nlev, size(gpt_flux_dn_dir, 3) ≠ ngpt])
        error_msg = "reduce: gpt_flux_dn_dir array incorrectly sized"
        return error_msg
      end
    end
    #
    # Output arrays
    #
    if associated(this.flux_up)
      if any([size(this.flux_up, 1), size(this.flux_up, 2)] ≠ [ncol,nlev])
        error_msg = "reduce: flux_up array incorrectly sized"
        return error_msg
      end
    end
    if associated(this.flux_dn)
      if any([size(this.flux_dn, 1), size(this.flux_dn, 2)] ≠ [ncol,nlev])
        error_msg = "reduce: flux_dn array incorrectly sized"
        return error_msg
      end
    end
    if associated(this.flux_net)
      if any([size(this.flux_net, 1), size(this.flux_net, 2)] ≠ [ncol,nlev])
        error_msg = "reduce: flux_net array incorrectly sized"
        return error_msg
      end
    end
    if associated(this.flux_dn_dir)
      if any([size(this.flux_dn_dir, 1), size(this.flux_dn_dir, 2)] ≠ [ncol,nlev])
        error_msg = "reduce: flux_dn_dir array incorrectly sized"
        return error_msg
      end
    end
    #
    # Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
    if associated(this.flux_dn_dir) && !present(gpt_flux_dn_dir)
      error_msg = "reduce: requesting direct downward flux but this hasn't been supplied"
      return
    end

    #
    # Broadband fluxes - call the kernels
    #
    if associated(this.flux_up    )
      sum_broadband!(ncol, nlev, ngpt, gpt_flux_up,     this.flux_up)
    end
    if associated(this.flux_dn    )
      sum_broadband!(ncol, nlev, ngpt, gpt_flux_dn,     this.flux_dn)
    end
    if associated(this.flux_dn_dir)
      sum_broadband!(ncol, nlev, ngpt, gpt_flux_dn_dir, this.flux_dn_dir)
    end

    if associated(this.flux_net)
      #
      #  Reuse down and up results if possible
      #
      if associated(this.flux_dn) && associated(this.flux_up)
        net_broadband!(ncol, nlev,      this.flux_dn, this.flux_up, this.flux_net)
      else
        net_broadband!(ncol, nlev, ngpt, gpt_flux_dn,  gpt_flux_up, this.flux_net)
      end
    end
  end # function
  # --------------------------------------------------------------------------------------
  #
  # Are any fluxes desired from this set of g-point fluxes? We can tell because memory will
  #   be allocated for output
  #
  # --------------------------------------------------------------------------------------
"""
    class(ty_fluxes_broadband), intent(in   ) :: this
    logical                                   :: are_desired_broadband
"""
  function are_desired_broadband(this::ty_fluxes_broadband)
    return any( [associated(this.flux_up), associated(this.flux_dn), associated(this.flux_dn_dir), associated(this.flux_net)] )
  end
  # --------------------------------------------------------------------------------------
end # module
