# This code is part of
# RRTM for GCM Applications - Parallel (RRTMGP)
#
# Eli Mlawer and Robert Pincus
# Andre Wehe and Jennifer Delamere
# email:  rrtmgp@aer.com
#
# Copyright 2015,  Atmospheric and Environmental Research and
# Regents of the University of Colorado.  All right reserved.
#
# Use and duplication is permitted under the terms of the
#    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
#
#
# This module is for packaging output quantities from RRTMGP based on spectral flux profiles
#    This implementation provides band-by-band flux profiles
#
module mo_fluxes_byband
  # use mo_rte_kind,      only: FT
  # use mo_fluxes,        only: ty_fluxes, ty_fluxes_broadband
  using ..mo_fluxes
  # use mo_optical_props, only: ty_optical_props
  using ..mo_optical_props
  # use mo_fluxes_byband_kernels,
  #                       only: sum_byband, net_byband
  using ..mo_fluxes_byband_kernels
  # implicit none

  # Output from radiation calculations
  #   Data components are pointers so results can be written directly into memory
  #   reduce() function accepts spectral flux profiles

mutable struct ty_fluxes_byband{FT} <: ty_fluxes{FT}
  fluxes_broadband
  # flux_up#::Array{FT,2}
  # flux_dn#::Array{FT,2}
  # flux_net#::Array{FT,2}
  # flux_dn_dir#::Array{FT,2}
  bnd_flux_up#::Array{FT,3}
  bnd_flux_dn#::Array{FT,3}
  bnd_flux_net#::Array{FT,3}
  bnd_flux_dn_dir#::Array{FT,3}
end

#   type, extends(ty_fluxes_broadband) :: ty_fluxes_byband
#     real(FT), dimension(:,:,:), pointer :: bnd_flux_up => NULL(),  # Band-by-band fluxes
#                                            bnd_flux_dn => NULL()    # (ncol, nlev, nband)
#     real(FT), dimension(:,:,:), pointer :: bnd_flux_net => NULL()   # Net (down - up)
#     real(FT), dimension(:,:,:), pointer :: bnd_flux_dn_dir => NULL() # Direct flux down
#   contains
#     procedure :: reduce => reduce_byband
#     procedure :: are_desired => are_desired_byband
#   end type ty_fluxes_byband
# contains
  # --------------------------------------------------------------------------------------
  function reduce_byband(this::ty_fluxes_byband, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir=nothing)# result(error_msg)
    # class(ty_fluxes_byband),           intent(inout) :: this
    # real(kind=FT), dimension(:,:,:),   intent(in   ) :: gpt_flux_up # Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    # real(kind=FT), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn # Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    # class(ty_optical_props),           intent(in   ) :: spectral_disc  #< derived type with spectral information
    # logical,                           intent(in   ) :: top_at_1
    # real(kind=FT), dimension(:,:,:), optional,
    #                                    intent(in   ) :: gpt_flux_dn_dir# Direct flux down
    # character(len=128)                               :: error_msg
    # # ------
    # integer :: ncol, nlev, ngpt, nbnd
    # integer, dimension(2, spectral_disc%get_nband()) :: band_lims
    # # ------
    ncol = size(gpt_flux_up, 1)
    nlev = size(gpt_flux_up, 2)
    ngpt = get_ngpt(spectral_disc)
    nbnd = get_nband(spectral_disc)
    band_lims = deepcopy(get_band_lims_gpoint(spectral_disc))

    # Compute broadband fluxes
    #   This also checks that input arrays are consistently sized
    #
    reduce_broadband!(this.fluxes_broadband, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir)

    if(size(gpt_flux_up, 3) ≠ ngpt)
      error("reduce: spectral discretization and g-point flux arrays have differing number of g-points")
    end

    # Check sizes of output arrays
    if(associated(this.bnd_flux_up))
      if(any([size(this.bnd_flux_up, 1) ≠ ncol,
              size(this.bnd_flux_up, 2) ≠ nlev,
              size(this.bnd_flux_up, 3) ≠ nbnd]))
        error("reduce: bnd_flux_up array incorrectly sized (can't compute net flux either)")
      end
    end
    if(associated(this.bnd_flux_dn))
      if(any([size(this.bnd_flux_dn, 1) ≠ ncol,
              size(this.bnd_flux_dn, 2) ≠ nlev,
              size(this.bnd_flux_dn, 3) ≠ nbnd]))
        error("reduce: bnd_flux_dn array incorrectly sized (can't compute net flux either)")
      end
    end
    if(associated(this.bnd_flux_dn_dir))
      if(any([size(this.bnd_flux_dn_dir, 1) ≠ ncol,
              size(this.bnd_flux_dn_dir, 2) ≠ nlev,
              size(this.bnd_flux_dn_dir, 3) ≠ nbnd]))
        error("reduce: bnd_flux_dn_dir array incorrectly sized")
      end
    end
    if(associated(this.bnd_flux_net))
      if(any([size(this.bnd_flux_net, 1) ≠ ncol,
              size(this.bnd_flux_net, 2) ≠ nlev,
              size(this.bnd_flux_net, 3) ≠ nbnd]))
        error("reduce: bnd_flux_net array incorrectly sized (can't compute net flux either)")
      end
    end
    #
    # Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
    if(associated(this.bnd_flux_dn_dir) && !present(gpt_flux_dn_dir))
      error("reduce: requesting bnd_flux_dn_dir but direct flux hasn't been supplied")
    end

    # -------
    #$acc enter data copyin(band_lims)
    # Band-by-band fluxes
    # Up flux
    if(associated(this.bnd_flux_up))
      this.bnd_flux_up     = sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_up)
    end

    # -------
    # Down flux
    if(associated(this.bnd_flux_dn))
      this.bnd_flux_dn     = sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn)
    end

    if(associated(this.bnd_flux_dn_dir))
      this.bnd_flux_dn_dir = sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn_dir)
    end

    # -------
    # Net flux
    #
    if(associated(this.bnd_flux_net))
      #
      #  Reuse down and up results if possible
      #
      if(associated(this.bnd_flux_dn) && associated(this.bnd_flux_up))
        this.bnd_flux_net = net_byband(ncol, nlev,       nbnd,                             this.bnd_flux_dn, this.bnd_flux_up)
      else
        this.bnd_flux_net = net_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn, gpt_flux_up)
      end
    end
    #$acc exit data delete(band_lims)
  end
  # --------------------------------------------------------------------------------------
  # Are any fluxes desired from this set of g-point fluxes? We can tell because memory will
  #   be allocated for output
  #
  function are_desired_byband(this::ty_fluxes_byband)
    # class(ty_fluxes_byband), intent(in   ) :: this
    # logical                                :: are_desired_byband

    return any([associated(this.bnd_flux_up),
                associated(this.bnd_flux_dn),
                associated(this.bnd_flux_dn_dir),
                associated(this.bnd_flux_net),
                are_desired(this.fluxes_broadband)])
  end

end # module
