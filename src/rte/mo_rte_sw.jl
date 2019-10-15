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
#  Contains a single routine to compute direct and diffuse fluxes of solar radiation given
#    atmospheric optical properties on a spectral grid
#    information about vertical ordering
#    boundary conditions
#      solar zenith angle, spectrally-resolved incident colimated flux, surface albedos for direct and diffuse radiation
#    optionally, a boundary condition for incident diffuse radiation
#
# It is the user's responsibility to ensure that boundary conditions (incident fluxes, surface albedos) are on the same
#   spectral grid as the optical properties.
#
# Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
#   whatever summary the user needs.
#
# The routine does error checking and choses which lower-level kernel to invoke based on
#   what kinds of optical properties are supplied
#
# -------------------------------------------------------------------------------------------------
module mo_rte_sw
#  use mo_rte_kind,      only: FT, wl
#  use mo_util_array,    only: any_vals_less_than, any_vals_outside
#  use mo_optical_props, only: ty_optical_props, &
#                              ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
#  use mo_fluxes,        only: ty_fluxes
#  use mo_rte_solver_kernels, &
#                        only: apply_BC, sw_solver_noscat, sw_solver_2stream
#  implicit none
#  private

#  public :: rte_sw

  using ..mo_util_array
  using ..mo_optical_props
  using ..mo_fluxes
  using ..mo_rte_solver_kernels
  using ..fortran_intrinsics

  export rte_sw!, expand_and_transpose

  # --------------------------------------------------
  function rte_sw!(atmos::ty_optical_props_arry,
                   top_at_1,
                   mu0::Array{FT},
                   inc_flux,
                   sfc_alb_dir,
                   sfc_alb_dif,
                   fluxes,
                   inc_flux_dif=nothing) where {FT<:AbstractFloat} #result(error_msg)
#    class(ty_optical_props_arry), intent(in   ) :: atmos           # Optical properties provided as arrays
#    logical,                      intent(in   ) :: top_at_1        # Is the top of the domain at index 1?
#                                                                   # (if not, ordering is bottom-to-top)
#    real(FT), dimension(:),       intent(in   ) :: mu0             # cosine of solar zenith angle (ncol)
#    real(FT), dimension(:,:),     intent(in   ) :: inc_flux,    &  # incident flux at top of domain [W/m2] (ncol, ngpt)
#                                                   sfc_alb_dir, &  # surface albedo for direct and
#                                                   sfc_alb_dif     # diffuse radiation (nband, ncol)
#    class(ty_fluxes),             intent(inout) :: fluxes          # Class describing output calculations
#    real(FT), dimension(:,:), optional, &
#                                  intent(in   ) :: inc_flux_dif    # incident diffuse flux at top of domain [W/m2] (ncol, ngpt)
#    character(len=128)                          :: error_msg       # If empty, calculation was successful
#    # --------------------------------
#    #
#    # Local variables
#    #
#    integer :: ncol, nlay, ngpt, nband
#    integer :: icol

#    real(FT), dimension(:,:,:), allocatable :: gpt_flux_up, gpt_flux_dn, gpt_flux_dir
#    real(FT), dimension(:,:),   allocatable :: sfc_alb_dir_gpt, sfc_alb_dif_gpt
    # ------------------------------------------------------------------------------------

    ncol  = get_ncol(atmos)
    nlay  = get_nlay(atmos)
    ngpt  = get_ngpt(atmos)
    nband = get_nband(atmos)

    # ------------------------------------------------------------------------------------
    #
    # Error checking -- consistency of sizes and validity of values
    #
    # --------------------------------
    if(!are_desired(fluxes))
      error("rte_sw: no space allocated for fluxes")
    end

    #
    # Sizes and values of input arrays
    #
    if( size(mu0)[1] !=  ncol )
      error("rte_sw: mu0 inconsistently sized")
    end

    if(any_vals_outside(mu0, FT(0), FT(1)))
      error("rte_sw: one or more mu0 <= 0 or > 1")
    end

    if(any([size(inc_flux,1), size(inc_flux,2)] .!= [ncol, ngpt]))
      error("rte_sw: inc_flux inconsistently sized")
    end

    if any_vals_less_than(inc_flux, FT(0))
      error("rte_sw: one or more inc_flux < 0")
    end

    if present(inc_flux_dif)
      if any([size(inc_flux_dif)[1], size(inc_flux_dif)[2]] .!= [ncol, ngpt])
        error("rte_sw: inc_flux_dif inconsistently sized")
      end
      if any_vals_less_than(inc_flux_dif, FT(0))
        error("rte_sw: one or more inc_flux_dif < 0")
      end
    end

    if any([size(sfc_alb_dir)[1], size(sfc_alb_dir)[2]] .!= [nband, ncol])
      error("rte_sw: sfc_alb_dir inconsistently sized")
    end
    if(any_vals_outside(sfc_alb_dir,  FT(0), FT(1)))
      error("rte_sw: sfc_alb_dir out of bounds [0,1]")
    end
    if(any([size(sfc_alb_dif)[1], size(sfc_alb_dif)[2]] .!= [nband, ncol]))
      error("rte_sw: sfc_alb_dif inconsistently sized")
    end
    if(any_vals_outside(sfc_alb_dif, FT(0), FT(1)))
      error("rte_sw: sfc_alb_dif out of bounds [0,1]")
    end

    # ------------------------------------------------------------------------------------
    gpt_flux_up  = Array{FT}(undef,ncol, nlay+1, ngpt)
    gpt_flux_dn  = Array{FT}(undef,ncol, nlay+1, ngpt)
    gpt_flux_dir = Array{FT}(undef,ncol, nlay+1, ngpt)

    # ------------------------------------------------------------------------------------
    # Lower boundary condition -- expand surface albedos by band to gpoints
    #   and switch dimension ordering
#    #$acc enter data create(sfc_alb_dir_gpt, sfc_alb_dif_gpt)
    sfc_alb_dir_gpt = expand_and_transpose(atmos, sfc_alb_dir)
    sfc_alb_dif_gpt = expand_and_transpose(atmos, sfc_alb_dif)

    # test_data(sfc_alb_dir_gpt, "sfc_alb_dir_gpt")
    # test_data(sfc_alb_dif_gpt, "sfc_alb_dif_gpt")

    # ------------------------------------------------------------------------------------
    #
    # Compute the radiative transfer...
    #
    #
    # Apply boundary conditions
    #   On input flux_dn is the diffuse component; the last action in each solver is to add
    #   direct and diffuse to represent the total, consistent with the LW
    #
#    #$acc enter data copyin(mu0)
#    #$acc enter data create(gpt_flux_up, gpt_flux_dn, gpt_flux_dir)

#    #$acc enter data copyin(inc_flux)
    gpt_flux_dir = apply_BC(ncol, nlay, ngpt, top_at_1,   inc_flux, mu0)
    # test_data(gpt_flux_dir, "BC_gpt_flux_dir")
#    #$acc exit data delete(inc_flux)
    if present(inc_flux_dif)
#      #$acc enter data copyin(inc_flux_dif)
      gpt_flux_dn  = apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dif)
#      #$acc exit data delete(inc_flux_dif)
    else
      gpt_flux_dn  = apply_BC(ncol, nlay, ngpt, top_at_1, FT)
    end
    # test_data(gpt_flux_dn, "BC_gpt_flux_dn")
#    select type (atmos)
      if ( atmos isa ty_optical_props_1scl )
        #
        # Direct beam only
        #
#        #$acc enter data copyin(atmos, atmos%tau)
        validate!(atmos)
        sw_solver_noscat!(ncol, nlay, ngpt, top_at_1,
                              atmos.tau, mu0,
                              gpt_flux_dir)
        #
        # No diffuse flux
        #
#        #gpt_flux_up = 0._wp
#        #gpt_flux_dn = 0._wp
#        #$acc exit data delete(atmos%tau, atmos)
      elseif ( atmos isa ty_optical_props_2str )
        #
        # two-stream calculation with scattering
        #
#        #$acc enter data copyin(atmos, atmos%tau, atmos%ssa, atmos%g)
        validate!(atmos)

        # test_data(atmos.tau, "atmos_tau")
        # test_data(atmos.ssa, "atmos_ssa")
        # test_data(atmos.g, "atmos_g")
        # test_data(gpt_flux_dn, "gpt_flux_dn_in")
        # test_data(gpt_flux_dir, "gpt_flux_dir_in")

        sw_solver_2stream!(ncol, nlay, ngpt, top_at_1,
                               atmos.tau, atmos.ssa, atmos.g, mu0,
                               sfc_alb_dir_gpt, sfc_alb_dif_gpt,
                               gpt_flux_up, gpt_flux_dn, gpt_flux_dir)

        # test_data(gpt_flux_up, "gpt_flux_up")
        # test_data(gpt_flux_dn, "gpt_flux_dn_out")
        # test_data(gpt_flux_dir, "gpt_flux_dir_out")

#        #$acc exit data delete(atmos%tau, atmos%ssa, atmos%g, atmos)
#        #$acc exit data delete(sfc_alb_dir_gpt, sfc_alb_dif_gpt)
      else                #class is (ty_optical_props_nstr)
        #
        # n-stream calculation
        #
        # not yet implemented so fail
        #
        error("sw_solver(...ty_optical_props_nstr...) not yet implemented")
      end

    #
    # ...and reduce spectral fluxes to desired output quantities
    #
    reduce!(fluxes,gpt_flux_up, gpt_flux_dn, atmos, top_at_1, gpt_flux_dir)
#    #$acc exit data delete(mu0)
#    #$acc exit data delete(gpt_flux_up, gpt_flux_dn, gpt_flux_dir)
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
  #
  function expand_and_transpose(ops::ty_optical_props,arr_in::Array{FT}) where FT
#    class(ty_optical_props),  intent(in ) :: ops
#    real(FT), dimension(:,:), intent(in ) :: arr_in  # (nband, ncol)
#    real(FT), dimension(:,:), intent(out) :: arr_out # (ncol, igpt)
#    # -------------
#    integer :: ncol, nband, ngpt
#    integer :: icol, iband, igpt
#    integer, dimension(2,ops%get_nband()) :: limits

    ncol  = size(arr_in,2)
    nband = get_nband(ops)
    ngpt  = get_ngpt(ops)
    arr_out = Array{FT}(undef,ncol, ngpt)
    limits = get_band_lims_gpoint(ops)
#    #$acc parallel loop collapse(2) copyin(arr_in, limits)
    for iband = 1:nband
      for icol = 1:ncol
        for igpt = limits[1, iband]:limits[2, iband]
          arr_out[icol, igpt] = arr_in[iband,icol]
        end
      end
    end
    return arr_out
  end
#-------------------------------------------------------------
end
