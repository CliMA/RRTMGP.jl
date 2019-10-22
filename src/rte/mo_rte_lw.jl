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
#    atmospheric optical properties, spectrally-resolved
#    information about vertical ordering
#    internal Planck source functions, defined per g-point on the same spectral grid at the atmosphere
#    boundary conditions: surface emissivity defined per band
#    optionally, a boundary condition for incident diffuse radiation
#    optionally, an integer number of angles at which to do Gaussian quadrature if scattering is neglected
#
# If optical properties are supplied via class ty_optical_props_1scl (absorption optical thickenss only)
#    then an emission/absorption solver is called
#    If optical properties are supplied via class ty_optical_props_2str fluxes are computed via
#    two-stream calculations and adding.
#
# It is the user's responsibility to ensure that emissivity is on the same
#   spectral grid as the optical properties.
#
# Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
#   whatever summary the user needs.
#
# The routine does error checking and choses which lower-level kernel to invoke based on
#   what kinds of optical properties are supplied
#
# -------------------------------------------------------------------------------------------------
module mo_rte_lw
#=  use mo_rte_kind,      only: FT, wl
  use mo_util_array,    only: any_vals_less_than, any_vals_outside
  use mo_optical_props, only: ty_optical_props, &
                              ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_source_functions,   &
                        only: ty_source_func_lw
  use mo_fluxes,        only: ty_fluxes
  use mo_rte_solver_kernels, &
                        only: apply_BC, lw_solver_noscat_GaussQuad, lw_solver_2stream
  implicit none
  private

  public :: rte_lw
contains =#


  using ..mo_util_array
  using ..mo_optical_props
  using ..mo_fluxes
  using ..mo_rte_solver_kernels
  using ..fortran_intrinsics

  export rte_lw!, expand_and_transpose
#  export rte_lw!, expand_and_transpose
  # --------------------------------------------------
  #
  # Interface using only optical properties and source functions as inputs; fluxes as outputs.
  #
  # --------------------------------------------------
#  function rte_lw(optical_props, top_at_1, &
#                  sources, sfc_emis,       &
#                  fluxes,                  &
#                  inc_flux, n_gauss_angles) result(error_msg)
  function rte_lw!(optical_props, top_at_1, 
                  sources, sfc_emis,      
                  fluxes,                  
                  inc_flux=nothing, n_gauss_angles=nothing) #result(error_msg)
#    class(ty_optical_props_arry), intent(in   ) :: optical_props     # Array of ty_optical_props. This type is abstract
#                                                                     # and needs to be made concrete, either as an array
#                                                                     # (class ty_optical_props_arry) or in some user-defined way
#    logical,                      intent(in   ) :: top_at_1          # Is the top of the domain at index 1?
#                                                                     # (if not, ordering is bottom-to-top)
#    type(ty_source_func_lw),      intent(in   ) :: sources
#    real(FT), dimension(:,:),     intent(in   ) :: sfc_emis    # emissivity at surface [] (nband, ncol)
#    class(ty_fluxes),             intent(inout) :: fluxes      # Array of ty_fluxes. Default computes broadband fluxes at all levels
#                                                               #   if output arrays are defined. Can be extended per user desires.
#    real(FT), dimension(:,:),   &
#              target, optional, intent(in   ) :: inc_flux    # incident flux at domain top [W/m2] (ncol, ngpts)
#    integer,          optional, intent(in   ) :: n_gauss_angles # Number of angles used in Gaussian quadrature
#                                                                # (no-scattering solution)
#    character(len=128)                        :: error_msg   # If empty, calculation was successful
    # --------------------------------
    #
    # Local variables
    #
    #integer :: ncol, nlay, ngpt, nband
    #integer :: n_quad_angs
    #integer :: icol, iband, igpt
    #real(FT), dimension(:,:,:), allocatable :: gpt_flux_up, gpt_flux_dn
    #real(FT), dimension(:,:),   allocatable :: sfc_emis_gpt
    # --------------------------------------------------
    #
    # Weights and angle secants for first order (k=1) Gaussian quadrature.
    #   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
    #   after Abramowitz & Stegun 1972, page 921
    #
    FT = eltype(optical_props.band_lims_wvn) 

    max_gauss_pts = Integer.(4)
        gauss_Ds  = reshape([1.66000000, 0.00000000, 0.00000000, 0.00000000,   # Diffusivity angle, not Gaussian angle
                             1.18350343, 2.81649655, 0.00000000, 0.00000000, 
                             1.09719858, 1.69338507, 4.70941630, 0.00000000, 
                             1.06056257, 1.38282560, 2.40148179, 7.15513024], 
                             max_gauss_pts, max_gauss_pts)

        gauss_wts = reshape([0.5000000000, 0.0000000000, 0.0000000000, 0.0000000000,
                             0.3180413817, 0.1819586183, 0.0000000000, 0.0000000000,
                             0.2009319137, 0.2292411064, 0.0698269799, 0.0000000000,
                             0.1355069134, 0.2034645680, 0.1298475476, 0.0311809710],
                             max_gauss_pts, max_gauss_pts)
    # ------------------------------------------------------------------------------------
    #
    # Error checking
    #   if inc_flux is present it has the right dimensions, is positive definite
    #
    # --------------------------------
    ncol  = get_ncol(optical_props)
    nlay  = get_nlay(optical_props)
    ngpt  = get_ngpt(optical_props)
    nband = get_nband(optical_props)
    error_msg = ""

    # ------------------------------------------------------------------------------------
    #
    # Error checking -- consistency of sizes and validity of values
    #
    # --------------------------------
    if(!are_desired(fluxes))
      error("rte_sw: no space allocated for fluxes")
    end
    #
    # Source functions
    #
    if any( [get_ncol(sources), get_nlay(sources), get_ngpt(sources)]  ≠ [ncol, nlay, ngpt]) 
      error("rte_lw: sources and optical properties inconsistently sized")
    end
    # Also need to validate

    #
    # Surface emissivity
    #
    if any( [size(sfc_emis,1), size(sfc_emis,2)] ≠ [nband, ncol])
      error("rte_lw: sfc_emis inconsistently sized")
    end

    if any_vals_outside(sfc_emis, FT.(0.), FT.(1.0)) 
      error("rte_lw: sfc_emis has values < 0 or > 1")
    end

    #
    # Incident flux, if present
    #
    if !(inc_flux isa Nothing)
      if any( [size(inc_flux,1), size(inc_flux,2)] ≠ [ncol, ngpt] )
        error("rte_lw: inc_flux inconsistently sized")
      end
      if any_vals_less_than(inc_flux, FT(0.0))
        error("rte_lw: inc_flux has values < 0")
      end
    end
    #
    # Number of quadrature points for no-scattering calculation
    #
    n_quad_angs = 1
    if !(n_gauss_angles isa Nothing)
      if (n_gauss_angles > max_gauss_pts)
        error("rte_lw: asking for too many quadrature points for no-scattering calculation")
      end
      if(n_gauss_angles < 1) 
        error("rte_lw: have to ask for at least one quadrature point for no-scattering calculation")
      end
      n_quad_angs = n_gauss_angles
    end
    #
    # Ensure values of tau, ssa, and g are reasonable
    #
    validate!(optical_props)
    # ------------------------------------------------------------------------------------
    #
    #    Lower boundary condition -- expand surface emissivity by band to gpoints
    #
#    allocate(gpt_flux_up (ncol, nlay+1, ngpt), gpt_flux_dn(ncol, nlay+1, ngpt))
#    allocate(sfc_emis_gpt(ncol,         ngpt))
  
    gpt_flux_up  = Array{FT}(undef,ncol,nlay+1,ngpt)
    gpt_flux_dn  = Array{FT}(undef,ncol,nlay+1,ngpt)
    sfc_emis_gpt = Array{FT}(undef,ncol,       ngpt)

#    ##$acc enter data copyin(sources, sources%lay_source, sources%lev_source_inc, sources%lev_source_dec, sources%sfc_source)
#    #$acc enter data copyin(optical_props)
#    #$acc enter data create(gpt_flux_dn, gpt_flux_up)
#    #$acc enter data create(sfc_emis_gpt)
    sfc_emis_gpt = expand_and_transpose(optical_props, sfc_emis)
    #
    #   Upper boundary condition
    #
    if !(inc_flux isa Nothing)
      #$acc enter data copyin(inc_flux)
      gpt_flux_dn = apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux)
      #$acc exit data delete(inc_flux)
    else
      #
      # Default is zero incident diffuse flux
      #
      gpt_flux_dn = apply_BC(ncol, nlay, ngpt, top_at_1, FT)
    end

    #
    # Compute the radiative transfer...
    #

    if optical_props isa ty_optical_props_1scl

      # No scattering two-stream calculation
      lw_solver_noscat_GaussQuad!(ncol, nlay, ngpt, top_at_1, 
                              n_quad_angs, gauss_Ds[1:n_quad_angs,n_quad_angs], gauss_wts[1:n_quad_angs,n_quad_angs], 
                              optical_props.tau,                                                  
                              sources.lay_source, sources.lev_source_inc, sources.lev_source_dec, 
                              sfc_emis_gpt, sources.sfc_source,  
                              gpt_flux_up, gpt_flux_dn)

    elseif optical_props isa ty_optical_props_2str
      # two-stream calculation with scattering
#        call lw_solver_2stream(ncol, nlay, ngpt, logical(top_at_1, wl), 
#                               optical_props%tau, optical_props%ssa, optical_props%g,              
#                               sources%lay_source, sources%lev_source_inc, sources%lev_source_dec, 
#                               sfc_emis_gpt, sources%sfc_source,       
#                               gpt_flux_up, gpt_flux_dn)
      lw_solver_2stream!(ncol, nlay, ngpt, top_at_1,
                                 optical_props.tau, optical_props.ssa, optical_props.g,
                                 sources.lay_source, sources.lev_source_inc, sources.lev_source_dec, sfc_emis_gpt, sources.sfc_src,
                                 gpt_flux_up, gpt_flux_dn)
    elseif optical_props isa ty_optical_props_nstr
      # n-stream calculation
      error("lw_solver(...ty_optical_props_nstr...) not yet implemented")
    end
#    error_msg = fluxes%reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1)

    reduce!(fluxes,gpt_flux_up, gpt_flux_dn, optical_props, top_at_1)
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
  #
  #--------------------------------------------------------------------------------------------------------------------
end #---Module
