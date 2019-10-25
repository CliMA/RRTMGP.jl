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
# Provides cloud optical properties as a function of effective radius for the RRTMGP bands
#   Based on Mie calculations for liquid
#     and results from doi:10.1175/JAS-D-12-039.1 for ice with variable surface roughness
#   Can use either look-up tables or Pade approximates according to which data has been loaded
#   Mike Iacono (AER) is the original author
#
# The class can be used as-is but is also intended as an example of how to extend the RTE framework
# -------------------------------------------------------------------------------------------------

module mo_cloud_optics

using ..mo_optical_props


type, extends(ty_optical_props), public :: ty_cloud_optics
  private
  #
  # Ice surface roughness category - needed for Yang (2013) ice optics parameterization
  #
  integer            :: icergh = 0  # (1 = none, 2 = medium, 3 = high)
  #
  # Lookup table information
  #
  # Upper and lower limits of the tables
  real(FT) :: radliq_lwr = FT(0), radliq_upr = FT(0)
  real(FT) :: radice_lwr = FT(0), radice_upr = FT(0)
  # How many steps in the table? (for convenience)
  integer  :: liq_nsteps = 0,        ice_nsteps = 0
  # How big is each step in the table?
  real(FT) :: liq_step_size = FT(0), ice_step_size = FT(0)
  #
  # The tables themselves.
  #
  real(FT), dimension(:,:    ), allocatable :: lut_extliq, lut_ssaliq, lut_asyliq # (nsize_liq, nbnd)
  real(FT), dimension(:,:,:  ), allocatable :: lut_extice, lut_ssaice, lut_asyice # (nsize_ice, nbnd, nrghice)

  #
  # Pade approximant coefficients
  #
  real(FT), dimension(:,:,:  ), allocatable :: pade_extliq                 # (nbnd, nsizereg, ncoeff_ext)
  real(FT), dimension(:,:,:  ), allocatable :: pade_ssaliq,  pade_asyliq   # (nbnd, nsizereg, ncoeff_ssa_g)
  real(FT), dimension(:,:,:,:), allocatable :: pade_extice                 # (nbnd, nsizereg, ncoeff_ext, nrghice)
  real(FT), dimension(:,:,:,:), allocatable :: pade_ssaice, pade_asyice    # (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  # Particle size regimes for Pade formulations
  real(FT), dimension(:), allocatable :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq  # (nbound)
  real(FT), dimension(:), allocatable :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice  # (nbound)
  # -----
# contains
#   generic,   public :: load  => load_lut, load_pade
#   procedure, public :: finalize
#   procedure, public :: cloud_optics
#   procedure, public :: get_min_radius_liq
#   procedure, public :: get_min_radius_ice
#   procedure, public :: get_max_radius_liq
#   procedure, public :: get_max_radius_ice
#   procedure, public :: get_num_ice_roughness_types
#   procedure, public :: set_ice_roughness
#   # Internal procedures
#   procedure, private :: load_lut
#   procedure, private :: load_pade
# end type ty_cloud_optics

# ------------------------------------------------------------------------------
#
# Routines to load data needed for cloud optics calculations. Two routines: one to load
#    lookup-tables and one for coefficients for Pade approximates
#
# ------------------------------------------------------------------------------
function load_lut(this::ty_cloud_optics, band_lims_wvn,
                  radliq_lwr, radliq_upr, radliq_fac,
                  radice_lwr, radice_upr, radice_fac,
                  lut_extliq, lut_ssaliq, lut_asyliq,
                  lut_extice, lut_ssaice, lut_asyice)
  # class(ty_cloud_optics),     intent(inout) :: this
  # real(FT), dimension(:,:),   intent(in   ) :: band_lims_wvn # Spectral discretization
  # # Lookup table interpolation constants
  # # Lower and upper bounds of the tables; also the constant for calculating interpolation indices for liquid
  # real(FT),                   intent(in   ) :: radliq_lwr, radliq_upr, radliq_fac
  # real(FT),                   intent(in   ) :: radice_lwr, radice_upr, radice_fac
  # # LUT coefficients
  # # Extinction, single-scattering albedo, and asymmetry parameter for liquid and ice respectively
  # real(FT), dimension(:,:),   intent(in)    :: lut_extliq, lut_ssaliq, lut_asyliq
  # real(FT), dimension(:,:,:), intent(in)    :: lut_extice, lut_ssaice, lut_asyice
  # character(len=128)    :: error_msg
  # # -------
  # #
  # # Local variables
  # #
  # integer               :: nbnd, nrghice, nsize_liq, nsize_ice

  init!(this, "RRTMGP cloud optics", band_lims_wvn)
  #
  # LUT coefficient dimensions
  #
  nsize_liq = size(lut_extliq,1)
  nsize_ice = size(lut_extice,1)
  nbnd      = size(lut_extliq,2)
  nrghice   = size(lut_extice,3)
  #
  # Error checking
  #   Can we check for consistency between table bounds and _fac?
  #
  if(nbnd ≠ get_nband(this))
    error("cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization")
  end
  if(size(lut_extice, 2) ≠ nbnd)
    error("cloud_optics%init(): array lut_extice has the wrong number of bands")
  end
  if(any([size(lut_ssaliq, 1), size(lut_ssaliq, 2)] ≠ [nsize_liq, nbnd]))
    error("cloud_optics%init(): array lut_ssaliq isn't consistently sized")
  end
  if(any([size(lut_asyliq, 1), size(lut_asyliq, 2)] ≠ [nsize_liq, nbnd]))
    error("cloud_optics%init(): array lut_asyliq isn't consistently sized")
  end
  if(any([size(lut_ssaice, 1), size(lut_ssaice, 2), size(lut_ssaice, 3)] ≠ [nsize_ice, nbnd, nrghice]))
    error("cloud_optics%init(): array lut_ssaice  isn't consistently sized")
  end
  if(any([size(lut_asyice, 1), size(lut_asyice, 2), size(lut_asyice, 3)] ≠ [nsize_ice, nbnd, nrghice]))
    error("cloud_optics%init(): array lut_asyice  isn't consistently sized")
  end

  this.liq_nsteps = nsize_liq
  this.ice_nsteps = nsize_ice
  this.liq_step_size = (radliq_upr - radliq_lwr)/FT(nsize_liq-1)
  this.ice_step_size = (radice_upr - radice_lwr)/FT(nsize_ice-1)
  # Allocate LUT coefficients
  this.lut_extliq = Array{FT}(undef, nsize_liq, nbnd)
  this.lut_ssaliq = Array{FT}(undef, nsize_liq, nbnd)
  this.lut_asyliq = Array{FT}(undef, nsize_liq, nbnd)
  this.lut_extice = Array{FT}(undef, nsize_ice, nbnd, nrghice)
  this.lut_ssaice = Array{FT}(undef, nsize_ice, nbnd, nrghice)
  this.lut_asyice = Array{FT}(undef, nsize_ice, nbnd, nrghice)

  # Load LUT constants
  this.radliq_lwr = radliq_lwr
  this.radliq_upr = radliq_upr
  this.radice_lwr = radice_lwr
  this.radice_upr = radice_upr

  # Load LUT coefficients
  this.lut_extliq .= lut_extliq
  this.lut_ssaliq .= lut_ssaliq
  this.lut_asyliq .= lut_asyliq
  this.lut_extice .= lut_extice
  this.lut_ssaice .= lut_ssaice
  this.lut_asyice .= lut_asyice
end
# ------------------------------------------------------------------------------
#
# Cloud optics initialization function - Pade
#
# ------------------------------------------------------------------------------
function load_pade(this, band_lims_wvn,
                   pade_extliq, pade_ssaliq, pade_asyliq,
                   pade_extice, pade_ssaice, pade_asyice,
                   pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq,
                   pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice)
  # class(ty_cloud_optics),       intent(inout) :: this          # cloud specification data
  # real(FT), dimension(:,:),     intent(in   ) :: band_lims_wvn # Spectral discretization
  # #
  # # Pade coefficients: extinction, single-scattering albedo, and asymmetry factor for liquid and ice
  # #
  # real(FT), dimension(:,:,:),   intent(in)    :: pade_extliq, pade_ssaliq, pade_asyliq
  # real(FT), dimension(:,:,:,:), intent(in)    :: pade_extice, pade_ssaice, pade_asyice
  # #
  # # Boundaries of size regimes. Liquid and ice are separate;
  # #   extinction is fit to different numbers of size bins than single-scattering albedo and asymmetry factor
  # #
  # real(FT),  dimension(:),       intent(in)    :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq
  # real(FT),  dimension(:),       intent(in)    :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice
  # character(len=128)    :: error_msg

  # integer               :: nbnd, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g, nbound

  # Pade coefficient dimensions
  nbnd         = size(pade_extliq, 1)
  nsizereg     = size(pade_extliq, 2)
  ncoeff_ext   = size(pade_extliq, 3)
  ncoeff_ssa_g = size(pade_ssaliq, 3)
  nrghice      = size(pade_extice, 4)
  nbound       = size(pade_sizreg_extliq)
  init!(this, "RRTMGP cloud optics", band_lims_wvn)
  #
  # Error checking
  #
  if(nbnd ≠ get_nband(this))
    error("cloud_optics%init(): number of bands inconsistent between lookup tables, spectral discretization")
  end
  if(any([size(pade_ssaliq, 1), size(pade_ssaliq, 2), size(pade_ssaliq, 3)] ≠ [nbnd, nsizereg, ncoeff_ssa_g]))
    error("cloud_optics%init(): array pade_ssaliq isn't consistently sized")
  end
  if(any([size(pade_asyliq, 1), size(pade_asyliq, 2), size(pade_asyliq, 3)] ≠ [nbnd, nsizereg, ncoeff_ssa_g]))
    error("cloud_optics%init(): array pade_asyliq isn't consistently sized")
  end
  if(any([size(pade_extice, 1), size(pade_extice, 2), size(pade_extice, 3), size(pade_extice, 4)] ≠ [nbnd, nsizereg, ncoeff_ext, nrghice]))
    error("cloud_optics%init(): array pade_extice isn't consistently sized")
  end
  if(any([size(pade_ssaice, 1), size(pade_ssaice, 2), size(pade_ssaice, 3), size(pade_ssaice, 4)] ≠ [nbnd, nsizereg, ncoeff_ssa_g, nrghice]))
    error("cloud_optics%init(): array pade_ssaice isn't consistently sized")
  end
  if(any([size(pade_asyice, 1), size(pade_asyice, 2), size(pade_asyice, 3), size(pade_asyice, 4)] ≠ [nbnd, nsizereg, ncoeff_ssa_g, nrghice]))
    error("cloud_optics%init(): array pade_asyice isn't consistently sized")
  end
  if(any([size(pade_sizreg_ssaliq), size(pade_sizreg_asyliq), size(pade_sizreg_extice), size(pade_sizreg_ssaice), size(pade_sizreg_asyice)] ≠ nbound))
    error("cloud_optics%init(): one or more Pade size regime arrays are inconsistently sized")
  end
  if(nsizereg ≠ 3)
    error("cloud_optics%init(): Expecting precisely three size regimes for Pade approximants")
  end

  #
  # Consistency among size regimes
  #
  this.radliq_lwr = pade_sizreg_extliq[1]
  this.radliq_upr = pade_sizreg_extliq[nbound]
  this.radice_lwr = pade_sizreg_extice[1]
  this.radice_upr = pade_sizreg_extice[nbound]

  if(any([pade_sizreg_ssaliq(1), pade_sizreg_asyliq(1)] < this%radliq_lwr))
    error_msg = "cloud_optics%init(): one or more Pade size regimes have inconsistent lowest values"
  if(any([pade_sizreg_ssaice(1), pade_sizreg_asyice(1)] < this%radice_lwr))
    error_msg = "cloud_optics%init(): one or more Pade size regimes have inconsistent lower values"

  if(any([pade_sizreg_ssaliq(nbound), pade_sizreg_asyliq(nbound)] > this%radliq_upr))
    error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radliq_upr"
  if(any([pade_sizreg_ssaice(nbound), pade_sizreg_asyice(nbound)] > this%radice_upr))
    error_msg = "cloud_optics%init(): one or more Pade size regimes have lowest value less than radice_upr"

  #
  # Allocate Pade coefficients
  #
  this.pade_extliq = Array{FT}(undef, nbnd, nsizereg, ncoeff_ext)
  this.pade_ssaliq = Array{FT}(undef, nbnd, nsizereg, ncoeff_ssa_g)
  this.pade_asyliq = Array{FT}(undef, nbnd, nsizereg, ncoeff_ssa_g)
  this.pade_extice = Array{FT}(undef, nbnd, nsizereg, ncoeff_ext,   nrghice)
  this.pade_ssaice = Array{FT}(undef, nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  this.pade_asyice = Array{FT}(undef, nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  #
  # Allocate Pade coefficient particle size regime boundaries
  #
  this.pade_sizreg_extliq = Array{FT}(undef, nbound)
  this.pade_sizreg_ssaliq = Array{FT}(undef, nbound)
  this.pade_sizreg_asyliq = Array{FT}(undef, nbound)
  this.pade_sizreg_extice = Array{FT}(undef, nbound)
  this.pade_sizreg_ssaice = Array{FT}(undef, nbound)
  this.pade_sizreg_asyice = Array{FT}(undef, nbound)
  #
  # Load data
  #
  this.pade_extliq .= pade_extliq
  this.pade_ssaliq .= pade_ssaliq
  this.pade_asyliq .= pade_asyliq
  this.pade_extice .= pade_extice
  this.pade_ssaice .= pade_ssaice
  this.pade_asyice .= pade_asyice
  this.pade_sizreg_extliq .= pade_sizreg_extliq
  this.pade_sizreg_ssaliq .= pade_sizreg_ssaliq
  this.pade_sizreg_asyliq .= pade_sizreg_asyliq
  this.pade_sizreg_extice .= pade_sizreg_extice
  this.pade_sizreg_ssaice .= pade_sizreg_ssaice
  this.pade_sizreg_asyice .= pade_sizreg_asyice
end

function finalize!(this::ty_cloud_optics{FT}) where FT
  # class(ty_cloud_optics), intent(inout) :: this

  this.radliq_lwr = FT(0)
  this.radliq_upr = FT(0)
  this.radice_lwr = FT(0)
  this.radice_upr = FT(0)

  # Lookup table cloud optics coefficients
  if allocated(this.lut_extliq)
    deallocate!(this.lut_extliq)
    deallocate!(this.lut_ssaliq)
    deallocate!(this.lut_asyliq)
    deallocate!(this.lut_extice)
    deallocate!(this.lut_ssaice)
    deallocate!(this.lut_asyice)
    this.liq_nsteps = 0
    this.ice_nsteps = 0
    this.liq_step_size = FT(0)
    this.ice_step_size = FT(0)
  end

  # Pade cloud optics coefficients
  if allocated(this.pade_extliq)
    deallocate!(this.pade_extliq)
    deallocate!(this.pade_ssaliq)
    deallocate!(this.pade_asyliq)
    deallocate!(this.pade_extice)
    deallocate!(this.pade_ssaice)
    deallocate!(this.pade_asyice)
    deallocate!(this.pade_sizreg_extliq)
    deallocate!(this.pade_sizreg_ssaliq)
    deallocate!(this.pade_sizreg_asyliq)
    deallocate!(this.pade_sizreg_extice)
    deallocate!(this.pade_sizreg_ssaice)
    deallocate!(this.pade_sizreg_asyice)
  end
end


#
# Derive cloud optical properties from provided cloud physical properties

#
# Compute single-scattering properties
#
function cloud_optics(this::ty_cloud_optics{FT},
                      ncol, nlay, nbnd, nrghice,
                      liqmsk, icemsk, clwp, ciwp, reliq, reice,
                      optical_props) where FT
  # class(ty_cloud_optics),
  #           intent(in   ) :: this
  # integer,  intent(in   ) :: ncol, nlay, nbnd
  # integer,  intent(in   ) :: nrghice              # number of ice roughness categories
  # logical,  intent(in   ) :: liqmsk(ncol,nlay),  # Cloud mask for liquid and ice clouds respectively
  #                            icemsk(ncol,nlay)
  # real(FT), intent(in   ) :: ciwp  (ncol,nlay),      # cloud ice water path
  #                            clwp  (ncol,nlay),      # cloud liquid water path
  #                            reice (ncol,nlay),       # cloud ice particle effective size (microns)
  #                            reliq (ncol,nlay)      # cloud liquid particle effective radius (microns)
  # class(ty_optical_props_arry),
  #           intent(inout) :: optical_props
  #                                            # Dimensions: (ncol,nlay,nbnd)

  # character(len=128)      :: error_msg
  # # ------- Local -------
  # type(ty_optical_props_2str) :: clouds_liq, clouds_ice
  # integer  :: nsizereg, ibnd, imom
  # # ----------------------------------------
  # #
  # # Error checking
  # #
  # # ----------------------------------------
  if !(allocated(this.lut_extliq) || allocated(this.pade_extliq))
    error("cloud optics: no data has been initialized")
  end

  if !bands_are_equal(this, optical_props)
    error("cloud optics: optical properties don't have the same band structure")
  end

  if get_nband(optical_props) ≠ get_ngpt(optical_props) 
    error("cloud optics: optical properties must be requested by band not g-points")
  end

  if this.icergh < 1 || this.icergh > get_num_ice_roughness_types(this)
     error("cloud optics: cloud ice surface roughness flag is out of bounds")
   end

  if any(liqmsk && (reliq < this.radliq_lwr || reliq > this.radliq_upr))
    error("cloud optics: liquid effective radius is out of bounds")
  end

  if any(icemsk && (reice < this.radice_lwr || reice > this.radice_upr))
    error("cloud optics: ice effective radius is out of bounds")
  end

  if any((liqmsk &&  clwp < FT(0)) || (icemsk &&  ciwp < FT(0)))
    error("cloud optics: negative clwp or ciwp where clouds are supposed to be")
  end

  # ----------------------------------------
  #
  # Compute cloud optical properties. Use lookup tables if available, Pade approximants if not.
  #
  alloc_2str!(clouds_liq, ncol,nlay, this)

  alloc_2str!(clouds_ice, ncol,nlay, this)
  if allocated(this.lut_extliq)
    #
    # Liquid
    #
    clouds_liq.tau = compute_from_table(ncol,nlay,nbnd,liqmsk,reliq,this.liq_nsteps,this.liq_step_size,this.radliq_lwr, this.lut_extliq)
    clouds_liq.ssa = compute_from_table(ncol,nlay,nbnd,liqmsk,reliq,this.liq_nsteps,this.liq_step_size,this.radliq_lwr, this.lut_ssaliq)
    clouds_liq.g   = compute_from_table(ncol,nlay,nbnd,liqmsk,reliq,this.liq_nsteps,this.liq_step_size,this.radliq_lwr, this.lut_asyliq)
    for ibnd = 1:nbnd
      for icol in 1:ncol
        for ilay in 1:nlay
          liqmsk[icol,ilay] && (clouds_liq.tau[icol,ilay,ibnd] = clouds_liq.tau[icol,ilay,ibnd] * clwp[icol,ilay])
        end
      end
    end
    #
    # Ice
    #
    clouds_ice%tau = compute_from_table(ncol,nlay,nbnd,icemsk,reice,this.ice_nsteps,this.ice_step_size,this.radice_lwr, this.lut_extice[:,:,this.icergh])
    clouds_ice%ssa = compute_from_table(ncol,nlay,nbnd,icemsk,reice,this.ice_nsteps,this.ice_step_size,this.radice_lwr, this.lut_ssaice[:,:,this.icergh])
    clouds_ice%g   = compute_from_table(ncol,nlay,nbnd,icemsk,reice,this.ice_nsteps,this.ice_step_size,this.radice_lwr, this.lut_asyice[:,:,this.icergh])
    for ibnd = 1:nbnd
      for icol in 1:ncol
        for ilay in 1:nlay
          icemsk[icol,ilay] && (clouds_ice.tau[icol,ilay,ibnd] = clouds_ice.tau[icol,ilay,ibnd] * ciwp[icol,ilay])
        end
      end
    end
  else
    #
    # Cloud optical properties from Pade coefficient method
    #
    # This assumes that all the Pade treaments have the same number of size regimes
    nsizereg = size(this.pade_extliq,2)
    #
    # Liquid
    #
    clouds_liq.tau =                    compute_from_pade(ncol,nlay,nbnd,liqmsk,reliq,nsizereg,this.pade_sizreg_extliq,2,3,this.pade_extliq)
    clouds_liq.ssa = FT(1) - max(FT(0), compute_from_pade(ncol,nlay,nbnd,liqmsk,reliq,nsizereg,this.pade_sizreg_ssaliq,2,2,this.pade_ssaliq))
    clouds_liq.g   =                    compute_from_pade(ncol,nlay,nbnd,liqmsk,reliq,nsizereg,this.pade_sizreg_asyliq,2,2,this.pade_asyliq)
    for ibnd = 1:nbnd
      for icol in 1:ncol
        for ilay in 1:nlay
          liqmsk[icol,ilay] && (clouds_liq.tau[icol,ilay,ibnd] = clouds_liq.tau[icol,ilay,ibnd] * clwp[icol,ilay])
        end
      end
    end
    #
    # Ice
    #
    clouds_ice.tau =                    compute_from_pade(ncol,nlay,nbnd,icemsk,reice,nsizereg,this.pade_sizreg_extice,2,3, this.pade_extice[:,:,:,this.icergh])
    clouds_ice.ssa = FT(1) - max(FT(0), compute_from_pade(ncol,nlay,nbnd,icemsk,reice,nsizereg,this.pade_sizreg_ssaice,2,2, this.pade_ssaice[:,:,:,this.icergh]))
    clouds_ice.g   =                    compute_from_pade(ncol,nlay,nbnd,icemsk,reice,nsizereg,this.pade_sizreg_asyice,2,2, this.pade_asyice[:,:,:,this.icergh])
    for ibnd = 1:nbnd
      for icol in 1:ncol
        for ilay in 1:nlay
          icemsk[icol,ilay] && (clouds_ice.tau[icol,ilay,ibnd] = clouds_ice.tau[icol,ilay,ibnd] * ciwp[icol,ilay])
        end
      end
    end
  end

  #
  # Combine liquid and ice contributions into total cloud optical properties
  #
 increment!(clouds_ice, clouds_liq)
 #
 # Copy total cloud properties onto outputs
 #
 if optical_props isa ty_optical_props_1scl
   optical_props.tau[1:ncol,1:nlay,1:nbnd] = clouds_liq.tau[1:ncol,1:nlay,1:nbnd] *
                                    (FT(1) - clouds_liq.ssa[1:ncol,1:nlay,1:nbnd])
 elseif optical_props isa ty_optical_props_2str
   optical_props.tau[1:ncol,1:nlay,1:nbnd] = clouds_liq.tau[1:ncol,1:nlay,1:nbnd]
   optical_props.ssa[1:ncol,1:nlay,1:nbnd] = clouds_liq.ssa[1:ncol,1:nlay,1:nbnd]
   optical_props.g[  1:ncol,1:nlay,1:nbnd] = clouds_liq.g[  1:ncol,1:nlay,1:nbnd]
 elseif optical_props isa ty_optical_props_nstr
   optical_props.tau[  1:ncol,1:nlay,1:nbnd] = clouds_liq.tau[1:ncol,1:nlay,1:nbnd]
   optical_props.ssa[  1:ncol,1:nlay,1:nbnd] = clouds_liq.ssa[1:ncol,1:nlay,1:nbnd]
   optical_props.p[  1,1:ncol,1:nlay,1:nbnd] = clouds_liq.g[  1:ncol,1:nlay,1:nbnd]
   for imom = 2:get_nmom(optical_props)
     optical_props.p[imom,1:ncol,1:nlay,1:nbnd] = clouds_liq.g[          1:ncol,1:nlay,1:nbnd] *
                                                  optical_props.p[imom-1,1:ncol,1:nlay,1:nbnd]
   end
 end

end

"""
    set_ice_roughness!

"""
function set_ice_roughness!(this::ty_cloud_optics, icergh)
  # class(ty_cloud_optics), intent(inout) :: this
  # integer,                intent(in   ) :: icergh
  # character(len=128)                    :: error_msg

  icergh < 1 && error("set_ice_roughness!: must be > 0")

  this.icergh = icergh
end

function get_num_ice_roughness_types(this::ty_cloud_optics)
  # class(ty_cloud_optics), intent(in   ) :: this
  # integer                               :: i

  i = 0
  allocated(this.pade_extice) && (i = size(this.pade_extice, 4))
  allocated(this.lut_extice ) && (i = size(this.lut_extice,  3))
  return i
end

get_min_radius_liq(this::ty_cloud_optics) = this.radliq_lwr

get_max_radius_liq(this::ty_cloud_optics) = this.radliq_upr

get_min_radius_ice(this::ty_cloud_optics) = this.radice_lwr

get_max_radius_ice(this::ty_cloud_optics) = this.radice_upr

#
# Ancillary functions
#
#
# Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
#   elements starting at "offset." The table's second dimension is band.
# Returns 0 where the mask is false.
# We could also try gather/scatter for efficiency
#
function compute_from_table(ncol, nlay, nbnd, mask, size, nsteps, step_size, offset, table)
  # integer,                          intent(in) :: ncol, nlay, nbnd, nsteps
  # logical,  dimension(ncol,  nlay), intent(in) :: mask
  # real(FT), dimension(ncol,  nlay), intent(in) :: size
  # real(FT),                         intent(in) :: step_size, offset
  # real(FT), dimension(nsteps,nbnd), intent(in) :: table
  # real(FT), dimension(ncol,nlay,nbnd)          :: compute_from_table
  # # ---------------------------
  # integer  :: icol, ilay, ibnd
  # integer  :: index
  # real(FT) :: fint
  # # ---------------------------
  for ilay = 1:nlay
    for icol = 1:ncol
      if(mask[icol,ilay])
        index = min(floor((size(icol,ilay) - offset)/step_size)+1, nsteps-1)
        fint = (size(icol,ilay) - offset)/step_size - (index-1)
        for ibnd = 1:nbnd
          compute_from_table[icol,ilay,ibnd] = table[index,  ibnd] +
                                       fint * (table[index+1,ibnd] - table[index,ibnd])
        end
      else
        for ibnd = 1:nbnd
          compute_from_table[icol,ilay,ibnd] = FT(0)
        end
      end
    end
  end
end

#
# Pade functions
#

function compute_from_pade(ncol, nlay, nbnd, mask, size, nsizes, size_bounds, m, n, pade_coeffs)
  # integer,                        intent(in) :: ncol, nlay, nbnd, nsizes
  # logical,  dimension(ncol,nlay), intent(in) :: mask
  # real(FT), dimension(ncol,nlay), intent(in) :: size
  # real(FT), dimension(nsizes+1),  intent(in) :: size_bounds
  # integer,                        intent(in) :: m, n
  # real(FT), dimension(nbnd,nsizes,0:m+n),
  #                                 intent(in) :: pade_coeffs
  # real(FT), dimension(ncol,nlay,nbnd)        :: compute_from_pade
  # # ---------------------------
  # integer  :: icol, ilay, ibnd, irad

  for ilay = 1:nlay
    for icol = 1:ncol
      if mask[icol,ilay]
        #
        # Finds index into size regime table
        # This works only if there are precisely three size regimes (four bounds) and it's
        #   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
        #
        irad = min(floor((size(icol,ilay) - size_bounds(2))/size_bounds(3))+2, 3)
        compute_from_pade[icol,ilay,1:nbnd] =
             pade_eval(nbnd, nsizes, m, n, irad, size(icol,ilay), pade_coeffs)
      else
        for ibnd = 1:nbnd
          compute_from_pade[icol,ilay,ibnd] = FT(0)
        end
      end
    end
  end

end

#
# Evaluate Pade approximant of order [m/n]
#
function pade_eval(nbnd, nrads, m, n, irad, re, pade_coeffs::AbstractArray{FT}) where FT
  # integer,                intent(in) :: nbnd, nrads, m, n, irad
  # real(FT), dimension(nbnd, nrads, 0:m+n),
  #                         intent(in) :: pade_coeffs
  # real(FT),               intent(in) :: re
  # real(FT), dimension(nbnd)         :: pade_eval

  # integer :: iband
  # real(FT) :: numer, denom
  # integer  :: i
  res = Array{FT}(undef, nbnd)

  for iband = 1:nbnd
    denom = pade_coeffs[iband,irad,n+m]
    for i = n-1+m:-1:1+m
      denom = pade_coeffs[iband,irad,i]+re*denom
    end
    denom =  FT(1)                     +re*denom

    numer = pade_coeffs[iband,irad,m]
    for i = m-1:-1:1
      numer = pade_coeffs[iband,irad,i]+re*numer
    end
    numer = pade_coeffs[iband,irad,0]  +re*numer

    res[iband] = numer/denom
  end
  return res
end

end
