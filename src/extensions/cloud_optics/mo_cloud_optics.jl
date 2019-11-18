"""
    mo_cloud_optics

Provides cloud optical properties as a function of effective radius for the RRTMGP bands
  Based on Mie calculations for liquid
    and results from doi:10.1175/JAS-D-12-039.1 for ice with variable surface roughness
  Can use either look-up tables or Pade approximates according to which data has been loaded
  Mike Iacono (AER) is the original author
The class can be used as-is but is also intended as an example of how to extend the RTE framework
"""
module mo_cloud_optics

using OffsetArrays
using ..mo_optical_props
using ..mo_util_array
using ..fortran_intrinsics

export set_ice_roughness!, load_lut!, load_Pade!
export cloud_optics!
export get_min_radius_liq,
       get_max_radius_liq,
       get_min_radius_ice,
       get_max_radius_ice

export ty_cloud_optics

export PadeMethod
export LookUpTable
export ty_cloud_optics_new

abstract type AbstractInterpolationMethod end

struct PadeMethod <: AbstractInterpolationMethod
  # Pade approximation coefficients
  ext
  ssa
  asy
  # Particle size regimes for Pade formulations
  sizreg_ext
  sizreg_ssa
  sizreg_asy
end

struct LookUpTable{FT,I} <: AbstractInterpolationMethod
  rad_lwr::FT
  rad_upr::FT
  rad_fac::FT
  nsteps::I
  step_size::FT
  ext::Array{FT}
  ssa::Array{FT}
  asy::Array{FT}
  function LookUpTable(rad_lwr::FT,
                       rad_upr::FT,
                       rad_fac::FT,
                       ext::Array{FT},
                       ssa::Array{FT},
                       asy::Array{FT}) where {FT<:AbstractFloat}
    nsteps = size(ext, 1)
    step_size = (rad_upr - rad_lwr)/FT(nsteps-1)
    return new{FT,Int}(rad_lwr,rad_upr,rad_fac,nsteps,step_size,ext,ssa,asy)
  end
end
get_min_radius(lut::LookUpTable) = lut.rad_lwr
get_max_radius(lut::LookUpTable) = lut.rad_upr

mutable struct ty_cloud_optics_new{FT, I} <: ty_optical_props{FT, I}
  base
  # Ice surface roughness category - needed for Yang (2013) ice optics parameterization
  icergh#::I
  liq#AbstractInterpolationMethod
  ice#AbstractInterpolationMethod
end
ty_cloud_optics_new(FT, I) = ty_cloud_optics_new{FT,I}(ntuple(i->nothing, 4)...)


mutable struct ty_cloud_optics{FT, I} <: ty_optical_props{FT, I}
  base
  # Ice surface roughness category - needed for Yang (2013) ice optics parameterization
  icergh#::I

  #
  # Lookup table information
  #
  radliq_lwr#::FT
  radice_lwr#::FT
  radliq_upr#::FT
  radice_upr#::FT
  liq_nsteps#::I
  ice_nsteps#::I

  liq_step_size#::FT
  ice_step_size#::FT

  #
  # The tables themselves.
  #
  lut_extliq # Array{FT,2} (nsize_liq, nbnd)
  lut_ssaliq # Array{FT,2} (nsize_liq, nbnd)
  lut_asyliq # Array{FT,2} (nsize_liq, nbnd)
  lut_extice # Array{FT,3} # (nsize_ice, nbnd, nrghice)
  lut_ssaice # Array{FT,3} # (nsize_ice, nbnd, nrghice)
  lut_asyice # Array{FT,3} # (nsize_ice, nbnd, nrghice)

  #
  # Pade approximant coefficients
  #
  pade_extliq # Array{FT,3} # (nbnd, nsizereg, ncoeff_ext)
  pade_ssaliq # Array{FT,3} # (nbnd, nsizereg, ncoeff_ssa_g)
  pade_asyliq # Array{FT,3} # (nbnd, nsizereg, ncoeff_ssa_g)
  pade_extice # Array{FT,4} # (nbnd, nsizereg, ncoeff_ext, nrghice)
  pade_ssaice # Array{FT,4} # (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  pade_asyice # Array{FT,4} # (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
  # Particle size regimes for Pade formulations
  pade_sizreg_extliq # Array{FT,1}  # (nbound)
  pade_sizreg_ssaliq # Array{FT,1}  # (nbound)
  pade_sizreg_asyliq # Array{FT,1}  # (nbound)
  pade_sizreg_extice # Array{FT,1}  # (nbound)
  pade_sizreg_ssaice # Array{FT,1}  # (nbound)
  pade_sizreg_asyice # Array{FT,1}  # (nbound)
end

ty_cloud_optics(FT, I) = ty_cloud_optics{FT,I}(ntuple(i->nothing, 28)...)

# type, extends(ty_optical_props), public :: ty_cloud_optics
#   # private
#   #
#   # Ice surface roughness category - needed for Yang (2013) ice optics parameterization
#   #
#   integer            :: icergh = 0  # (1 = none, 2 = medium, 3 = high)
#   #
#   # Lookup table information
#   #
#   # Upper and lower limits of the tables
#   real(FT) :: radliq_lwr = FT(0), radliq_upr = FT(0)
#   real(FT) :: radice_lwr = FT(0), radice_upr = FT(0)
#   # How many steps in the table? (for convenience)
#   integer  :: liq_nsteps = 0,        ice_nsteps = 0
#   # How big is each step in the table?
#   real(FT) :: liq_step_size = FT(0), ice_step_size = FT(0)
#   #
#   # The tables themselves.
#   #
#   real(FT), dimension(:,:    ), allocatable :: lut_extliq, lut_ssaliq, lut_asyliq # (nsize_liq, nbnd)
#   real(FT), dimension(:,:,:  ), allocatable :: lut_extice, lut_ssaice, lut_asyice # (nsize_ice, nbnd, nrghice)

#   #
#   # Pade approximant coefficients
#   #
#   real(FT), dimension(:,:,:  ), allocatable :: pade_extliq                 # (nbnd, nsizereg, ncoeff_ext)
#   real(FT), dimension(:,:,:  ), allocatable :: pade_ssaliq,  pade_asyliq   # (nbnd, nsizereg, ncoeff_ssa_g)
#   real(FT), dimension(:,:,:,:), allocatable :: pade_extice                 # (nbnd, nsizereg, ncoeff_ext, nrghice)
#   real(FT), dimension(:,:,:,:), allocatable :: pade_ssaice, pade_asyice    # (nbnd, nsizereg, ncoeff_ssa_g, nrghice)
#   # Particle size regimes for Pade formulations
#   real(FT), dimension(:), allocatable :: pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq  # (nbound)
#   real(FT), dimension(:), allocatable :: pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice  # (nbound)
#   # -----
# # contains
# #   generic,   public :: load  => load_lut, load_pade
# #   procedure, public :: finalize
# #   procedure, public :: cloud_optics
# #   procedure, public :: get_min_radius_liq
# #   procedure, public :: get_min_radius_ice
# #   procedure, public :: get_max_radius_liq
# #   procedure, public :: get_max_radius_ice
# #   procedure, public :: get_num_ice_roughness_types
# #   procedure, public :: set_ice_roughness
# #   # Internal procedures
# #   procedure, private :: load_lut
# #   procedure, private :: load_pade
# # end type ty_cloud_optics

# ------------------------------------------------------------------------------
#
# Routines to load data needed for cloud optics calculations. Two routines: one to load
#    lookup-tables and one for coefficients for Pade approximates
#
# ------------------------------------------------------------------------------
function load_lut!(this::ty_cloud_optics_new{FT},
                   band_lims_wvn,
                   liq::LookUpTable,
                   ice::LookUpTable) where FT
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
  # # -------
  # #
  # # Local variables
  # #
  # integer               :: nbnd, nrghice, nsize_liq, nsize_ice

  this.base = ty_optical_props_base("RRTMGP cloud optics", band_lims_wvn)
  #
  # LUT coefficient dimensions
  #
  nsize_liq = size(liq.ext,1)
  nsize_ice = size(ice.ext,1)
  nbnd      = size(liq.ext,2)
  nrghice   = size(ice.ext,3)
  #
  # Error checking
  #   Can we check for consistency between table bounds and _fac?
  #
  @assert nbnd == get_nband(this)
  @assert size(ice.ext, 2) == nbnd
  @assert all(size(liq.ssa) .== (nsize_liq, nbnd))
  @assert all(size(liq.asy) .== (nsize_liq, nbnd))
  @assert all(size(ice.ssa) .== (nsize_ice, nbnd, nrghice))
  @assert all(size(ice.asy) .== (nsize_ice, nbnd, nrghice))

  this.liq = liq
  this.ice = ice
end
# ------------------------------------------------------------------------------
#
# Cloud optics initialization function - Pade
#
# ------------------------------------------------------------------------------
function load_Pade!(this::ty_cloud_optics{FT},
                    band_lims_wvn,
                    pade_extliq,
                    pade_ssaliq,
                    pade_asyliq,
                    pade_extice,
                    pade_ssaice,
                    pade_asyice,
                    pade_sizreg_extliq,
                    pade_sizreg_ssaliq,
                    pade_sizreg_asyliq,
                    pade_sizreg_extice,
                    pade_sizreg_ssaice,
                    pade_sizreg_asyice) where FT
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

  # integer               :: nbnd, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g, nbound

  # Pade coefficient dimensions
  nbnd         = size(pade_extliq, 1)
  nsizereg     = size(pade_extliq, 2)
  ncoeff_ext   = size(pade_extliq, 3)
  ncoeff_ssa_g = size(pade_ssaliq, 3)
  nrghice      = size(pade_extice, 4)
  nbound       = length(pade_sizreg_extliq)
  this.base = ty_optical_props_base("RRTMGP cloud optics", band_lims_wvn)

  # Error checking
  @assert nbnd == get_nband(this)
  @assert all(size(pade_ssaliq) .== (nbnd, nsizereg, ncoeff_ssa_g))
  @assert all(size(pade_asyliq) .== (nbnd, nsizereg, ncoeff_ssa_g))
  @assert all(size(pade_extice) .== (nbnd, nsizereg, ncoeff_ext, nrghice))
  @assert all(size(pade_ssaice) .== (nbnd, nsizereg, ncoeff_ssa_g, nrghice))
  @assert all(size(pade_asyice) .== (nbnd, nsizereg, ncoeff_ssa_g, nrghice))
  @assert length(pade_sizreg_ssaliq) == nbound
  @assert length(pade_sizreg_asyliq) == nbound
  @assert length(pade_sizreg_extice) == nbound
  @assert length(pade_sizreg_ssaice) == nbound
  @assert length(pade_sizreg_asyice) == nbound
  @assert nsizereg == 3

  #
  # Consistency among size regimes
  #
  this.radliq_lwr = pade_sizreg_extliq[1]
  this.radliq_upr = pade_sizreg_extliq[nbound]
  this.radice_lwr = pade_sizreg_extice[1]
  this.radice_upr = pade_sizreg_extice[nbound]

  @assert !any([pade_sizreg_ssaliq[1], pade_sizreg_asyliq[1]] .< this.radliq_lwr)
  @assert !any([pade_sizreg_ssaice[1], pade_sizreg_asyice[1]] .< this.radice_lwr)
  @assert !any([pade_sizreg_ssaliq[nbound], pade_sizreg_asyliq[nbound]] .> this.radliq_upr)
  @assert !any([pade_sizreg_ssaice[nbound], pade_sizreg_asyice[nbound]] .> this.radice_upr)

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


#
# Derive cloud optical properties from provided cloud physical properties

#
# Compute single-scattering properties
#
function cloud_optics!(this::ty_cloud_optics{FT},
                      clwp, ciwp, reliq, reice,
                      optical_props) where FT
    # class(ty_cloud_optics),
    #           intent(in   ) :: this
    # real(wp), intent(in   ) :: clwp  (:,:),    ! cloud ice water path    (units?)
    #                            ciwp  (:,:),    ! cloud liquid water path (units?)
    #                            reliq (:,:),    ! cloud ice particle effective size (microns)
    #                            reice (:,:)      ! cloud liquid particle effective radius (microns)
    # class(ty_optical_props_arry),
    #           intent(inout) :: optical_props
    #                                            ! Dimensions: (ncol,nlay,nbnd)

    # ! ------- Local -------
    # logical(wl), dimension(size(clwp,1), size(clwp,2)) :: liqmsk, icemsk
    # real(wp),    dimension(size(clwp,1), size(clwp,2), this%get_nband()) ::
    #             ltau, ltaussa, ltaussag, itau, itaussa, itaussag

  # # ------- Local -------
  # type(ty_optical_props_2str) :: clouds_liq, clouds_ice
  # integer  :: nsizereg, ibnd, imom


  ncol = size(clwp,1)
  nlay = size(clwp,2)
  nbnd = get_nband(this)
  liqmsk = BitArray(clwp .> FT(0))
  icemsk = BitArray(ciwp .> FT(0))

  ltau, ltaussa, ltaussag, itau, itaussa, itaussag = ntuple(i-> Array{FT}(undef, size(clwp,1), size(clwp,2), get_nband(this)), 6)

  # Error checking
  @assert allocated(this.lut_extliq) || allocated(this.pade_extliq)
  @assert bands_are_equal(this, optical_props)
  @assert get_nband(optical_props) == get_ngpt(optical_props)
  @assert !(this.icergh < 1 || this.icergh > get_num_ice_roughness_types(this))
  @assert !any_vals_outside(reliq, liqmsk, this.radliq_lwr, this.radliq_upr)
  @assert !any_vals_outside(reice, icemsk, this.radice_lwr, this.radice_upr)
  any(liqmsk) && @assert !any_vals_less_than(clwp, liqmsk, FT(0))
  any(icemsk) && @assert !any_vals_less_than(ciwp, icemsk, FT(0))

  #
  # Compute cloud optical properties. Use lookup tables if available, Pade approximants if not.
  #

  #
  # Cloud optical properties from Pade coefficient method
  #   Hard coded assumptions: order of approximants, three size regimes
  #
  nsizereg = size(this.pade_extliq,2)
  compute_all_from_pade!(ncol, nlay, nbnd, nsizereg,
                            liqmsk, clwp, reliq,
                            2, 3, this.pade_sizreg_extliq, this.pade_extliq,
                            2, 2, this.pade_sizreg_ssaliq, this.pade_ssaliq,
                            2, 2, this.pade_sizreg_asyliq, this.pade_asyliq,
                            ltau, ltaussa, ltaussag)
  compute_all_from_pade!(ncol, nlay, nbnd, nsizereg,
                            icemsk, ciwp, reice,
                            2, 3, this.pade_sizreg_extice, this.pade_extice[:,:,:,this.icergh],
                            2, 2, this.pade_sizreg_ssaice, this.pade_ssaice[:,:,:,this.icergh],
                            2, 2, this.pade_sizreg_asyice, this.pade_asyice[:,:,:,this.icergh],
                            itau, itaussa, itaussag)

  @assert optical_props isa ty_optical_props_1scl || optical_props isa ty_optical_props_2str
  #
  # Copy total cloud properties onto outputs
  #
  if optical_props isa ty_optical_props_1scl
    for ibnd = 1:nbnd
      for ilay = 1:nlay
        for icol = 1:ncol
          # Absorption optical depth  = (1-ssa) * tau = tau - taussa
          optical_props.tau[icol,ilay,ibnd] = (ltau[icol,ilay,ibnd] - ltaussa[icol,ilay,ibnd]) +
                                              (itau[icol,ilay,ibnd] - itaussa[icol,ilay,ibnd])
        end
      end
    end
  elseif optical_props isa ty_optical_props_2str
    for ibnd = 1:nbnd
      for ilay = 1:nlay
        for icol = 1:ncol
          tau    = ltau[icol,ilay,ibnd]    + itau[icol,ilay,ibnd]
          taussa = ltaussa[icol,ilay,ibnd] + itaussa[icol,ilay,ibnd]
          optical_props.g[icol,ilay,ibnd]  = (ltaussag[icol,ilay,ibnd] + itaussag[icol,ilay,ibnd]) / max(eps(FT), taussa...)
          optical_props.ssa[icol,ilay,ibnd] = taussa/max(eps(FT), tau...)
          optical_props.tau[icol,ilay,ibnd] = tau
        end
      end
    end
  end

end
function cloud_optics!(this::ty_cloud_optics_new{FT},
                      clwp, ciwp, reliq, reice,
                      optical_props) where FT
    # class(ty_cloud_optics),
    #           intent(in   ) :: this
    # real(wp), intent(in   ) :: clwp  (:,:),    ! cloud ice water path    (units?)
    #                            ciwp  (:,:),    ! cloud liquid water path (units?)
    #                            reliq (:,:),    ! cloud ice particle effective size (microns)
    #                            reice (:,:)      ! cloud liquid particle effective radius (microns)
    # class(ty_optical_props_arry),
    #           intent(inout) :: optical_props
    #                                            ! Dimensions: (ncol,nlay,nbnd)

    # ! ------- Local -------
    # logical(wl), dimension(size(clwp,1), size(clwp,2)) :: liqmsk, icemsk
    # real(wp),    dimension(size(clwp,1), size(clwp,2), this%get_nband()) ::
    #             ltau, ltaussa, ltaussag, itau, itaussa, itaussag

  # # ------- Local -------
  # type(ty_optical_props_2str) :: clouds_liq, clouds_ice
  # integer  :: nsizereg, ibnd, imom


  ncol = size(clwp,1)
  nlay = size(clwp,2)
  nbnd = get_nband(this)
  liqmsk = BitArray(clwp .> FT(0))
  icemsk = BitArray(ciwp .> FT(0))

  ltau, ltaussa, ltaussag, itau, itaussa, itaussag = ntuple(i-> Array{FT}(undef, size(clwp,1), size(clwp,2), get_nband(this)), 6)

  # Error checking
  @assert allocated(this.liq.ext)
  @assert bands_are_equal(this, optical_props)
  @assert get_nband(optical_props) == get_ngpt(optical_props)
  @assert !(this.icergh < 1 || this.icergh > get_num_ice_roughness_types(this))
  @assert !any_vals_outside(reliq, liqmsk, this.liq.rad_lwr, this.liq.rad_upr)
  @assert !any_vals_outside(reice, icemsk, this.ice.rad_lwr, this.ice.rad_upr)
  any(liqmsk) && @assert !any_vals_less_than(clwp, liqmsk, FT(0))
  any(icemsk) && @assert !any_vals_less_than(ciwp, icemsk, FT(0))

  #
  # Compute cloud optical properties. Use lookup tables if available, Pade approximants if not.
  #

  #
  # Liquid
  #
  compute_all_from_table!(ncol, nlay, nbnd, liqmsk, clwp, reliq,
                              this.liq.nsteps,
                              this.liq.step_size,
                              this.liq.rad_lwr,
                              this.liq.ext,
                              this.liq.ssa,
                              this.liq.asy,
                              ltau, ltaussa, ltaussag)
  #
  # Ice
  #
  compute_all_from_table!(ncol, nlay, nbnd, icemsk, ciwp, reice,
                              this.ice.nsteps,
                              this.ice.step_size,
                              this.ice.rad_lwr,
                              this.ice.ext[:,:,this.icergh],
                              this.ice.ssa[:,:,this.icergh],
                              this.ice.asy[:,:,this.icergh],
                              itau, itaussa, itaussag)

  @assert optical_props isa ty_optical_props_1scl || optical_props isa ty_optical_props_2str
  #
  # Copy total cloud properties onto outputs
  #
  if optical_props isa ty_optical_props_1scl
    for ibnd = 1:nbnd
      for ilay = 1:nlay
        for icol = 1:ncol
          # Absorption optical depth  = (1-ssa) * tau = tau - taussa
          optical_props.tau[icol,ilay,ibnd] = (ltau[icol,ilay,ibnd] - ltaussa[icol,ilay,ibnd]) +
                                              (itau[icol,ilay,ibnd] - itaussa[icol,ilay,ibnd])
        end
      end
    end
  elseif optical_props isa ty_optical_props_2str
    for ibnd = 1:nbnd
      for ilay = 1:nlay
        for icol = 1:ncol
          tau    = ltau[icol,ilay,ibnd]    + itau[icol,ilay,ibnd]
          taussa = ltaussa[icol,ilay,ibnd] + itaussa[icol,ilay,ibnd]
          optical_props.g[icol,ilay,ibnd]  = (ltaussag[icol,ilay,ibnd] + itaussag[icol,ilay,ibnd]) / max(eps(FT), taussa...)
          optical_props.ssa[icol,ilay,ibnd] = taussa/max(eps(FT), tau...)
          optical_props.tau[icol,ilay,ibnd] = tau
        end
      end
    end
  end

end

"""
    set_ice_roughness!

"""
function set_ice_roughness!(this::ty_optical_props, icergh)
  @assert icergh > 0
  this.icergh = icergh
end

get_num_ice_roughness_types(this::ty_cloud_optics) =
  allocated(this.pade_extice) ? size(this.pade_extice, 4) : 0
get_num_ice_roughness_types(this::ty_cloud_optics_new) =
  allocated(this.ice.ext ) ? size(this.ice.ext,  3) : 0

get_min_radius_liq(this::ty_cloud_optics) = this.radliq_lwr
get_max_radius_liq(this::ty_cloud_optics) = this.radliq_upr
get_min_radius_ice(this::ty_cloud_optics) = this.radice_lwr
get_max_radius_ice(this::ty_cloud_optics) = this.radice_upr
get_min_radius_liq(this::ty_cloud_optics_new) = get_min_radius(this.liq)
get_max_radius_liq(this::ty_cloud_optics_new) = get_max_radius(this.liq)
get_min_radius_ice(this::ty_cloud_optics_new) = get_min_radius(this.ice)
get_max_radius_ice(this::ty_cloud_optics_new) = get_max_radius(this.ice)


#
# Ancillary functions
#

#
# Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
#   elements starting at "offset." The table's second dimension is band.
# Returns 0 where the mask is false.
# We could also try gather/scatter for efficiency
#
function compute_all_from_table!(ncol, nlay, nbnd, mask, lwp, re,
                                  nsteps, step_size, offset,
                                  tau_table, ssa_table, asy_table,
                                  tau::Array{FT}, taussa, taussag) where FT
  # integer,                                intent(in) :: ncol, nlay, nbnd, nsteps
  # logical(wl), dimension(ncol,nlay),      intent(in) :: mask
  # real(wp),    dimension(ncol,nlay),      intent(in) :: lwp, re
  # real(wp),                               intent(in) :: step_size, offset
  # real(wp),    dimension(nsteps,   nbnd), intent(in) :: tau_table, ssa_table, asy_table
  # real(wp),    dimension(ncol,nlay,nbnd)             :: tau, taussa, taussag
  # # ---------------------------
  # integer  :: icol, ilay, ibnd
  # integer  :: index
  # real(wp) :: fint
  # real(wp) :: t, ts, tsg  # tau, tau*ssa, tau*ssa*g
  # # ---------------------------
  #$acc parallel loop gang vector default(present) collapse(3)
  for ibnd = 1:nbnd
    for ilay = 1:nlay
      for icol = 1:ncol
        if mask[icol,ilay]
          index = convert(Int, min(floor((re[icol,ilay] - offset)/step_size)+1, nsteps-1))
          fint = (re[icol,ilay] - offset)/step_size - (index-1)
          t   = lwp[icol,ilay] *                  (tau_table[index,  ibnd] + fint * (tau_table[index+1,ibnd] - tau_table[index,ibnd]))
          ts  = t              *                  (ssa_table[index,  ibnd] + fint * (ssa_table[index+1,ibnd] - ssa_table[index,ibnd]))
          taussag[icol,ilay,ibnd] = ts          * (asy_table[index,  ibnd] + fint * (asy_table[index+1,ibnd] - asy_table[index,ibnd]))
          taussa[icol,ilay,ibnd] = ts
          tau[icol,ilay,ibnd] = t
        else
          tau[icol,ilay,ibnd] = FT(0)
          taussa[icol,ilay,ibnd] = FT(0)
          taussag[icol,ilay,ibnd] = FT(0)
        end
      end
    end
  end
end

#
# Pade functions
#

function compute_all_from_pade!(ncol, nlay, nbnd, nsizes,
                                 mask, lwp, re,
                                 m_ext, n_ext, re_bounds_ext, coeffs_ext,
                                 m_ssa, n_ssa, re_bounds_ssa, coeffs_ssa,
                                 m_asy, n_asy, re_bounds_asy, coeffs_asy,
                                 tau::Array{FT}, taussa, taussag) where FT
  # integer,                        intent(in) :: ncol, nlay, nbnd, nsizes
  # logical(wl),
  #           dimension(ncol,nlay), intent(in) :: mask
  # real(wp), dimension(ncol,nlay), intent(in) :: lwp, re
  # real(wp), dimension(nsizes+1),  intent(in) :: re_bounds_ext, re_bounds_ssa, re_bounds_asy
  # integer,                        intent(in) :: m_ext, n_ext
  # real(wp), dimension(nbnd,nsizes,0:m_ext+n_ext),
  #                                 intent(in) :: coeffs_ext
  # integer,                        intent(in) :: m_ssa, n_ssa
  # real(wp), dimension(nbnd,nsizes,0:m_ssa+n_ssa),
  #                                 intent(in) :: coeffs_ssa
  # integer,                        intent(in) :: m_asy, n_asy
  # real(wp), dimension(nbnd,nsizes,0:m_asy+n_asy),
  #                                 intent(in) :: coeffs_asy
  # real(wp), dimension(ncol,nlay,nbnd)        :: tau, taussa, taussag
  # # ---------------------------
  # integer  :: icol, ilay, ibnd, irad, count
  # real(wp) :: t, ts

  #$acc parallel loop gang vector default(present) collapse(3)
  for ibnd = 1:nbnd
    for ilay = 1:nlay
      for icol = 1:ncol
        if mask[icol,ilay]
          #
          # Finds index into size regime table
          # This works only if there are precisely three size regimes (four bounds) and it's
          #   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
          #
          irad = convert(Int, min(floor((re[icol,ilay] - re_bounds_ext[2])/re_bounds_ext[3])+2, 3))
          t   = lwp[icol,ilay] * pade_eval(ibnd, nbnd, nsizes, m_ext, n_ext, irad, re[icol,ilay], coeffs_ext)

          irad = convert(Int,min(floor((re[icol,ilay] - re_bounds_ssa[2])/re_bounds_ssa[3])+2, 3))
          # Pade approximants for co-albedo can sometimes be negative
          ts   = t               * (FT(1) - max(FT(0), pade_eval(ibnd, nbnd, nsizes, m_ssa, n_ssa, irad, re[icol,ilay], coeffs_ssa)))
          irad = convert(Int, min(floor((re[icol,ilay] - re_bounds_asy[2])/re_bounds_asy[3])+2, 3))
          taussag[icol,ilay,ibnd] = ts             * pade_eval(ibnd, nbnd, nsizes, m_asy, n_asy, irad, re[icol,ilay], coeffs_asy)

          taussa[icol,ilay,ibnd] = ts
          tau[icol,ilay,ibnd] = t
        else
          tau[icol,ilay,ibnd] = FT(0)
          taussa[icol,ilay,ibnd] = FT(0)
          taussag[icol,ilay,ibnd] = FT(0)
        end
      end
    end
  end

end

#
# Ancillary functions
#
#

#
# Pade functions
#

#
# Evaluate Pade approximant of order [m/n]
#
function pade_eval(nbnd, nrads, m, n, irad, re, pade_coeffs)
  # integer,                intent(in) :: nbnd, nrads, m, n, irad
  # real(wp), dimension(nbnd, nrads, 0:m+n), &
  #                         intent(in) :: pade_coeffs
  # real(wp),               intent(in) :: re
  # real(wp), dimension(nbnd)          :: pade_eval_nbnd

  # integer :: iband
  # real(wp) :: numer, denom
  # integer  :: i
  FT = eltype(pade_coeffs)
  res = Vector{FT}(undef, nbnd)

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

#
# Evaluate Pade approximant of order [m/n]
#
function pade_eval(iband, nbnd, nrads, m, n, irad, re, pade_coeffs_array)
  # !$acc routine seq
  # !
  # integer,                intent(in) :: iband, nbnd, nrads, m, n, irad
  # real(wp), dimension(nbnd, nrads, 0:m+n), &
  #                         intent(in) :: pade_coeffs
  # real(wp),               intent(in) :: re
  # real(wp)                           :: pade_eval_1

  # real(wp) :: numer, denom
  # integer  :: i
  FT = eltype(re)
  pade_coeffs = OffsetArray{FT}(undef, 1:nbnd, 1:nrads, 0:m+n)
  pade_coeffs[:,:,:] = pade_coeffs_array[:,:,:]

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

  pade_eval_1 = numer/denom
  return pade_eval_1
end

end
