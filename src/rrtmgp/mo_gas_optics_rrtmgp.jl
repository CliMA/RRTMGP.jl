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
#
# Class for computing spectrally-resolved gas optical properties and source functions
#   given atmopsheric physical properties (profiles of temperature, pressure, and gas concentrations)
#   The class must be initialized with data (provided as a netCDF file) before being used.
#
# Two variants apply to internal Planck sources (longwave radiation in the Earth's atmosphere) and to
#   external stellar radiation (shortwave radiation in the Earth's atmosphere).
#   The variant is chosen based on what information is supplied during initialization.
#   (It might make more sense to define two sub-classes)
#
# -------------------------------------------------------------------------------------------------
module mo_gas_optics_rrtmgp
  # use mo_rte_kind,           only: wp, wl
  # use mo_rrtmgp_constants,   only: avogad, m_dry, m_h2o, grav
  using ..mo_rrtmgp_constants
  # use mo_util_array,         only: zero_array, any_vals_less_than, any_vals_outside
  using ..mo_util_array
  # use mo_optical_props,      only: ty_optical_props
  using ..mo_optical_props
  # use mo_source_functions,   only: ty_source_func_lw
  using ..mo_source_functions
  # use mo_gas_optics_kernels, only: interpolation,
  #                                  compute_tau_absorption, compute_tau_rayleigh, compute_Planck_source,
  #                                  combine_and_reorder_2str, combine_and_reorder_nstr
  using ..mo_gas_optics_kernels

  using ..mo_util_string
  # use mo_gas_concentrations, only: ty_gas_concs
  using ..mo_gas_concentrations
  # use mo_optical_props,      only: ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  using ..mo_optical_props
  # use mo_gas_optics,         only: ty_gas_optics
  # using mo_gas_optics # only defines abstract interfaces
  using ..mo_util_reorder
  export gas_optics!, ty_gas_optics_rrtmgp

  # -------------------------------------------------------------------------------------------------
  # type, extends(ty_gas_optics), public :: ty_gas_optics_rrtmgp
  struct ty_gas_optics_rrtmgp{T,I} <: ty_optical_props{T,I}
    # private
    #
    # RRTMGP computes absorption in each band arising from
    #   two major species in each band, which are combined to make
    #     a relative mixing ratio eta and a total column amount (col_mix)
    #   contributions from zero or more minor species whose concentrations
    #     may be scaled by other components of the atmosphere
    #
    # Absorption coefficients are interpolated from tables on a pressure/temperature/(eta) grid
    #
    # ------------------------------------
    # Interpolation variables: Temperature and pressure grids
    #
    # real(wp),      dimension(:),     allocatable :: press_ref,  press_ref_log, temp_ref
    press_ref::Vector{T}
    press_ref_log::Vector{T}
    temp_ref::Vector{T}
    #
    # Derived and stored for convenience:
    #   Min and max for temperature and pressure intepolation grids
    #   difference in ln pressure between consecutive reference levels
    #   log of reference pressure separating the lower and upper atmosphere
    #
    # real(wp) :: press_ref_min, press_ref_max, temp_ref_min,  temp_ref_max
    press_ref_min::T
    press_ref_max::T
    temp_ref_min::T
    temp_ref_max::T

    # real(wp) :: press_ref_log_delta, temp_ref_delta, press_ref_trop_log
    press_ref_log_delta::T
    temp_ref_delta::T
    press_ref_trop_log::T
    # ------------------------------------
    # Major absorbers ("key species")
    #   Each unique set of major species is called a flavor.
    #
    # Names  and reference volume mixing ratios of major gases
    #
    # character(32), dimension(:),  allocatable :: gas_names     # gas names
    gas_names::Vector{String}     # gas names
    # real(wp), dimension(:,:,:),   allocatable :: vmr_ref       # vmr_ref(lower or upper atmosphere, gas, temp)
    vmr_ref::Array{T,3}       # vmr_ref(lower or upper atmosphere, gas, temp)
    #
    # Which two gases are in each flavor? By index
    #
    # integer,  dimension(:,:),     allocatable :: flavor        # major species pair; (2,nflav)
    flavor::Array{I, 2}        # major species pair; (2,nflav)
    #
    # Which flavor for each g-point? One each for lower, upper atmosphere
    #
    # integer,  dimension(:,:),     allocatable :: gpoint_flavor # flavor = gpoint_flavor(2, g-point)
    gpoint_flavor::Array{I, 2} # flavor = gpoint_flavor(2, g-point)
    #
    # Major gas absorption coefficients
    #
    # real(wp), dimension(:,:,:,:), allocatable :: kmajor        #  kmajor(g-point,eta,pressure,temperature)
    kmajor::Array{T,4}        #  kmajor(g-point,eta,pressure,temperature)
    #
    # ------------------------------------
    # Minor species, independently for upper and lower atmospheres
    #   Array extents in the n_minor dimension will differ between upper and lower atmospheres
    #   Each contribution has starting and ending g-points
    #
    # integer, dimension(:,:), allocatable :: minor_limits_gpt_lower,
    #                                         minor_limits_gpt_upper
    minor_limits_gpt_lower::Array{I,2}
    minor_limits_gpt_upper::Array{I,2}
    #
    # Minor gas contributions might be scaled by other gas amounts; if so we need to know
    #   the total density and whether the contribution is scaled by the partner gas
    #   or its complement (i.e. all other gases)
    # Water vapor self- and foreign continua work like this, as do
    #   all collision-induced abosption pairs
    #
    # logical(wl), dimension(:), allocatable :: minor_scales_with_density_lower,
    #                                           minor_scales_with_density_upper
    minor_scales_with_density_lower::Vector{Bool}
    minor_scales_with_density_upper::Vector{Bool}

    # logical(wl), dimension(:), allocatable :: scale_by_complement_lower, scale_by_complement_upper
    scale_by_complement_lower::Vector{Bool}
    scale_by_complement_upper::Vector{Bool}
    # integer,     dimension(:), allocatable :: idx_minor_lower,           idx_minor_upper
    idx_minor_lower::Vector{I}
    idx_minor_upper::Vector{I}
    # integer,     dimension(:), allocatable :: idx_minor_scaling_lower,   idx_minor_scaling_upper
    idx_minor_scaling_lower::Vector{I}
    idx_minor_scaling_upper::Vector{I}
    #
    # Index into table of absorption coefficients
    #
    # integer, dimension(:), allocatable :: kminor_start_lower,        kminor_start_upper
    kminor_start_lower::Vector{I}
    kminor_start_upper::Vector{I}
    #
    # The absorption coefficients themselves
    #
    # real(wp), dimension(:,:,:), allocatable :: kminor_lower, kminor_upper # kminor_lower(n_minor,eta,temperature)
    kminor_lower::Array{T,3}
    kminor_upper::Array{T,3} # kminor_lower(n_minor,eta,temperature)
    #
    # -----------------------------------------------------------------------------------
    #
    # Rayleigh scattering coefficients
    #
    # real(wp), dimension(:,:,:,:), allocatable :: krayl # krayl(g-point,eta,temperature,upper/lower atmosphere)
    krayl::Array{T, 4} # krayl(g-point,eta,temperature,upper/lower atmosphere)
    #
    # -----------------------------------------------------------------------------------
    # Planck function spectral mapping
    #   Allocated only when gas optics object is internal-source
    #
    # real(wp), dimension(:,:,:,:), allocatable :: planck_frac   # stored fraction of Planck irradiance in band for given g-point
    #                                                            # planck_frac(g-point, eta, pressure, temperature)
    planck_frac::Array{T, 4}   # stored fraction of Planck irradiance in band for given g-point
                                                               # planck_frac(g-point, eta, pressure, temperature)
    # real(wp), dimension(:,:),     allocatable :: totplnk       # integrated Planck irradiance by band; (Planck temperatures,band)
    totplnk::Array{T,2}       # integrated Planck irradiance by band; (Planck temperatures,band)
    # real(wp)                                  :: totplnk_delta # temperature steps in totplnk
    totplnk_delta::T # temperature steps in totplnk
    # -----------------------------------------------------------------------------------
    # Solar source function spectral mapping
    #   Allocated only when gas optics object is external-source
    #
    # real(wp), dimension(:), allocatable :: solar_src # incoming solar irradiance(g-point)
    solar_src::Vector{T} # incoming solar irradiance(g-point)
    #
    # -----------------------------------------------------------------------------------
    # Ancillary
    # -----------------------------------------------------------------------------------
    # Index into %gas_names -- is this a key species in any band?
    # logical, dimension(:), allocatable :: is_key
    is_key::Vector{Bool}
    # -----------------------------------------------------------------------------------

  # contains
  #   # Type-bound procedures
  #   # Public procedures
  #   # public interface
  #   generic,   public :: load       => load_int,       load_ext
  #   procedure, public :: source_is_internal
  #   procedure, public :: source_is_external
  #   procedure, public :: get_ngas
  #   procedure, public :: get_gases
  #   procedure, public :: get_press_min
  #   procedure, public :: get_press_max
  #   procedure, public :: get_temp_min
  #   procedure, public :: get_temp_max
  #   # Internal procedures
  #   procedure, private :: load_int
  #   procedure, private :: load_ext
  #   procedure, public  :: gas_optics_int
  #   procedure, public  :: gas_optics_ext
  #   procedure, private :: check_key_species_present
  #   procedure, private :: get_minor_list
  #   # Interpolation table dimensions
  #   procedure, private :: get_nflav
  #   procedure, private :: get_neta
  #   procedure, private :: get_npres
  #   procedure, private :: get_ntemp
  #   procedure, private :: get_nPlanckTemp
  end

  # -------------------------------------------------------------------------------------------------
  #
  # col_dry is the number of molecules per cm-2 of dry air
  #
  # public :: get_col_dry # Utility function, not type-bound
  export get_col_dry

#   interface check_range
#     module procedure check_range_1D, check_range_2D, check_range_3D
#   end interface check_range

#   interface check_extent
#     module procedure check_extent_1D, check_extent_2D, check_extent_3D
#     module procedure check_extent_4D, check_extent_5D, check_extent_6D
#   end interface check_extent
# contains
  # --------------------------------------------------------------------------------------
  #
  # Public procedures
  #
  # --------------------------------------------------------------------------------------
  #
  # Two functions to define array sizes needed by gas_optics()
  #
  function get_ngas(this::ty_gas_optics_rrtmgp)
    # return the number of gases registered in the spectral configuration
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # integer                                        :: get_ngas

    return size(this.gas_names)
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # return the number of distinct major gas pairs in the spectral bands (referred to as
  # "flavors" - all bands have a flavor even if there is one or no major gas)
  #
  function get_nflav(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # integer                                 :: get_nflav

    return size(this.flavor, 2)
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Compute gas optical depth and Planck source functions,
  #  given temperature, pressure, and composition
  #
  function gas_optics!(this::ty_gas_optics_rrtmgp,
                          play, plev, tlay, tsfc, gas_desc::ty_gas_concs,
                          optical_props::ty_optical_props_arry, sources::ty_source_func_lw,
                          col_dry=nothing, tlev=nothing) # result(error_msg)
    # inputs
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # real(wp), dimension(:,:), intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
    #                                            plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
    #                                            tlay      # layer temperatures [K]; (ncol,nlay)
    # real(wp), dimension(:),   intent(in   ) :: tsfc      # surface skin temperatures [K]; (ncol)
    # type(ty_gas_concs),       intent(in   ) :: gas_desc  # Gas volume mixing ratios
    # # output
    # class(ty_optical_props_arry),
    #                           intent(inout) :: optical_props # Optical properties
    # class(ty_source_func_lw    ),
    #                           intent(inout) :: sources       # Planck sources
    # character(len=128)                      :: error_msg
    # # Optional inputs
    # real(wp), dimension(:,:),   intent(in   ),
    #                        optional, target :: col_dry,   # Column dry amount; dim(ncol,nlay)
    #                                            tlev        # level temperatures [K]; (ncol,nlay+1)
    # ----------------------------------------------------------
    # Local variables
    # Interpolation coefficients for use in source function
    # integer,     dimension(size(play,dim=1), size(play,dim=2)) :: jtemp, jpress
    # logical(wl), dimension(size(play,dim=1), size(play,dim=2)) :: tropo
    # real(wp),    dimension(2,2,2,get_nflav(this),size(play,dim=1), size(play,dim=2)) :: fmajor
    # integer,     dimension(2,    get_nflav(this),size(play,dim=1), size(play,dim=2)) :: jeta

    # integer :: ncol, nlay, ngpt, nband, ngas, nflav
    # # ----------------------------------------------------------
    DT = eltype(play)

    jpress = Array{Int}(undef, size(play))
    jtemp = Array{Int}(undef, size(play))
    tropo = Array{Bool}(undef, size(play))
    fmajor = Array{DT}(undef, 2,2,2,get_nflav(this),size(play)...)
    jeta = Array{Int}(undef, 2,    get_nflav(this), size(play)...)

    ncol  = size(play, 1)
    nlay  = size(play, 2)
    ngpt  = get_ngpt(this)
    nband = get_nband(this)
    #
    # Gas optics
    #
    #$acc enter data create(jtemp, jpress, tropo, fmajor, jeta)
    error_msg = compute_gas_taus(this,
                                 ncol, nlay, ngpt, nband,
                                 play, plev, tlay, gas_desc,
                                 optical_props,
                                 jtemp, jpress, jeta, tropo, fmajor,
                                 col_dry)
    error_msg  ≠ "" && return error_msg

    # ----------------------------------------------------------
    #
    # External source -- check arrays sizes and values
    # input data sizes and values
    #
    error_msg = check_extent(tsfc, ncol, "tsfc")
    error_msg  ≠ "" && return error_msg
    error_msg = check_range(tsfc, this%temp_ref_min,  this%temp_ref_max,  "tsfc")
    error_msg  ≠ "" && return error_msg
    if present(tlev)
      error_msg = check_extent(tlev, ncol, nlay+1, "tlev")
      error_msg  ≠ "" && return error_msg
      error_msg = check_range(tlev, this%temp_ref_min, this%temp_ref_max, "tlev")
      error_msg  ≠ "" && return error_msg
    end

    #
    #   output extents
    #
    if any([get_ncol(sources), get_nlay(sources), get_ngpt(sources)] .≠ [ncol, nlay, ngpt])
      error_msg = "gas_optics%gas_optics: source function arrays inconsistently sized"
    end
    error_msg  ≠ "" && return error_msg

    #
    # Interpolate source function
    #
    error_msg = source(this,
                       ncol, nlay, nband, ngpt,
                       play, plev, tlay, tsfc,
                       jtemp, jpress, jeta, tropo, fmajor,
                       sources,
                       tlev)
    #$acc exit data delete(jtemp, jpress, tropo, fmajor, jeta)
  end
  #------------------------------------------------------------------------------------------
  #
  # Compute gas optical depth given temperature, pressure, and composition
  #
  function gas_optics!(this::ty_gas_optics_rrtmgp,
                          play, plev, tlay, gas_desc::ty_gas_concs,    # mandatory inputs
                          optical_props::ty_optical_props_arry, toa_src,        # mandatory outputs
                          col_dry=nothing) # result(error_msg)      # optional input

    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # real(wp), dimension(:,:), intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
    #                                            plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
    #                                            tlay      # layer temperatures [K]; (ncol,nlay)
    # type(ty_gas_concs),       intent(in   ) :: gas_desc  # Gas volume mixing ratios
    # # output
    # class(ty_optical_props_arry),
    #                           intent(inout) :: optical_props
    # real(wp), dimension(:,:), intent(  out) :: toa_src     # Incoming solar irradiance(ncol,ngpt)
    # character(len=128)                      :: error_msg

    # Optional inputs
    # real(wp), dimension(:,:), intent(in   ),
    #                        optional, target :: col_dry # Column dry amount; dim(ncol,nlay)
    # ----------------------------------------------------------
    # Local variables
    # Interpolation coefficients for use in source function
    # integer,     dimension(size(play,dim=1), size(play,dim=2)) :: jtemp, jpress
    # logical(wl), dimension(size(play,dim=1), size(play,dim=2)) :: tropo
    # real(wp),    dimension(2,2,2,get_nflav(this),size(play,dim=1), size(play,dim=2)) :: fmajor
    # integer,     dimension(2,    get_nflav(this),size(play,dim=1), size(play,dim=2)) :: jeta

    # integer :: ncol, nlay, ngpt, nband, ngas, nflav
    # integer :: igpt, icol
    # # ----------------------------------------------------------

    DT = eltype(play)

    jpress = Array{Int}(undef, size(play))
    jtemp = Array{Int}(undef, size(play))
    tropo = Array{Bool}(undef, size(play))
    fmajor = Array{DT}(undef, 2,2,2,get_nflav(this),size(play)...)
    jeta = Array{Int}(undef, 2,    get_nflav(this), size(play)...)

    ncol  = size(play, 1)
    nlay  = size(play, 2)
    ngpt  = get_ngpt(this)
    nband = get_nband(this)
    ngas  = get_ngas(this)
    nflav = get_nflav(this)
    #
    # Gas optics
    #
    #$acc enter data create(jtemp, jpress, tropo, fmajor, jeta)
    error_msg = compute_gas_taus(this,
                                 ncol, nlay, ngpt, nband,
                                 play, plev, tlay, gas_desc,
                                 optical_props,
                                 jtemp, jpress, jeta, tropo, fmajor,
                                 col_dry)
    #$acc exit data delete(jtemp, jpress, tropo, fmajor, jeta)
    error_msg  ≠ "" && return error_msg

    # ----------------------------------------------------------
    #
    # External source function is constant
    #
    error_msg = check_extent(toa_src,     ncol,         ngpt, "toa_src")
    error_msg  ≠ "" && return error_msg
    #$acc parallel loop collapse(2)
    for igpt in 1:ngpt
       for icol in 1:ncol
          toa_src[icol,igpt] = this.solar_src[igpt]
       end
    end
  end
  #------------------------------------------------------------------------------------------
  #
  # Returns optical properties and interpolation coefficients
  #
  function compute_gas_taus(this::ty_gas_optics_rrtmgp,
                            ncol, nlay, ngpt, nband,
                            play, plev, tlay, gas_desc::ty_gas_concs,
                            optical_props::ty_optical_props_arry,
                            jtemp, jpress, jeta, tropo, fmajor,
                            col_dry=nothing) # result(error_msg)

    # class(ty_gas_optics_rrtmgp),
    #                                   intent(in   ) :: this
    # integer,                          intent(in   ) :: ncol, nlay, ngpt, nband
    # real(wp), dimension(:,:),         intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
    #                                                    plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
    #                                                    tlay      # layer temperatures [K]; (ncol,nlay)
    # type(ty_gas_concs),               intent(in   ) :: gas_desc  # Gas volume mixing ratios
    # class(ty_optical_props_arry),     intent(inout) :: optical_props #inout because components are allocated
    # # Interpolation coefficients for use in internal source function
    # integer,     dimension(                       ncol, nlay), intent(  out) :: jtemp, jpress
    # integer,     dimension(2,    get_nflav(this),ncol, nlay), intent(  out) :: jeta
    # logical(wl), dimension(                       ncol, nlay), intent(  out) :: tropo
    # real(wp),    dimension(2,2,2,get_nflav(this),ncol, nlay), intent(  out) :: fmajor
    # character(len=128)                                         :: error_msg

    # Optional inputs
    # real(wp), dimension(:,:), intent(in   ),
    #                        optional, target :: col_dry # Column dry amount; dim(ncol,nlay)
    # ----------------------------------------------------------
    # Local variables
    # real(wp), dimension(ngpt,nlay,ncol) :: tau, tau_rayleigh  # absorption, Rayleigh scattering optical depths
    # # integer :: igas, idx_h2o # index of some gases
    # # Number of molecules per cm^2
    # real(wp), dimension(ncol,nlay), target  :: col_dry_arr
    # real(wp), dimension(:,:),       pointer :: col_dry_wk
    # #
    # # Interpolation variables used in major gas but not elsewhere, so don't need exporting
    # #
    # real(wp), dimension(ncol,nlay,  this%get_ngas()) :: vmr     # volume mixing ratios
    # real(wp), dimension(ncol,nlay,0:this%get_ngas()) :: col_gas # column amounts for each gas, plus col_dry
    # real(wp), dimension(2,    get_nflav(this),ncol,nlay) :: col_mix # combination of major species's column amounts
    #                                                      # index(1) : reference temperature level
    #                                                      # index(2) : flavor
    #                                                      # index(3) : layer
    # real(wp), dimension(2,2,  get_nflav(this),ncol,nlay) :: fminor # interpolation fractions for minor species
    #                                                       # index(1) : reference eta level (temperature dependent)
    #                                                       # index(2) : reference temperature level
    #                                                       # index(3) : flavor
    #                                                       # index(4) : layer
    # integer :: ngas, nflav, neta, npres, ntemp
    # integer :: nminorlower, nminorklower,nminorupper, nminorkupper
    # logical :: use_rayl
    # ----------------------------------------------------------
    DT = eltype(play)
    tau = Array{DT}(undef, ngpt,nlay,ncol)          # absorption, Rayleigh scattering optical depths
    tau_rayleigh = Array{DT}(undef, ngpt,nlay,ncol) # absorption, Rayleigh scattering optical depths
    col_dry_arr = Array{DT}(undef, ncol, nlay)
    col_dry_wk = Array{DT}(undef, ncol, nlay)

    vmr     = Array{DT}(undef, ncol,nlay,  get_ngas(this)) # volume mixing ratios
    col_gas = Array{DT}(undef, ncol,nlay,0:get_ngas(this)) # column amounts for each gas, plus col_dry
    col_mix = Array{DT}(undef, 2,    get_nflav(this),ncol,nlay) # combination of major species's column amounts
    fminor  = Array{DT}(undef, 2,2,  get_nflav(this),ncol,nlay) # interpolation fractions for minor species

    #
    # Error checking
    #
    use_rayl = allocated(this.krayl)
    error_msg = ""
    # Check for initialization
    if !is_initialized(this)
      error_msg = "ERROR: spectral configuration not loaded"
      return error_msg
    end
    #
    # Check for presence of key species in ty_gas_concs; return error if any key species are not present
    #
    error_msg = check_key_species_present(this, gas_desc)
    error_msg ≠ "" && return error_msg

    #
    # Check input data sizes and values
    #
    error_msg = check_extent(play, ncol, nlay,   "play")
    error_msg  ≠ "" && return error_msg
    error_msg = check_extent(plev, ncol, nlay+1, "plev")
    error_msg  ≠ "" && return error_msg
    error_msg = check_extent(tlay, ncol, nlay,   "tlay")
    error_msg  ≠ "" && return error_msg
    error_msg = check_range(play, this.press_ref_min,this.press_ref_max, "play")
    error_msg  ≠ "" && return error_msg
    error_msg = check_range(plev, this.press_ref_min, this.press_ref_max, "plev")
    error_msg  ≠ "" && return error_msg
    error_msg = check_range(tlay, this.temp_ref_min,  this.temp_ref_max,  "tlay")
    error_msg  ≠ "" && return error_msg
    if present(col_dry)
      error_msg = check_extent(col_dry, ncol, nlay, "col_dry")
      error_msg  ≠ "" && return error_msg
      error_msg = check_range(col_dry, DT(0), floatmax(DT), "col_dry")
      error_msg  ≠ "" && return error_msg
    end

    # ----------------------------------------------------------
    ngas  = get_ngas(this)
    nflav = get_nflav(this)
    neta  = get_neta(this)
    npres = get_npres(this)
    ntemp = get_ntemp(this)
    # number of minor contributors, total num absorption coeffs
    nminorlower  = size(this.minor_scales_with_density_lower)
    nminorklower = size(this.kminor_lower, 1)
    nminorupper  = size(this.minor_scales_with_density_upper)
    nminorkupper = size(this.kminor_upper, 1)
    #
    # Fill out the array of volume mixing ratios
    #
    for igas in 1:ngas
      #
      # Get vmr if  gas is provided in ty_gas_concs
      #
      if any(lower_case(this.gas_names[igas]) == gas_desc.gas_name[:])
         error_msg = gas_desc.get_vmr(this.gas_names[igas], vmr[:,:,igas])
         error_msg ≠ "" && return error_msg
      end
    end

    #
    # Compute dry air column amounts (number of molecule per cm^2) if user hasn't provided them
    #
    idx_h2o = string_loc_in_array("h2o", this.gas_names)
    if present(col_dry)
      col_dry_wk = col_dry
    else
      col_dry_arr = get_col_dry(vmr[:,:,idx_h2o], plev, tlay) # dry air column amounts computation
      col_dry_wk = col_dry_arr
    end
    #
    # compute column gas amounts [molec/cm^2]
    #
    col_gas[1:ncol,1:nlay,0] .= col_dry_wk[1:ncol,1:nlay]
    for igas = 1:ngas
      col_gas[1:ncol,1:nlay,igas] .= vmr[1:ncol,1:nlay,igas] .* col_dry_wk[1:ncol,1:nlay]
    end

    #
    # ---- calculate gas optical depths ----
    #
    #$acc enter data create(jtemp, jpress, jeta, tropo, fmajor)
    #$acc enter data create(tau, tau_rayleigh)
    #$acc enter data create(col_mix, fminor)
    #$acc enter data copyin(play, tlay, col_gas)
    #$acc enter data copyin(this)
    #$acc enter data copyin(this%gpoint_flavor)
    zero_array!(ngpt, nlay, ncol, tau)
    interpolation!(
            ncol,nlay,                        # problem dimensions
            ngas, nflav, neta, npres, ntemp,  # interpolation dimensions
            this.flavor,
            this.press_ref_log,
            this.temp_ref,
            this.press_ref_log_delta,
            this.temp_ref_min,
            this.temp_ref_delta,
            this.press_ref_trop_log,
            this.vmr_ref,
            play,
            tlay,
            col_gas,
            jtemp,         # outputs
            fmajor,fminor,
            col_mix,
            tropo,
            jeta,jpress)
    compute_tau_absorption!(
            ncol,nlay,nband,ngpt,                      # dimensions
            ngas,nflav,neta,npres,ntemp,
            nminorlower, nminorklower,                # number of minor contributors, total num absorption coeffs
            nminorupper, nminorkupper,
            idx_h2o,
            this.gpoint_flavor,
            get_band_lims_gpoint(this),
            this.kmajor,
            this.kminor_lower,
            this.kminor_upper,
            this.minor_limits_gpt_lower,
            this.minor_limits_gpt_upper,
            this.minor_scales_with_density_lower,
            this.minor_scales_with_density_upper,
            this.scale_by_complement_lower,
            this.scale_by_complement_upper,
            this.idx_minor_lower,
            this.idx_minor_upper,
            this.idx_minor_scaling_lower,
            this.idx_minor_scaling_upper,
            this.kminor_start_lower,
            this.kminor_start_upper,
            tropo,
            col_mix,fmajor,fminor,
            play,tlay,col_gas,
            jeta,jtemp,jpress,
            tau)
    if allocated(this.krayl)
      #$acc enter data attach(col_dry_wk) copyin(this%krayl)
      compute_tau_rayleigh!(          #Rayleigh scattering optical depths
            ncol,nlay,nband,ngpt,
            ngas,nflav,neta,npres,ntemp,  # dimensions
            this.gpoint_flavor,
            get_band_lims_gpoint(this),
            this.krayl,                   # inputs from object
            idx_h2o, col_dry_wk,col_gas,
            fminor,jeta,tropo,jtemp,      # local input
            tau_rayleigh)
      #$acc exit data detach(col_dry_wk) delete(this%krayl)
    end
    error_msg ≠ "" && return error_msg

    # Combine optical depths and reorder for radiative transfer solver.
    combine_and_reorder!(tau, tau_rayleigh, allocated(this.krayl), optical_props)
    #$acc exit data delete(tau, tau_rayleigh)
    #$acc exit data delete(play, tlay, col_gas)
    #$acc exit data delete(col_mix, fminor)
    #$acc exit data delete(this%gpoint_flavor)
    #$acc exit data copyout(jtemp, jpress, jeta, tropo, fmajor)
  end
  #------------------------------------------------------------------------------------------
  #
  # Compute Planck source functions at layer centers and levels
  #
  function source(this::ty_gas_optics_rrtmgp,
                  ncol, nlay, nbnd, ngpt,
                  play, plev, tlay, tsfc,
                  jtemp, jpress, jeta, tropo, fmajor,
                  sources::ty_source_func_lw,          # Planck sources
                  tlev)                                # optional input
                  #result(error_msg)
    # # inputs
    # class(ty_gas_optics_rrtmgp),    intent(in ) :: this
    # integer,                               intent(in   ) :: ncol, nlay, nbnd, ngpt
    # real(wp), dimension(ncol,nlay),        intent(in   ) :: play   # layer pressures [Pa, mb]
    # real(wp), dimension(ncol,nlay+1),      intent(in   ) :: plev   # level pressures [Pa, mb]
    # real(wp), dimension(ncol,nlay),        intent(in   ) :: tlay   # layer temperatures [K]
    # real(wp), dimension(ncol),             intent(in   ) :: tsfc   # surface skin temperatures [K]
    # # Interplation coefficients
    # integer,     dimension(ncol,nlay),     intent(in   ) :: jtemp, jpress
    # logical(wl), dimension(ncol,nlay),     intent(in   ) :: tropo
    # real(wp),    dimension(2,2,2,get_nflav(this),ncol,nlay),
    #                                        intent(in   ) :: fmajor
    # integer,     dimension(2,    get_nflav(this),ncol,nlay),
    #                                        intent(in   ) :: jeta
    # class(ty_source_func_lw    ),          intent(inout) :: sources
    # real(wp), dimension(ncol,nlay+1),      intent(in   ),
    #                                   optional, target :: tlev          # level temperatures [K]
    # character(len=128)                                 :: error_msg
    # ----------------------------------------------------------
    # integer                                      :: icol, ilay, igpt
    # real(wp), dimension(ngpt,nlay,ncol)          :: lay_source_t, lev_source_inc_t, lev_source_dec_t
    # real(wp), dimension(ngpt,     ncol)          :: sfc_source_t
    # # Variables for temperature at layer edges [K] (ncol, nlay+1)
    # real(wp), dimension(   ncol,nlay+1), target  :: tlev_arr
    # real(wp), dimension(:,:),            pointer :: tlev_wk

    DT
    lay_source_t = Array{DT}(undef, ngpt,nlay,ncol)
    lev_source_inc_t = Array{DT}(undef, ngpt,nlay,ncol)
    lev_source_dec_t = Array{DT}(undef, ngpt,nlay,ncol)
    sfc_source_t = Array{DT}(undef, ngpt,ncol)
    tlev_arr = Array{DT}(undef, ncol,nlay+1)

    # ----------------------------------------------------------
    error_msg = ""
    #
    # Source function needs temperature at interfaces/levels and at layer centers
    #
    if present(tlev)
      #   Users might have provided these
      tlev_wk = tlev
    else
      tlev_wk = tlev_arr
      #
      # Interpolate temperature to levels if not provided
      #   Interpolation and extrapolation at boundaries is weighted by pressure
      #
      for icol = 1:ncol
         tlev_arr[icol,1] = tlay[icol,1] +
                             (plev[icol,1]-play[icol,1])*(tlay[icol,2]-tlay[icol,1]) /
                                                         (play[icol,2]-play[icol,1])
      end
      for ilay in 2:nlay
        for icol in 1:ncol
           tlev_arr[icol,ilay] = (play[icol,ilay-1]*tlay[icol,ilay-1]*(plev[icol,ilay]-play[icol,ilay]) +
                                  play[icol,ilay]*tlay[icol,ilay]*(play[icol,ilay-1]-plev[icol,ilay])) /
                                  (plev[icol,ilay]*(play[icol,ilay-1] - play[icol,ilay]))
        end
      end
      for icol = 1:ncol
         tlev_arr[icol,nlay+1] = tlay[icol,nlay] +
                                 (plev[icol,nlay+1]-play[icol,nlay])*(tlay[icol,nlay]-tlay(icol,nlay-1)) /
                                                                       (play(icol,nlay)-play(icol,nlay-1))
      end
    end

    #-------------------------------------------------------------------
    # Compute internal (Planck) source functions at layers and levels,
    #  which depend on mapping from spectral space that creates k-distribution.
    #$acc enter data copyin(sources)
    #$acc enter data create(sources%lay_source, sources%lev_source_inc, sources%lev_source_dec, sources%sfc_source)
    #$acc enter data create(sfc_source_t, lay_source_t, lev_source_inc_t, lev_source_dec_t) attach(tlev_wk)
    compute_Planck_source!(ncol, nlay, nbnd, ngpt,
                get_nflav(this), get_neta(this), get_npres(this), get_ntemp(this), get_nPlanckTemp(this),
                tlay, tlev_wk, tsfc, fmerge(1,nlay,play[1,1] > play[1,nlay]),
                fmajor, jeta, tropo, jtemp, jpress,
                get_gpoint_bands(this), get_band_lims_gpoint(this), this.planck_frac, this.temp_ref_min,
                this.totplnk_delta, this.totplnk, this.gpoint_flavor,
                sfc_source_t, lay_source_t, lev_source_inc_t, lev_source_dec_t)
    #$acc parallel loop collapse(2)
    for igpt in 1:ngpt
      for icol in 1:ncol
        sources.sfc_source[icol,igpt] = sfc_source_t[igpt,icol]
      end
    end
    reorder123x321!(lay_source_t, sources.lay_source)
    reorder123x321!(lev_source_inc_t, sources.lev_source_inc)
    reorder123x321!(lev_source_dec_t, sources.lev_source_dec)
    #$acc exit data delete(sfc_source_t, lay_source_t, lev_source_inc_t, lev_source_dec_t) detach(tlev_wk)
    #$acc exit data copyout(sources%lay_source, sources%lev_source_inc, sources%lev_source_dec, sources%sfc_source)
    #$acc exit data copyout(sources)
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Initialization
  #
  #--------------------------------------------------------------------------------------------------------------------
  # Initialize object based on data read from netCDF file however the user desires.
  #  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  # This interface is for the internal-sources object -- includes Plank functions and fractions
  #
  function load!(this::ty_gas_optics_rrtmgp, available_gases::ty_gas_concs, gas_names, key_species,
                    band2gpt, band_lims_wavenum,
                    press_ref, press_ref_trop, temp_ref,
                    temp_ref_p, temp_ref_t, vmr_ref,
                    kmajor, kminor_lower, kminor_upper,
                    gas_minor,identifier_minor,
                    minor_gases_lower, minor_gases_upper,
                    minor_limits_gpt_lower, minor_limits_gpt_upper,
                    minor_scales_with_density_lower,
                    minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper,
                    scale_by_complement_lower,
                    scale_by_complement_upper,
                    kminor_start_lower,
                    kminor_start_upper,
                    totplnk, planck_frac, rayl_lower, rayl_upper) #result(err_message)
    # class(ty_gas_optics_rrtmgp),     intent(inout) :: this
    # class(ty_gas_concs),                    intent(in   ) :: available_gases # Which gases does the host model have available?
    # character(len=*),   dimension(:),       intent(in   ) :: gas_names
    # integer,            dimension(:,:,:),   intent(in   ) :: key_species
    # integer,            dimension(:,:),     intent(in   ) :: band2gpt
    # real(wp),           dimension(:,:),     intent(in   ) :: band_lims_wavenum
    # real(wp),           dimension(:),       intent(in   ) :: press_ref, temp_ref
    # real(wp),                               intent(in   ) :: press_ref_trop, temp_ref_p, temp_ref_t
    # real(wp),           dimension(:,:,:),   intent(in   ) :: vmr_ref
    # real(wp),           dimension(:,:,:,:), intent(in   ) :: kmajor
    # real(wp),           dimension(:,:,:),   intent(in   ) :: kminor_lower, kminor_upper
    # real(wp),           dimension(:,:),     intent(in   ) :: totplnk
    # real(wp),           dimension(:,:,:,:), intent(in   ) :: planck_frac
    # real(wp),           dimension(:,:,:),   intent(in   ),
    #                                           allocatable :: rayl_lower, rayl_upper
    # character(len=*),   dimension(:),       intent(in   ) :: gas_minor,identifier_minor
    # character(len=*),   dimension(:),       intent(in   ) :: minor_gases_lower,
    #                                                          minor_gases_upper
    # integer,            dimension(:,:),     intent(in   ) :: minor_limits_gpt_lower,
    #                                                          minor_limits_gpt_upper
    # logical(wl),        dimension(:),       intent(in   ) :: minor_scales_with_density_lower,
    #                                                          minor_scales_with_density_upper
    # character(len=*),   dimension(:),       intent(in   ) :: scaling_gas_lower,
    #                                                          scaling_gas_upper
    # logical(wl),        dimension(:),       intent(in   ) :: scale_by_complement_lower,
    #                                                          scale_by_complement_upper
    # integer,            dimension(:),       intent(in   ) :: kminor_start_lower,
    #                                                          kminor_start_upper
    # character(len = 128) :: err_message
    # # ----
    init_abs_coeffs!(this,
                     available_gases,
                     gas_names, key_species,
                     band2gpt, band_lims_wavenum,
                     press_ref, temp_ref,
                     press_ref_trop, temp_ref_p, temp_ref_t,
                     vmr_ref,
                     kmajor, kminor_lower, kminor_upper,
                     gas_minor,identifier_minor,
                     minor_gases_lower, minor_gases_upper,
                     minor_limits_gpt_lower,
                     minor_limits_gpt_upper,
                     minor_scales_with_density_lower,
                     minor_scales_with_density_upper,
                     scaling_gas_lower, scaling_gas_upper,
                     scale_by_complement_lower,
                     scale_by_complement_upper,
                     kminor_start_lower,
                     kminor_start_upper,
                     rayl_lower, rayl_upper)
    # Planck function tables
    #
    this.totplnk = totplnk
    this.planck_frac = planck_frac
    # Temperature steps for Planck function interpolation
    #   Assumes that temperature minimum and max are the same for the absorption coefficient grid and the
    #   Planck grid and the Planck grid is equally spaced
    this.totplnk_delta =  (this.temp_ref_max-this.temp_ref_min) / (size(this.totplnk, 1)-1)
  end

  #--------------------------------------------------------------------------------------------------------------------
  #
  # Initialize object based on data read from netCDF file however the user desires.
  #  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  # This interface is for the external-sources object -- includes TOA source function table
  #
  function load!(this::ty_gas_optics_rrtmgp, available_gases::ty_gas_concs, gas_names, key_species,
                    band2gpt, band_lims_wavenum,
                    press_ref, press_ref_trop, temp_ref,
                    temp_ref_p, temp_ref_t, vmr_ref,
                    kmajor, kminor_lower, kminor_upper,
                    gas_minor,identifier_minor,
                    minor_gases_lower, minor_gases_upper,
                    minor_limits_gpt_lower, minor_limits_gpt_upper,
                    minor_scales_with_density_lower,
                    minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper,
                    scale_by_complement_lower,
                    scale_by_complement_upper,
                    kminor_start_lower,
                    kminor_start_upper,
                    solar_src, rayl_lower, rayl_upper)
    # class(ty_gas_optics_rrtmgp), intent(inout) :: this
    # class(ty_gas_concs),                intent(in   ) :: available_gases # Which gases does the host model have available?
    # character(len=*),
    #           dimension(:),       intent(in) :: gas_names
    # integer,  dimension(:,:,:),   intent(in) :: key_species
    # integer,  dimension(:,:),     intent(in) :: band2gpt
    # real(wp), dimension(:,:),     intent(in) :: band_lims_wavenum
    # real(wp), dimension(:),       intent(in) :: press_ref, temp_ref
    # real(wp),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
    # real(wp), dimension(:,:,:),   intent(in) :: vmr_ref
    # real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    # real(wp), dimension(:,:,:),   intent(in) :: kminor_lower, kminor_upper
    # character(len=*),   dimension(:),
    #                               intent(in) :: gas_minor,
    #                                             identifier_minor
    # character(len=*),   dimension(:),
    #                               intent(in) :: minor_gases_lower,
    #                                             minor_gases_upper
    # integer,  dimension(:,:),     intent(in) ::
    #                                             minor_limits_gpt_lower,
    #                                             minor_limits_gpt_upper
    # logical(wl), dimension(:),    intent(in) ::
    #                                             minor_scales_with_density_lower,
    #                                             minor_scales_with_density_upper
    # character(len=*),   dimension(:),intent(in) ::
    #                                             scaling_gas_lower,
    #                                             scaling_gas_upper
    # logical(wl), dimension(:),    intent(in) ::
    #                                             scale_by_complement_lower,
    #                                             scale_by_complement_upper
    # integer,  dimension(:),       intent(in) ::
    #                                             kminor_start_lower,
    #                                             kminor_start_upper
    # real(wp), dimension(:),       intent(in), allocatable :: solar_src
    #                                                         # allocatable status to change when solar source is present in file
    # real(wp), dimension(:,:,:), intent(in), allocatable :: rayl_lower, rayl_upper
    # character(len = 128) err_message
    # ----
    init_abs_coeffs!(this,
                     available_gases,
                     gas_names, key_species,
                     band2gpt, band_lims_wavenum,
                     press_ref, temp_ref,
                     press_ref_trop, temp_ref_p, temp_ref_t,
                     vmr_ref,
                     kmajor, kminor_lower, kminor_upper,
                     gas_minor,identifier_minor,
                     minor_gases_lower, minor_gases_upper,
                     minor_limits_gpt_lower,
                     minor_limits_gpt_upper,
                     minor_scales_with_density_lower,
                     minor_scales_with_density_upper,
                     scaling_gas_lower, scaling_gas_upper,
                     scale_by_complement_lower,
                     scale_by_complement_upper,
                     kminor_start_lower,
                     kminor_start_upper,
                     rayl_lower, rayl_upper)
    #
    # Solar source table init
    #
    this.solar_src = solar_src

  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Initialize absorption coefficient arrays,
  #   including Rayleigh scattering tables if provided (allocated)
  #
  function init_abs_coeffs!(this::ty_gas_optics_rrtmgp,
                           available_gases::ty_gas_concs,
                           gas_names, key_species,
                           band2gpt, band_lims_wavenum,
                           press_ref, temp_ref,
                           press_ref_trop, temp_ref_p, temp_ref_t,
                           vmr_ref,
                           kmajor, kminor_lower, kminor_upper,
                           gas_minor,identifier_minor,
                           minor_gases_lower, minor_gases_upper,
                           minor_limits_gpt_lower,
                           minor_limits_gpt_upper,
                           minor_scales_with_density_lower,
                           minor_scales_with_density_upper,
                           scaling_gas_lower, scaling_gas_upper,
                           scale_by_complement_lower,
                           scale_by_complement_upper,
                           kminor_start_lower,
                           kminor_start_upper,
                           rayl_lower, rayl_upper)
    # class(ty_gas_optics_rrtmgp), intent(inout) :: this
    # class(ty_gas_concs),                intent(in   ) :: available_gases
    # character(len=*),
    #           dimension(:),       intent(in) :: gas_names
    # integer,  dimension(:,:,:),   intent(in) :: key_species
    # integer,  dimension(:,:),     intent(in) :: band2gpt
    # real(wp), dimension(:,:),     intent(in) :: band_lims_wavenum
    # real(wp), dimension(:),       intent(in) :: press_ref, temp_ref
    # real(wp),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
    # real(wp), dimension(:,:,:),   intent(in) :: vmr_ref
    # real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    # real(wp), dimension(:,:,:),   intent(in) :: kminor_lower, kminor_upper
    # character(len=*),   dimension(:),
    #                               intent(in) :: gas_minor,
    #                                             identifier_minor
    # character(len=*),   dimension(:),
    #                               intent(in) :: minor_gases_lower,
    #                                             minor_gases_upper
    # integer,  dimension(:,:),     intent(in) :: minor_limits_gpt_lower,
    #                                             minor_limits_gpt_upper
    # logical(wl), dimension(:),    intent(in) :: minor_scales_with_density_lower,
    #                                             minor_scales_with_density_upper
    # character(len=*),   dimension(:),
    #                               intent(in) :: scaling_gas_lower,
    #                                             scaling_gas_upper
    # logical(wl), dimension(:),    intent(in) :: scale_by_complement_lower,
    #                                             scale_by_complement_upper
    # integer,  dimension(:),       intent(in) :: kminor_start_lower,
    #                                             kminor_start_upper
    # real(wp), dimension(:,:,:),   intent(in),
    #                              allocatable :: rayl_lower, rayl_upper
    # character(len=128)                       :: err_message
    # --------------------------------------------------------------------------
    # logical,  dimension(:),     allocatable :: gas_is_present
    # logical,  dimension(:),     allocatable :: key_species_present_init
    # integer,  dimension(:,:,:), allocatable :: key_species_red
    # real(wp), dimension(:,:,:), allocatable :: vmr_ref_red
    # character(len=256),
    #           dimension(:),     allocatable :: minor_gases_lower_red,
    #                                            minor_gases_upper_red
    # character(len=256),
    #           dimension(:),     allocatable :: scaling_gas_lower_red,
    #                                            scaling_gas_upper_red
    # integer :: i, j, idx
    # integer :: ngas
    # # --------------------------------------
    DT = eltype(kmajor)

    init!(this, band_lims_wavenum, band2gpt)
    #
    # Which gases known to the gas optics are present in the host model (available_gases)?
    #
    ngas = size(gas_names)
    gas_is_present = Vector{Bool}(undef, ngas)

    for i in 1:ngas
      gas_is_present[i] = string_in_array(gas_names[i], available_gases.gas_name)
    end
    #
    # Now the number of gases is the union of those known to the k-distribution and provided
    #   by the host model
    #
    ngas = count(gas_is_present)
    #
    # Initialize the gas optics object, keeping only those gases known to the
    #   gas optics and also present in the host model
    #
    this.gas_names = pack(gas_names, mask=gas_is_present)

    vmr_ref_red = Array{DT}(undef, size(vmr_ref, 1),0:ngas, size(vmr_ref, 3))

    # Gas 0 is used in single-key species method, set to 1.0 (col_dry)
    vmr_ref_red[:,0,:] = vmr_ref[:,1,:]
    for i = 1:ngas
      idx = string_loc_in_array(this.gas_names[i], gas_names)
      vmr_ref_red[:,i,:] = vmr_ref[:,idx+1,:]
    end
    move_alloc!(vmr_ref_red, this.vmr_ref)
    #
    # Reduce minor arrays so variables only contain minor gases that are available
    # Reduce size of minor Arrays
    #
    reduce_minor_arrays!(available_gases,
                             gas_names,
                             gas_minor,identifier_minor,
                             kminor_lower,
                             minor_gases_lower,
                             minor_limits_gpt_lower,
                             minor_scales_with_density_lower,
                             scaling_gas_lower,
                             scale_by_complement_lower,
                             kminor_start_lower,
                             this.kminor_lower,
                             minor_gases_lower_red,
                             this.minor_limits_gpt_lower,
                             this.minor_scales_with_density_lower,
                             scaling_gas_lower_red,
                             this.scale_by_complement_lower,
                             this.kminor_start_lower)
    reduce_minor_arrays!(available_gases,
                             gas_names,
                             gas_minor,identifier_minor,
                             kminor_upper,
                             minor_gases_upper,
                             minor_limits_gpt_upper,
                             minor_scales_with_density_upper,
                             scaling_gas_upper,
                             scale_by_complement_upper,
                             kminor_start_upper,
                             this.kminor_upper,
                             minor_gases_upper_red,
                             this.minor_limits_gpt_upper,
                             this.minor_scales_with_density_upper,
                             scaling_gas_upper_red,
                             this.scale_by_complement_upper,
                             this.kminor_start_upper)

    # Arrays not reduced by the presence, or lack thereof, of a gas
    this.press_ref = press_ref
    this.temp_ref  = temp_ref
    this.kmajor    = kmajor
    DT = eltype(kmajor)

    # TODO: Check if .neqv. is the same as ≠
    if allocated(rayl_lower) ≠ allocated(rayl_upper)
      error("rayl_lower and rayl_upper must have the same allocation status")
    end
    if allocated(rayl_lower)
      this.krayl = Array{DT}(undef, size(rayl_lower,1),size(rayl_lower,2),size(rayl_lower,3),2)
      this.krayl[:,:,:,1] = rayl_lower
      this.krayl[:,:,:,2] = rayl_upper
    end

    # ---- post processing ----
    # Incoming coefficients file has units of Pa
    this.press_ref[:] .= this.press_ref[:]

    # creates log reference pressure
    this.press_ref_log = Array{DT}(undef, size(this.press_ref))
    this.press_ref_log[:] = log(this.press_ref[:])

    # log scale of reference pressure
    this.press_ref_trop_log = log(press_ref_trop)

    # Get index of gas (if present) for determining col_gas
    create_idx_minor!(this.gas_names, gas_minor, identifier_minor, minor_gases_lower_red, this.idx_minor_lower)
    create_idx_minor!(this.gas_names, gas_minor, identifier_minor, minor_gases_upper_red, this.idx_minor_upper)
    # Get index of gas (if present) that has special treatment in density scaling
    create_idx_minor_scaling!(this.gas_names, scaling_gas_lower_red, this.idx_minor_scaling_lower)
    create_idx_minor_scaling!(this.gas_names, scaling_gas_upper_red, this.idx_minor_scaling_upper)

    # create flavor list
    # Reduce (remap) key_species list; checks that all key gases are present in incoming
    create_key_species_reduce!(gas_names,this.gas_names, key_species,key_species_red,key_species_present_init)
    check_key_species_present_init(gas_names,key_species_present_init)

    # create flavor list
    create_flavor!(key_species_red, this.flavor)
    # create gpoint_flavor list
    create_gpoint_flavor!(key_species_red, get_gpoint_bands(this), this.flavor, this.gpoint_flavor)

    # minimum, maximum reference temperature, pressure -- assumes low-to-high ordering
    #   for T, high-to-low ordering for p
    this.temp_ref_min  = this.temp_ref[1]
    this.temp_ref_max  = this.temp_ref[length(this.temp_ref)]
    this.press_ref_min = this.press_ref[length(this.press_ref)]
    this.press_ref_max = this.press_ref[1]

    # creates press_ref_log, temp_ref_delta
    this.press_ref_log_delta = (log(this.press_ref_min)-log(this.press_ref_max))/(length(this.press_ref)-1)
    this.temp_ref_delta      = (this.temp_ref_max-this.temp_ref_min)/(length(this.temp_ref)-1)

    # Which species are key in one or more bands?
    #   this%flavor is an index into this%gas_names
    #
    allocated(this.is_key) && deallocate!(this.is_key) # Shouldn't ever happen...
    this.is_key = Array{Bool}(undef, get_ngas(this))
    this.is_key[:] = false
    for j in 1:size(this.flavor, 2)
      for i in 1:size(this.flavor, 1) # should be 2
        if this.flavor[i,j] ≠ 0
          this.is_key[this.flavor[i,j]] = true
        end
      end
    end

  end
  # ----------------------------------------------------------------------------------------------------
  function check_key_species_present_init(gas_names, key_species_present_init) # result(err_message)
    # logical,          dimension(:), intent(in) :: key_species_present_init
    # character(len=*), dimension(:), intent(in) :: gas_names
    # character(len=128)                             :: err_message

    # integer :: i

    for i in 1:length(key_species_present_init)
      if !key_species_present_init[i]
        error("gas_optics: required gases" * trim(gas_names[i]) * " are not provided")
      end
    end
  end
  #------------------------------------------------------------------------------------------
  #
  # Ensure that every key gas required by the k-distribution is
  #    present in the gas concentration object
  #
  function check_key_species_present(this::ty_gas_optics_rrtmgp, gas_desc::ty_gas_concs) # result(error_msg)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # class(ty_gas_concs),                intent(in) :: gas_desc
    # character(len=128)                             :: error_msg

    # # Local variables
    # character(len=32), dimension(count(this%is_key(:)  )) :: key_gas_names
    # integer                                               :: igas
    # # --------------------------------------
    error_msg = ""
    key_gas_names = pack(this.gas_names, mask=this.is_key)
    for igas = 1:length(key_gas_names)
      if !string_in_array(key_gas_names[igas], gas_desc.gas_name)
        error_msg = " " * trim(lower_case(key_gas_names[igas])) * trim(error_msg)
      end
    end
    if len_trim(error_msg) > 0
      error_msg = "gas_optics: required gases" * trim(error_msg) * " are not provided"
    end
    return error_msg
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Function to define names of key and minor gases to be used by gas_optics().
  # The final list gases includes those that are defined in gas_optics_specification
  # and are provided in ty_gas_concs.
  #
  function get_minor_list(this::ty_gas_optics_rrtmgp, gas_desc::ty_gas_concs, ngas, names_spec)
    # class(ty_gas_optics_rrtmgp), intent(in)       :: this
    # class(ty_gas_concs), intent(in)                      :: gas_desc
    # integer, intent(in)                                  :: ngas
    # character(32), dimension(ngas), intent(in)           :: names_spec

    # # List of minor gases to be used in gas_optics()
    # character(len=32), dimension(:), allocatable         :: get_minor_list
    # # Logical flag for minor species in specification (T = minor; F = not minor)
    # logical, dimension(size(names_spec))                 :: gas_is_present
    # integer                                              :: igas, icnt

    allocated(get_minor_list) && deallocate!(get_minor_list)
    for igas = 1:get_ngas(this)
      gas_is_present[igas] = string_in_array(names_spec[igas], gas_desc.gas_name)
    end
    icnt = count(gas_is_present)
    get_minor_list = Vector{String}(undef, icnt)
    get_minor_list[:] .= pack(this.gas_names, mask=gas_is_present)
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Inquiry functions
  #
  #--------------------------------------------------------------------------------------------------------------------
  #
  # return true if initialized for internal sources, false otherwise
  #
  function source_is_internal(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # logical                          :: source_is_internal
    return allocated(this.totplnk) && allocated(this.planck_frac)
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # return true if initialized for external sources, false otherwise
  #
  function source_is_external(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # logical                          :: source_is_external
    return allocated(this.solar_src)
  end

  #--------------------------------------------------------------------------------------------------------------------
  #
  # return the gas names
  #
  function get_gases(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # character(32), dimension(get_ngas(this))     :: get_gases

    return this.gas_names
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # return the minimum pressure on the interpolation grids
  #
  function get_press_min(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # real(wp)                                       :: get_press_min

    return this.press_ref_min
  end

  #--------------------------------------------------------------------------------------------------------------------
  #
  # return the maximum pressure on the interpolation grids
  #
  function get_press_max(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # real(wp)                                       :: get_press_max

    return this.press_ref_max
  end

  #--------------------------------------------------------------------------------------------------------------------
  #
  # return the minimum temparature on the interpolation grids
  #
  function get_temp_min(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # real(wp)                                       :: get_temp_min

    return this.temp_ref_min
  end

  #--------------------------------------------------------------------------------------------------------------------
  #
  # return the maximum temparature on the interpolation grids
  #
  function get_temp_max(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # real(wp)                                       :: get_temp_max

    return this.temp_ref_max
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Utility function, provided for user convenience
  # computes column amounts of dry air using hydrostatic equation
  #
  function get_col_dry(vmr_h2o, plev, tlay, latitude) # result(col_dry)
    # # input
    # real(wp), dimension(:,:), intent(in) :: vmr_h2o  # volume mixing ratio of water vapor to dry air; (ncol,nlay)
    # real(wp), dimension(:,:), intent(in) :: plev     # Layer boundary pressures [Pa] (ncol,nlay+1)
    # real(wp), dimension(:,:), intent(in) :: tlay     # Layer temperatures [K] (ncol,nlay)
    # real(wp), dimension(:),   optional,
    #                           intent(in) :: latitude # Latitude [degrees] (ncol)
    # # output
    # real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: col_dry # Column dry amount (ncol,nlay)
    # ------------------------------------------------
    # # first and second term of Helmert formula
    # real(wp), parameter :: helmert1 = 9.80665_wp
    # real(wp), parameter :: helmert2 = 0.02586_wp
    # # local variables
    # real(wp), dimension(size(tlay,dim=1)                 ) :: g0 # (ncol)
    # real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: delta_plev # (ncol,nlay)
    # real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: m_air # average mass of air; (ncol,nlay)
    # integer :: nlev, nlay
    # # ------------------------------------------------
    DT = eltype(plev)

    # first and second term of Helmert formula
    helmert1 = DT(9.80665)
    helmert2 = DT(0.02586)
    # local variables
    g0         = Array{DT}(undef, size(tlay,1)             ) # (ncol)
    delta_plev = Array{DT}(undef, size(tlay,1),size(tlay,2)) # (ncol,nlay)
    m_air      = Array{DT}(undef, size(tlay,1),size(tlay,2)) # average mass of air; (ncol,nlay)
    integer :: nlev, nlay
    # ------------------------------------------------
    nlay = size(tlay, 2)
    nlev = size(plev, 2)

    if present(latitude)
      g0[:] .= helmert1 - helmert2 * cos(DT(2) * π * latitude[:] / DT(180)) # acceleration due to gravity [m/s^2]
    else
      g0[:] = grav
    end
    delta_plev[:,:] .= abs(plev[:,1:nlev-1] - plev[:,2:nlev])

    # Get average mass of moist air per mole of moist air
    m_air[:,:] .= (m_dry .+ m_h2o .* vmr_h2o[:,:]) ./ (1 .+ vmr_h2o[:,:])

    # Hydrostatic equation
    col_dry[:,:] = DT(10) .* delta_plev[:,:] .* avogad ./ (DT(1000)*m_air[:,:]*DT(100)*spread(g0[:],dim=2, ncopies=nlay))
    col_dry[:,:] = col_dry[:,:] ./ (DT(1) .+ vmr_h2o[:,:])
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Internal procedures
  #
  #--------------------------------------------------------------------------------------------------------------------
  function rewrite_key_species_pair(key_species_pair)
    # (0,0) becomes (2,2) -- because absorption coefficients for these g-points will be 0.
    # integer, dimension(2) :: rewrite_key_species_pair
    # integer, dimension(2), intent(in) :: key_species_pair
    result = key_species_pair
    if all(key_species_pair[:] .== [0,0])
      result[:] .= [2,2]
    end
  end

  # ---------------------------------------------------------------------------------------
  # true is key_species_pair exists in key_species_list
  function key_species_pair_exists(key_species_list, key_species_pair)
    # logical                             :: key_species_pair_exists
    # integer, dimension(:,:), intent(in) :: key_species_list
    # integer, dimension(2),   intent(in) :: key_species_pair
    # integer :: i
    for i=1:size(key_species_list,2)
      if all(key_species_list[:,i] .== key_species_pair[:])
        result = true
        return
      end
    end
    result = false
  end
  # ---------------------------------------------------------------------------------------
  # create flavor list --
  #   an unordered array of extent (2,:) containing all possible pairs of key species
  #   used in either upper or lower atmos
  #
  function create_flavor(key_species, flavor)
    # integer, dimension(:,:,:), intent(in) :: key_species
    # integer, dimension(:,:), allocatable, intent(out) :: flavor
    # integer, dimension(2,size(key_species,3)*2) :: key_species_list
    # integer :: ibnd, iatm, i, iflavor

    key_species_list = Array{Int}(undef, 2,size(key_species,3)*2)

    # prepare list of key_species
    i = 1
    for ibnd=1:size(key_species,3)
      for iatm=1:size(key_species,1)
        key_species_list[:,i] .= key_species[:,iatm,ibnd]
        i = i + 1
      end
    end
    # rewrite single key_species pairs
    for i=1:size(key_species_list,2)
        key_species_list[:,i] = rewrite_key_species_pair(key_species_list[:,i])
    end
    # count unique key species pairs
    iflavor = 0
    for i=1:size(key_species_list,2)
      if !key_species_pair_exists(key_species_list[:,1:i-1],key_species_list[:,i])
        iflavor = iflavor + 1
      end
    end
    # fill flavors
    flavor = Array{Int}(undef, 2,iflavor)
    iflavor = 0
    for i=1:size(key_species_list,2)
      if !key_species_pair_exists(key_species_list[:,1:i-1],key_species_list[:,i])
        iflavor = iflavor + 1
        flavor[:,iflavor] = key_species_list[:,i]
      end
    end
  end
  # ---------------------------------------------------------------------------------------
  #
  # create index list for extracting col_gas needed for minor gas optical depth calculations
  #
  function create_idx_minor!(gas_names, gas_minor, identifier_minor, minor_gases_atm, idx_minor_atm)
    # character(len=*), dimension(:), intent(in) :: gas_names
    # character(len=*), dimension(:), intent(in) ::
    #                                               gas_minor,
    #                                               identifier_minor
    # character(len=*), dimension(:), intent(in) :: minor_gases_atm
    # integer, dimension(:), allocatable,
    #                                intent(out) :: idx_minor_atm

    # # local
    # integer :: imnr
    # integer :: idx_mnr
    idx_minor_atm = Vector{Int}(undef, size(minor_gases_atm,1))
    for imnr = 1:size(minor_gases_atm,1) # loop over minor absorbers in each band
          # Find identifying string for minor species in list of possible identifiers (e.g. h2o_slf)
          idx_mnr     = string_loc_in_array(minor_gases_atm[imnr], identifier_minor)
          # Find name of gas associated with minor species identifier (e.g. h2o)
          idx_minor_atm[imnr] = string_loc_in_array(gas_minor[idx_mnr],    gas_names)
    end

  end

  # ---------------------------------------------------------------------------------------
  #
  # create index for special treatment in density scaling of minor gases
  #
  function create_idx_minor_scaling!(gas_names,
    scaling_gas_atm, idx_minor_scaling_atm)
    # character(len=*), dimension(:), intent(in) :: gas_names
    # character(len=*), dimension(:), intent(in) :: scaling_gas_atm
    # integer, dimension(:), allocatable,
    #                                intent(out) :: idx_minor_scaling_atm

    # # local
    # integer :: imnr
    idx_minor_scaling_atm = Vector{Int}(undef, size(scaling_gas_atm,1))
    for imnr = 1:size(scaling_gas_atm,1) # loop over minor absorbers in each band
          # This will be -1 if there's no interacting gas
          idx_minor_scaling_atm[imnr] = string_loc_in_array(scaling_gas_atm[imnr], gas_names)
    end

  end
  # ---------------------------------------------------------------------------------------
  function create_key_species_reduce!(gas_names, gas_names_red,
    key_species,key_species_red,key_species_present_init)
    # character(len=*),
    #           dimension(:),       intent(in) :: gas_names
    # character(len=*),
    #           dimension(:),       intent(in) :: gas_names_red
    # integer,  dimension(:,:,:),   intent(in) :: key_species
    # integer,  dimension(:,:,:), allocatable, intent(out) :: key_species_red

    # logical, dimension(:), allocatable, intent(out) :: key_species_present_init
    # integer :: ip, ia, it, np, na, nt

    np = size(key_species,1)
    na = size(key_species,2)
    nt = size(key_species,3)
    key_species_red = Array{Int}(undef, size(key_species,1),
                                        size(key_species,2),
                                        size(key_species,3))
    key_species_present_init = Vector{Bool}(undef, size(gas_names))
    key_species_present_init = true

    for ip = 1:np
      for ia = 1:na
        for it = 1:nt
          if key_species[ip,ia,it] ≠ 0
            key_species_red[ip,ia,it] = string_loc_in_array(gas_names[key_species[ip,ia,it]],gas_names_red)
            if key_species_red[ip,ia,it] == -1
              key_species_present_init[key_species[ip,ia,it]] = false
            end
          else
            key_species_red[ip,ia,it] = key_species[ip,ia,it]
          end
        end
      end
    end

  end

# ---------------------------------------------------------------------------------------
  function reduce_minor_arrays(available_gases::ty_gas_concs,
                           gas_names,
                           gas_minor,identifier_minor,
                           kminor_atm,
                           minor_gases_atm,
                           minor_limits_gpt_atm,
                           minor_scales_with_density_atm,
                           scaling_gas_atm,
                           scale_by_complement_atm,
                           kminor_start_atm,
                           kminor_atm_red,
                           minor_gases_atm_red,
                           minor_limits_gpt_atm_red,
                           minor_scales_with_density_atm_red,
                           scaling_gas_atm_red,
                           scale_by_complement_atm_red,
                           kminor_start_atm_red)

    # class(ty_gas_concs),                intent(in   ) :: available_gases
    # character(len=*), dimension(:),     intent(in) :: gas_names
    # real(wp),         dimension(:,:,:), intent(in) :: kminor_atm
    # character(len=*), dimension(:),     intent(in) :: gas_minor,
    #                                                   identifier_minor
    # character(len=*), dimension(:),     intent(in) :: minor_gases_atm
    # integer,          dimension(:,:),   intent(in) :: minor_limits_gpt_atm
    # logical(wl),      dimension(:),     intent(in) :: minor_scales_with_density_atm
    # character(len=*), dimension(:),     intent(in) :: scaling_gas_atm
    # logical(wl),      dimension(:),     intent(in) :: scale_by_complement_atm
    # integer,          dimension(:),     intent(in) :: kminor_start_atm
    # real(wp),         dimension(:,:,:), allocatable,
    #                                     intent(out) :: kminor_atm_red
    # character(len=*), dimension(:), allocatable,
    #                                     intent(out) :: minor_gases_atm_red
    # integer,          dimension(:,:), allocatable,
    #                                     intent(out) :: minor_limits_gpt_atm_red
    # logical(wl),      dimension(:),    allocatable,
    #                                     intent(out) ::minor_scales_with_density_atm_red
    # character(len=*), dimension(:), allocatable,
    #                                     intent(out) ::scaling_gas_atm_red
    # logical(wl),      dimension(:), allocatable, intent(out) ::
    #                                             scale_by_complement_atm_red
    # integer,          dimension(:), allocatable, intent(out) ::
    #                                             kminor_start_atm_red

    # # Local variables
    # integer :: i, j
    # integer :: idx_mnr, nm, tot_g, red_nm
    # integer :: icnt, n_elim, ng
    # logical, dimension(:), allocatable :: gas_is_present
    DT = eltype(kminor_atm)

    nm = size(minor_gases_atm)
    tot_g=0
    gas_is_present = Vector{Bool}(undef, nm)
    for i = 1:length(minor_gases_atm)
      idx_mnr = string_loc_in_array(minor_gases_atm[i], identifier_minor)
      gas_is_present[i] = string_in_array(gas_minor[idx_mnr],available_gases.gas_name)
      if gas_is_present[i]
        tot_g = tot_g + (minor_limits_gpt_atm[2,i]-minor_limits_gpt_atm[1,i]+1)
      end
    end
    red_nm = count(gas_is_present)

    if red_nm == nm
      kminor_atm_red = kminor_atm
      minor_gases_atm_red = minor_gases_atm
      minor_limits_gpt_atm_red = minor_limits_gpt_atm
      minor_scales_with_density_atm_red = minor_scales_with_density_atm
      scaling_gas_atm_red = scaling_gas_atm
      scale_by_complement_atm_red = scale_by_complement_atm
      kminor_start_atm_red = kminor_start_atm
    else
      minor_gases_atm_red= pack(minor_gases_atm, mask=gas_is_present)
      minor_scales_with_density_atm_red = pack(minor_scales_with_density_atm,
        mask=gas_is_present)
      scaling_gas_atm_red = pack(scaling_gas_atm,
        mask=gas_is_present)
      scale_by_complement_atm_red = pack(scale_by_complement_atm,
        mask=gas_is_present)
      kminor_start_atm_red = pack(kminor_start_atm,
        mask=gas_is_present)

      minor_limits_gpt_atm_red = Array{Int}(undef, 2, red_nm)
      kminor_atm_red = Array{DT}(undef, tot_g, size(kminor_atm,2), size(kminor_atm,3))

      icnt = 0
      n_elim = 0
      for i = 1:nm
        ng = minor_limits_gpt_atm[2,i]-minor_limits_gpt_atm[1,i]+1
        if gas_is_present[i]
          icnt = icnt + 1
          minor_limits_gpt_atm_red[1:2,icnt] = minor_limits_gpt_atm[1:2,i]
          kminor_start_atm_red[icnt] = kminor_start_atm[i]-n_elim
          for j = 1:ng
            kminor_atm_red[kminor_start_atm_red[icnt]+j-1,:,:] = kminor_atm[kminor_start_atm[i]+j-1,:,:]
          end
        else
          n_elim = n_elim + ng
        end
      end
    end

  end

# ---------------------------------------------------------------------------------------
  # returns flavor index; -1 if not found
  function key_species_pair2flavor(flavor, key_species_pair)
    # integer :: key_species_pair2flavor
    # integer, dimension(:,:), intent(in) :: flavor
    # integer, dimension(2), intent(in) :: key_species_pair
    # integer :: iflav
    for iflav=1:size(flavor,2)
      if all(key_species_pair[:] == flavor[:,iflav])
        return iflav
      end
    end
    return -1
  end

  # ---------------------------------------------------------------------------------------
  #
  # create gpoint_flavor list
  #   a map pointing from each g-point to the corresponding entry in the "flavor list"
  #
  function create_gpoint_flavor(key_species, gpt2band, flavor, gpoint_flavor)
    # integer, dimension(:,:,:), intent(in) :: key_species
    # integer, dimension(:), intent(in) :: gpt2band
    # integer, dimension(:,:), intent(in) :: flavor
    # integer, dimension(:,:), intent(out), allocatable :: gpoint_flavor
    # integer :: ngpt, igpt, iatm
    ngpt = size(gpt2band)
    gpoint_flavor = Array{Int}(undef, 2,ngpt)
    for igpt=1:ngpt
      for iatm=1:2
        gpoint_flavor[iatm,igpt] = key_species_pair2flavor(
          flavor,
          rewrite_key_species_pair(key_species[:,iatm,gpt2band[igpt]])
        )
      end
    end
  end

 #--------------------------------------------------------------------------------------------------------------------
 #
 # Utility function to combine optical depths from gas absorption and Rayleigh scattering
 #   (and reorder them for convenience, while we're at it)
 #
 function combine_and_reorder(tau, tau_rayleigh, has_rayleigh, optical_props::ty_optical_props_arry)
    # real(wp), dimension(:,:,:),   intent(in) :: tau
    # real(wp), dimension(:,:,:),   intent(in) :: tau_rayleigh
    # logical,                      intent(in) :: has_rayleigh
    # class(ty_optical_props_arry), intent(inout) :: optical_props

    # integer :: ncol, nlay, ngpt, nmom

    ncol = size(tau, 3)
    nlay = size(tau, 2)
    ngpt = size(tau, 1)
    #$acc enter data copyin(optical_props)
    if !has_rayleigh
      # index reorder (ngpt, nlay, ncol) -> (ncol,nlay,gpt)
      #$acc enter data copyin(tau)
      #$acc enter data create(optical_props%tau)
      reorder123x321!(tau, optical_props.tau)
        if optical_props isa ty_optical_props_2str
          #$acc enter data create(optical_props%ssa, optical_props%g)
          zero_array!(     ncol,nlay,ngpt,optical_props.ssa)
          zero_array!(     ncol,nlay,ngpt,optical_props.g  )
          #$acc exit data copyout(optical_props%ssa, optical_props%g)
        elseif optical_props isa ty_optical_props_nstr # We ought to be able to combine this with above
          nmom = size(optical_props.p, 1)
          #$acc enter data create(optical_props%ssa, optical_props%p)
          zero_array!(     ncol,nlay,ngpt,optical_props.ssa)
          zero_array!(nmom,ncol,nlay,ngpt,optical_props.p  )
          #$acc exit data copyout(optical_props%ssa, optical_props%p)
        end
      #$acc exit data copyout(optical_props%tau)
      #$acc exit data delete(tau)
    else
      # combine optical depth and rayleigh scattering
      #$acc enter data copyin(tau, tau_rayleigh)
        if optical_props isa ty_optical_props_1scl
          # User is asking for absorption optical depth
          #$acc enter data create(optical_props%tau)
          reorder123x321!(tau, optical_props.tau)
          #$acc exit data copyout(optical_props%tau)
        elseif optical_props isa ty_optical_props_2str
          #$acc enter data create(optical_props%tau, optical_props%ssa, optical_props%g)
          combine_and_reorder_2str!(ncol, nlay, ngpt,       tau, tau_rayleigh,
                                        optical_props.tau, optical_props.ssa, optical_props.g)
          #$acc exit data copyout(optical_props%tau, optical_props%ssa, optical_props%g)
        elseif optical_props isa ty_optical_props_nstr # We ought to be able to combine this with above
          nmom = size(optical_props.p, 1)
          #$acc enter data create(optical_props%tau, optical_props%ssa, optical_props%p)
          combine_and_reorder_nstr!(ncol, nlay, ngpt, nmom, tau, tau_rayleigh,
                                        optical_props.tau, optical_props.ssa, optical_props.p)
          #$acc exit data copyout(optical_props%tau, optical_props%ssa, optical_props%p)
        end
      #$acc exit data delete(tau, tau_rayleigh)
    end
    #$acc exit data copyout(optical_props)
  end

  #--------------------------------------------------------------------------------------------------------------------
  # Sizes of tables: pressure, temperate, eta (mixing fraction)
  #   Equivalent routines for the number of gases and flavors (get_ngas(), get_nflav()) are defined above because they're
  #   used in function defintions
  # Table kmajor has dimensions (ngpt, neta, npres, ntemp)
  #--------------------------------------------------------------------------------------------------------------------
  #
  # return extent of eta dimension
  #
  function get_neta(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # integer                          :: get_neta

    return size(this.kmajor,2)
  end
  # --------------------------------------------------------------------------------------
  #
  # return the number of pressures in reference profile
  #   absorption coefficient table is one bigger since a pressure is repeated in upper/lower atmos
  #
  function get_npres(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # integer                          :: get_npres

    return size(this.kmajor,3)-1
  end
  # --------------------------------------------------------------------------------------
  #
  # return the number of temperatures
  #
  function get_ntemp(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # integer                          :: get_ntemp

    return size(this.kmajor,4)
  end
  # --------------------------------------------------------------------------------------
  #
  # return the number of temperatures for Planck function
  #
  function get_nPlanckTemp(this::ty_gas_optics_rrtmgp)
    # class(ty_gas_optics_rrtmgp), intent(in) :: this
    # integer                          :: get_nPlanckTemp

    return size(this.totplnk,1) # dimensions are Planck-temperature, band
  end
  #--------------------------------------------------------------------------------------------------------------------
  # Generic procedures for checking sizes, limits
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Extents
  #
  # --------------------------------------------------------------------------------------
  function check_extent_1d(array, n1, label)
    # real(wp), dimension(:          ), intent(in) :: array
    # integer,                          intent(in) :: n1
    # character(len=*),                 intent(in) :: label
    # character(len=128)                           :: check_extent_1d

    s = ""
    if size(array,1) ≠ n1
      s = trim(label) * " has incorrect size."
    end
    return s
  end
  # --------------------------------------------------------------------------------------
  function check_extent_2d(array, n1, n2, label)
    # real(wp), dimension(:,:        ), intent(in) :: array
    # integer,                          intent(in) :: n1, n2
    # character(len=*),                 intent(in) :: label
    # character(len=128)                           :: check_extent_2d

    s = ""
    if(size(array,1) ≠ n1 || size(array,2) ≠ n2 )
      s = trim(label) * " has incorrect size."
    end
    return s
  end
  # --------------------------------------------------------------------------------------
  function check_extent_3d(array, n1, n2, n3, label)
    # real(wp), dimension(:,:,:      ), intent(in) :: array
    # integer,                          intent(in) :: n1, n2, n3
    # character(len=*),                 intent(in) :: label
    # character(len=128)                           :: check_extent_3d

    s = ""
    if(size(array,1) ≠ n1 || size(array,2) ≠ n2 || size(array,3) ≠ n3)
      s = trim(label) * " has incorrect size."
    end
    return s
  end
  # --------------------------------------------------------------------------------------
  function check_extent_4d(array, n1, n2, n3, n4, label)
    # real(wp), dimension(:,:,:,:    ), intent(in) :: array
    # integer,                          intent(in) :: n1, n2, n3, n4
    # character(len=*),                 intent(in) :: label
    # character(len=128)                           :: check_extent_4d

    s = ""
    if(size(array,1) ≠ n1 || size(array,2) ≠ n2 || size(array,3) ≠ n3 ||
       size(array,4) ≠ n4)
      s = trim(label) * " has incorrect size."
    end
    return s
  end
  # --------------------------------------------------------------------------------------
  function check_extent_5d(array, n1, n2, n3, n4, n5, label)
    # real(wp), dimension(:,:,:,:,:  ), intent(in) :: array
    # integer,                          intent(in) :: n1, n2, n3, n4, n5
    # character(len=*),                 intent(in) :: label
    # character(len=128)                           :: check_extent_5d

    s = ""
    if(size(array,1) ≠ n1 || size(array,2) ≠ n2 || size(array,3) ≠ n3 ||
       size(array,4) ≠ n4 || size(array,5) ≠ n5)
      s = trim(label) * " has incorrect size."
    end
    return s
  end
  # --------------------------------------------------------------------------------------
  function check_extent_6d(array, n1, n2, n3, n4, n5, n6, label)
    # real(wp), dimension(:,:,:,:,:,:), intent(in) :: array
    # integer,                          intent(in) :: n1, n2, n3, n4, n5, n6
    # character(len=*),                 intent(in) :: label
    # character(len=128)                           :: check_extent_6d

    s = ""
    if(size(array,1) ≠ n1 || size(array,2) ≠ n2 || size(array,3) ≠ n3 ||
       size(array,4) ≠ n4 || size(array,5) ≠ n5 || size(array,6) ≠ n6 )
      s = trim(label) * " has incorrect size."
    end
    return s
  end
  # --------------------------------------------------------------------------------------
  #
  # Values
  #
  # --------------------------------------------------------------------------------------
  function check_range_1D(val, minV, maxV, label)
    # real(wp), dimension(:),     intent(in) :: val
    # real(wp),                   intent(in) :: minV, maxV
    # character(len=*),           intent(in) :: label
    # character(len=128)                     :: check_range_1D

    s = ""
    if(any(val < minV) || any(val > maxV))
      s = trim(label) * " values out of range."
    end
    return s
  end
  # --------------------------------------------------------------------------------------
  function check_range_2D(val, minV, maxV, label)
    # real(wp), dimension(:,:),   intent(in) :: val
    # real(wp),                   intent(in) :: minV, maxV
    # character(len=*),           intent(in) :: label
    # character(len=128)                     :: check_range_2D

    s = ""
    if(any(val < minV) || any(val > maxV))
      s = trim(label) * " values out of range."
    end
    return s
  end
  # --------------------------------------------------------------------------------------
  function check_range_3D(val, minV, maxV, label)
    # real(wp), dimension(:,:,:), intent(in) :: val
    # real(wp),                   intent(in) :: minV, maxV
    # character(len=*),           intent(in) :: label
    # character(len=128)                     :: check_range_3D

    s = ""
    if(any(val < minV) || any(val > maxV))
      s = trim(label) * " values out of range."
    end
    return s
  end
  #------------------------------------------------------------------------------------------
end
