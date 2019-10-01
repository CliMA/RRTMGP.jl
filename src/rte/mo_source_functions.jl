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
# Encapsulate source function arrays for longwave/lw/internal sources
#    and shortwave/sw/external source.
#
# -------------------------------------------------------------------------------------------------
module mo_source_functions

  using ..mo_optical_props

  # -------------------------------------------------------------------------------------------------
  #
  # Type for longwave sources: computed at layer center, at layer edges using
  #   spectral mapping in each direction separately, and at the surface


  export ty_source_func_lw
  struct ty_source_func_lw{T, I} <: ty_optical_props{T, I}
    band2gpt::Array{T,2}        # (begin g-point, end g-point) = band2gpt(2,band)
    gpt2band::Array{I,1}        # band = gpt2band(g-point)
    band_lims_wvn::Array{T,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
    name::String
    tau::Array{T,3}
    #
    lay_source     # Planck source at layer average temperature [W/m2] (ncol, nlay, ngpt)
    lev_source_inc # Planck source at layer edge in increasing ilay direction [W/m2] (ncol, nlay+1, ngpt)
    lev_source_dec # Planck source at layer edge in decreasing ilay direction [W/m2] (ncol, nlay+1, ngpt)
    sfc_source
  end
  export ty_source_func_sw
  struct ty_source_func_sw{T, I} <: ty_optical_props{T, I}
    band2gpt::Array{T,2}        # (begin g-point, end g-point) = band2gpt(2,band)
    gpt2band::Array{I,1}        # band = gpt2band(g-point)
    band_lims_wvn::Array{T,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
    name::String
    tau::Array{T,3}
    #
    toa_source
    lev_source_inc
    lev_source_dec
  end
  # type, extends(ty_optical_props), public :: ty_source_func_lw
  #   real(wp), allocatable, dimension(:,:,:) :: lay_source,       # Planck source at layer average temperature [W/m2] (ncol, nlay, ngpt)
  #                                              lev_source_inc,   # Planck source at layer edge in increasing ilay direction [W/m2] (ncol, nlay+1, ngpt)
  #                                              lev_source_dec     # Planck source at layer edge in decreasing ilay direction [W/m2] (ncol, nlay+1, ngpt)
  #                                                                 # in increasing/decreasing ilay direction
  #                                                                 # Includes spectral weighting that accounts for state-dependent
  #                                                                 # frequency to g-space mapping
  #   real(wp), allocatable, dimension(:,:  ) :: sfc_source
  # contains
  #   generic,   public :: alloc => alloc_lw, copy_and_alloc_lw
  #   procedure, private:: alloc_lw
  #   procedure, private:: copy_and_alloc_lw
  #   procedure, public :: is_allocated => is_allocated_lw
  #   procedure, public :: finalize => finalize_lw
  #   procedure, public :: get_subset => get_subset_range_lw
  #   procedure, public :: get_ncol => get_ncol_lw
  #   procedure, public :: get_nlay => get_nlay_lw
  #   # validate?
  # end type ty_source_func_lw
  # -------------------------------------------------------------------------------------------------
  #
  # Type for shortave sources: top-of-domain spectrally-resolved flux
  #
  # type, extends(ty_optical_props), public :: ty_source_func_sw
  #   real(wp), allocatable, dimension(:,:  ) :: toa_source
  # contains
  #   generic,   public :: alloc => alloc_sw, copy_and_alloc_sw
  #   procedure, private:: alloc_sw
  #   procedure, private:: copy_and_alloc_sw
  #   procedure, public :: is_allocated => is_allocated_sw
  #   procedure, public :: finalize => finalize_sw
  #   procedure, public :: get_subset => get_subset_range_sw
  #   procedure, public :: get_ncol => get_ncol_sw
  #   # validate?
  # end type ty_source_func_sw
  # -------------------------------------------------------------------------------------------------
# contains
  # ------------------------------------------------------------------------------------------
  #
  #  Routines for initialization, validity checking, finalization
  #
  # ------------------------------------------------------------------------------------------
  #
  # Longwave
  #
  # ------------------------------------------------------------------------------------------
  is_allocated_lw(this::ty_source_func_lw) = is_initialized(this) && allocated(this.sfc_source)
    # class(ty_source_func_lw), intent(in) :: this
    # logical                              :: is_allocated_lw

  # --------------------------------------------------------------
  function alloc_lw!(this::ty_source_func_lw{DT}, ncol, nlay) where DT
    # class(ty_source_func_lw),    intent(inout) :: this
    # integer,                     intent(in   ) :: ncol, nlay
    # character(len = 128)                       :: err_message

    # integer :: ngpt

    !is_initialized(this) && error("source_func_lw%alloc: not initialized so can't allocate")
    any([ncol, nlay] <= 0) && error("source_func_lw%alloc: must provide positive extents for ncol, nlay")

    allocated(this.sfc_source) && deallocate!(this.sfc_source)
    allocated(this.lay_source) && deallocate!(this.lay_source)
    allocated(this.lev_source_inc) && deallocate!(this.lev_source_inc)
    allocated(this.lev_source_dec) && deallocate!(this.lev_source_dec)
    ngpt = get_ngpt(this)
    this.sfc_source = Array(undef, ncol,ngpt)
    this.lay_source = Array(undef, ncol,nlay,ngpt)
    this.lev_source_inc = Array(undef, ncol,nlay,ngpt)
    this.lev_source_dec = Array(undef, ncol,nlay,ngpt)
  end

  # --------------------------------------------------------------
  function copy_and_alloc_lw(this::ty_source_func_lw{DT}, ncol, nlay, spectral_desc) where DT
    # class(ty_source_func_lw),    intent(inout) :: this
    # integer,                     intent(in   ) :: ncol, nlay
    # class(ty_optical_props ),    intent(in   ) :: spectral_desc
    # character(len = 128)                       :: err_message

    !is_initialized(spectral_desc) && error("source_func_lw%alloc: spectral_desc not initialized")
    finalize!(this)
    init!(this, spectral_desc)
    alloc!(this, ncol,nlay)
  end
  # ------------------------------------------------------------------------------------------
  #
  # Shortwave
  #
  # ------------------------------------------------------------------------------------------
  function is_allocated_sw(this::ty_source_func_sw)
    # class(ty_source_func_sw), intent(in) :: this
    # logical                              :: is_allocated_sw

    return is_initialized(this) && allocated(this.toa_source)
  end
  # --------------------------------------------------------------
  function alloc_sw!(this::ty_source_func_sw{DT}, ncol) where DT
    # class(ty_source_func_sw),    intent(inout) :: this
    # integer,                     intent(in   ) :: ncol
    # character(len = 128)                       :: err_message

    !is_initialized(this) && error("source_func_sw%alloc: not initialized so can't allocate")
    ncol <= 0 && error("source_func_sw%alloc: must provide positive extents for ncol")
    allocated(this.toa_source) && deallocate!(this.toa_source)

    this.toa_source = Array(undef, ncol, get_ngpt(this))
  end
  # --------------------------------------------------------------
  function copy_and_alloc_sw!(this::ty_source_func_sw, ncol, spectral_desc)
    # class(ty_source_func_sw),    intent(inout) :: this
    # integer,                     intent(in   ) :: ncol
    # class(ty_optical_props ),    intent(in   ) :: spectral_desc
    # character(len = 128)                       :: err_message

    !is_initialized(spectral_desc) &&  error("source_func_sw%alloc: spectral_desc not initialized")
    init!(this, spectral_desc)
    alloc!(this, ncol)
  end
  # ------------------------------------------------------------------------------------------
  #
  # Finalization (memory deallocation)
  #
  # ------------------------------------------------------------------------------------------
  function finalize_lw!(this::ty_source_func_lw)
    # class(ty_source_func_lw),    intent(inout) :: this

    allocated(this.lay_source    ) && deallocate!(this.lay_source)
    allocated(this.lev_source_inc) && deallocate!(this.lev_source_inc)
    allocated(this.lev_source_dec) && deallocate!(this.lev_source_dec)
    allocated(this.sfc_source    ) && deallocate!(this.sfc_source)
    finalize!(this)
  end
  # --------------------------------------------------------------
  function finalize_sw!(this::ty_source_func_lw)
    # class(ty_source_func_sw),    intent(inout) :: this

    allocated(this.toa_source    ) && deallocate!(this.toa_source)
    finalize!(this)
  end
  # ------------------------------------------------------------------------------------------
  #
  #  Routines for finding the problem size
  #
  # ------------------------------------------------------------------------------------------
  function get_ncol_lw(this::ty_source_func_lw)
    # class(ty_source_func_lw), intent(in) :: this
    # integer :: get_ncol_lw

    if is_allocated(this)
      get_ncol_lw = size(this.lay_source,1)
    else
      get_ncol_lw = 0
    end
  end
  # --------------------------------------------------------------
  function get_nlay_lw(this)
    # class(ty_source_func_lw), intent(in) :: this
    # integer :: get_nlay_lw

    if is_allocated(this)
      get_nlay_lw = size(this.lay_source,2)
    else
      get_nlay_lw = 0
    end
  end
  # --------------------------------------------------------------
  function get_ncol_sw(this)
    # class(ty_source_func_sw), intent(in) :: this
    # integer :: get_ncol_sw

    if is_allocated(this)
      get_ncol_sw = size(this.toa_source,1)
    else
      get_ncol_sw = 0
    end
  end
  # ------------------------------------------------------------------------------------------
  #
  #  Routines for subsetting
  #
  # ------------------------------------------------------------------------------------------
  function get_subset_range_lw(full::ty_source_func_lw, start, n, subset::ty_source_func_lw)
    # class(ty_source_func_lw), intent(inout) :: full
    # integer,                  intent(in   ) :: start, n
    # class(ty_source_func_lw), intent(inout) :: subset
    # character(128)                          :: err_message

    !is_allocated(full) && error("source_func_lw%subset: Asking for a subset of unallocated data")
    if start < 1 || start + n-1 > get_ncol(full)
       error("optical_props%subset: Asking for columns outside range")
     end

    #
    # Could check to see if subset is correctly sized, has consistent spectral discretization
    #
    is_allocated(subset) && finalize!(subset)
    alloc!(subset, n, get_nlay(full), full)

    subset.sfc_source[1:n,  :] = full.sfc_source[start:start+n-1,  :]
    subset.lay_source[1:n,:,:] = full.lay_source[start:start+n-1,:,:]
    subset.lev_source_inc[1:n,:,:] = full.lev_source_inc[start:start+n-1,:,:]
    subset.lev_source_dec[1:n,:,:] = full.lev_source_dec[start:start+n-1,:,:]
  end
  # ------------------------------------------------------------------------------------------
  function get_subset_range_sw(full, start, n, subset)
    # class(ty_source_func_sw), intent(inout) :: full
    # integer,                  intent(in   ) :: start, n
    # class(ty_source_func_sw), intent(inout) :: subset
    # character(128)                          :: err_message

    !is_allocated(full) && error("source_func_sw%subset: Asking for a subset of unallocated data")
    if (start < 1 || start + n-1 > get_ncol(full))
       error("optical_props%subset: Asking for columns outside range")
    end

    #
    # Could check to see if subset is correctly sized, has consistent spectral discretization
    #
    is_allocated(subset) && finalize!(subset)
    # Seems like I should be able to call "alloc" generically but the compilers are complaining
    copy_and_alloc_sw!(subset, n, full)

    subset.toa_source[1:n,  :] = full.toa_source[start:start+n-1,  :]
  end
end # module
