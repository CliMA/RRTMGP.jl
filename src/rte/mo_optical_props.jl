"""
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
# Encapsulate optical properties defined on a spectral grid of N bands.
#   The bands are described by their limiting wavenumbers. They need not be contiguous or complete.
#   A band may contain more than one spectral sub-point (g-point) in which case a mapping must be supplied.
#   A name may be provided and will be prepended to error messages.
#   The base class (ty_optical_props) encapsulates only this spectral discretization and must be initialized
#      with the spectral information before use.
#
#   Optical properties may be represented as arrays with dimensions ncol, nlay, ngpt
#   (abstract class ty_optical_props_arry).
#   The type holds arrays depending on how much information is needed
#   There are three possibilites
#      ty_optical_props_1scl holds absorption optical depth tau, used in calculations accounting for extinction and emission
#      ty_optical_props_2str holds extincion optical depth tau, single-scattering albedo (ssa), and asymmetry parameter g.
#      ty_optical_props_nstr holds extincion optical depth tau, ssa, and phase function moments p with leading dimension nmom.
#   These classes must be allocated before use. Initialization and allocation can be combined.
#   The classes have a validate() function that checks all arrays for valid values (e.g. tau > 0.)
#
# Optical properties can be delta-scaled (though this is currently implemented only for two-stream arrays)
#
# Optical properties can increment or "add themselves to" a set of properties represented with arrays
#   as long as both sets have the same underlying band structure. Properties defined by band
#   may be added to properties defined by g-point; the same value is assumed for all g-points with each band.
#
# Subsets of optical properties held as arrays may be extracted along the column dimension.
#
# -------------------------------------------------------------------------------------------------
"""
module mo_optical_props

using mo_util_array
using mo_optical_props_kernels

  # -------------------------------------------------------------------------------------------------
  #
  # Base class for optical properties
  #   Describes the spectral discretization including the wavenumber limits
  #   of each band (spectral region) and the mapping between g-points and bands
  #
  # -------------------------------------------------------------------------------------------------

abstract type ty_optical_props{T,I} end

struct ty_optical_props_1scl{T,I} <: ty_optical_props{T,I}
  band2gpt::Array{T,2}        # (begin g-point, end g-point) = band2gpt(2,band)
  gpt2band::Array{I,1}        # band = gpt2band(g-point)
  band_lims_wvn::Array{T,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  name::String
  tau::Array{T,3}
end

struct ty_optical_props_2str{T,I} <: ty_optical_props{T,I}
  band2gpt::Array{T,2}        # (begin g-point, end g-point) = band2gpt(2,band)
  gpt2band::Array{I,1}        # band = gpt2band(g-point)
  band_lims_wvn::Array{T,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  name::String
  tau::Array{T,3}
  ssa::Array{T,3}
  g::Array{T,3}
end

struct ty_optical_props_nstr{T,I} <: ty_optical_props{T,I}
  band2gpt::Array{T,2}        # (begin g-point, end g-point) = band2gpt(2,band)
  gpt2band::Array{I,1}        # band = gpt2band(g-point)
  band_lims_wvn::Array{T,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  name::String
  tau::Array{T,3}
  ssa::Array{T,3}
  p::Array{T,4}
end

  # -------------------------------------------------------------------------------------------------
  #
  #  Routines for the base class: initialization, validity checking, finalization
  #
  # -------------------------------------------------------------------------------------------------
  #
  # Base class: Initialization
  #   Values are assumed to be defined in bands a mapping between bands and g-points is provided
  #
  # -------------------------------------------------------------------------------------------------

"""
    init_base!(...)


    class(ty_optical_props),    intent(inout) :: this
    real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn
    integer,  dimension(:,:), &
                      optional, intent(in   ) :: band_lims_gpt
    character(len=*), optional, intent(in   ) :: name
    character(len = 128)                      :: err_message

    integer :: iband
    integer, dimension(2, size(band_lims_wvn, 2)) :: band_lims_gpt_lcl
"""
  function init!(this::ty_optical_props, band_lims_wvn, band_lims_gpt=nothing, name=nothing)
    # -------------------------
    #
    # Error checking -- are the arrays the size we expect, contain positive values?
    #
    err_message = ""
    if size(band_lims_wvn,1) ≠ 2
      err_message = "init(optical_props): band_lims_wvn 1st dim should be 2"
    end
    if any(band_lims_wvn < DT(0))
      err_message = "init(optical_props): band_lims_wvn has values <  0., respectively"
    end
    if length(strip(err_message)) > 0
      return
    end
    if band_lims_gpt ≠ nothing

      if size(band_lims_gpt, 1) ≠ 2
        err_message = "init(optical_props): band_lims_gpt 1st dim should be 2"
      end

      if size(band_lims_gpt,2) ≠ size(band_lims_wvn,2)
        err_message = "init(optical_props): band_lims_gpt and band_lims_wvn sized inconsistently"
      end

      if any(band_lims_gpt < 1)
        err_message = "init(optical_props): band_lims_gpt has values < 1"
      end

      if length(strip(err_message)) > 0
        return
      end

      band_lims_gpt_lcl[:,:] = band_lims_gpt[:,:]
    else
      #
      # Assume that values are defined by band, one g-point per band
      #
      for iband in 1:size(band_lims_wvn, 2)
        band_lims_gpt_lcl[1:2,iband] = iband
      end
    end
    #
    # Assignment
    #
    this.band2gpt = Array(undef, 2,size(band_lims_wvn,2))
    this.band_lims_wvn = Array(undef, 2,size(band_lims_wvn,2)))
    this.band2gpt      = band_lims_gpt_lcl
    this.band_lims_wvn = band_lims_wvn
    if name ≠ nothing
      this.name = strip(name)
    end

    #
    # Make a map between g-points and bands
    #   Efficient only when g-point indexes start at 1 and are contiguous.
    #

    this.gpt2band = Array(undef, maxval(band_lims_gpt_lcl)))
    for iband in 1:size(band_lims_gpt_lcl, dim=2)
      this.gpt2band[band_lims_gpt_lcl[1,iband]:band_lims_gpt_lcl[2,iband]] = iband
    end
    return err_message
  end

"""
    init_base_from_copy!(...)

    class(ty_optical_props),    intent(inout) :: this
    class(ty_optical_props),    intent(in   ) :: spectral_desc
    character(len = 128)                        :: err_message
"""
  function init_base_from_copy!(this::ty_optical_props, spectral_desc::ty_optical_props)
    err_message = init!(this,get_band_lims_wavenumber(spectral_desc), get_band_lims_gpoint(spectral_desc))
    return err_message
  end
  #-------------------------------------------------------------------------------------------------
  #
  # Base class: return true if initialized, false otherwise
  #
  # -------------------------------------------------------------------------------------------------
"""
    is_initialized(...)

    class(ty_optical_props), intent(in) :: this
    logical                             :: is_initialized_base
"""
  is_initialized(this::ty_optical_props) = length(this.band2gpt)>0
  #-------------------------------------------------------------------------------------------------
  #
  # Base class: finalize (deallocate memory)
  #
  # -------------------------------------------------------------------------------------------------
"""
    finalize!(...)

    class(ty_optical_props),    intent(inout) :: this
"""
  function finalize!(this::ty_optical_props)
    this.band2gpt .= 0
    this.gpt2band .= 0
    this.band_lims_wvn .= 0
    this.name = ""
  end
  # ------------------------------------------------------------------------------------------
  #
  #  Routines for array classes: initialization, allocation, and finalization
  #    Initialization and allocation can be combined by supplying either
  #
  # ------------------------------------------------------------------------------------------
  #
  # Straight allocation routines
  #
  # --- 1 scalar ------------------------------------------------------------------------
"""
    alloc_only_1scl!(...)

    class(ty_optical_props_1scl) :: this
    integer,          intent(in) :: ncol, nlay
    character(len=128)           :: err_message
"""
  function alloc!(this::ty_optical_props_1scl, ncol, nlay)
    err_message = ""
    if any([ncol, nlay] .<= 0)
      err_message = "alloc_only_1scl!: must provide positive extents for ncol, nlay"
    else
      this.tau = Array(undef, ncol, nlay, get_ngpt(this))
    end
    return err_message
  end

  # --- 2 stream ------------------------------------------------------------------------
"""
    alloc_only_2str!(...)

    class(ty_optical_props_2str)    :: this
    integer,             intent(in) :: ncol, nlay
    character(len=128)              :: err_message
"""

  function alloc!(this::ty_optical_props_2str, ncol, nlay)
    err_message = ""
    if any([ncol, nlay] .<= 0)
      err_message = "alloc_only_2str!: must provide positive extents for ncol, nlay"
    else
      this.tau = Array(undef, ncol,nlay,get_ngpt(this))
    end
    this.ssa = Array(undef, ncol,nlay,get_ngpt(this))
    this.g = Array(undef, ncol,nlay,get_ngpt(this))
    return err_message
  end

  # --- n stream ------------------------------------------------------------------------

"""
    alloc!()

    class(ty_optical_props_nstr)    :: this
    integer,             intent(in) :: nmom # number of moments
    integer,             intent(in) :: ncol, nlay
    character(len=128)              :: err_message
"""
  function alloc!(this::ty_optical_props_nstr, nmom, ncol, nlay)
    err_message = ""
    if any([ncol, nlay] .<= 0)
      err_message = "optical_props%alloc: must provide positive extents for ncol, nlay"
    else
      this.tau = Array(undef, ncol,nlay,get_ngpt(this))
    end
    this.ssa = Array(undef, ncol,nlay,get_ngpt(this))
    this.p = Array(undef, nmom,ncol,nlay,get_ngpt(this))
    return err_message
  end
  # ------------------------------------------------------------------------------------------
  #
  # Combined allocation/initialization routines
  #
  # ------------------------------------------------------------------------------------------
  #
  # Initialization by specifying band limits and possibly g-point/band mapping
  #
  # ---------------------------------------------------------------------------

"""
    init_and_alloc_1scl!(...)

    class(ty_optical_props_1scl)             :: this
    integer,                      intent(in) :: ncol, nlay
    real(wp), dimension(:,:),     intent(in) :: band_lims_wvn
    integer,  dimension(:,:), &
                      optional,   intent(in) :: band_lims_gpt
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""
  function init_and_alloc!(this::ty_optical_props_1scl, ncol, nlay, band_lims_wvn, band_lims_gpt, name)
    err_message = init!(this, band_lims_wvn, band_lims_gpt, name)
    err_message ≠ "" && return err_message
    err_message = alloc!(this, ncol, nlay)
  end
  # ---------------------------------------------------------------------------
"""
    init_and_alloc!(...)

    class(ty_optical_props_2str)             :: this
    integer,                      intent(in) :: ncol, nlay
    real(wp), dimension(:,:),     intent(in) :: band_lims_wvn
    integer,  dimension(:,:), &
                      optional,   intent(in) :: band_lims_gpt
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""
  function init_and_alloc!(this::ty_optical_props_2str, ncol, nlay, band_lims_wvn, band_lims_gpt, name)
    err_message = init!(this, band_lims_wvn, band_lims_gpt, name)
    err_message ≠ "" && return err_message
    err_message = alloc!(this, ncol, nlay)
  end
  # ---------------------------------------------------------------------------

"""
    init_and_alloc!(...)

    class(ty_optical_props_nstr)             :: this
    integer,                      intent(in) :: nmom, ncol, nlay
    real(wp), dimension(:,:),     intent(in) :: band_lims_wvn
    integer,  dimension(:,:), &
                      optional,   intent(in) :: band_lims_gpt
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""

  function init_and_alloc!(this::ty_optical_props_nstr, nmom, ncol, nlay, band_lims_wvn, band_lims_gpt, name)
    err_message = init!(this, band_lims_wvn, band_lims_gpt, name)
    err_message ≠ "" && return err_message
    err_message = alloc!(this, nmom, ncol, nlay)
  end

  #-------------------------------------------------------------------------------------------------
  #
  # Initialization from an existing spectral discretization/ty_optical_props
  #
  #-------------------------------------------------------------------------------------------------
"""
    copy_and_alloc!(...)

    class(ty_optical_props_1scl)             :: this
    integer,                      intent(in) :: ncol, nlay
    class(ty_optical_props     ), intent(in) :: spectral_desc
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""
  function copy_and_alloc!(this::ty_optical_props_1scl, ncol, nlay, spectral_desc, name)
    err_message = ""
    is_initialized(this) && finalize!(this)
    err_message = init!(this, get_band_lims_wavenumber(spectral_desc),
                              get_band_lims_gpoint(spectral_desc), name)
    err_message ≠ "" && return err_message
    err_message = alloc!(this, ncol, nlay)
  end

"""
    copy_and_alloc!(...)

    class(ty_optical_props_2str)             :: this
    integer,                      intent(in) :: ncol, nlay
    class(ty_optical_props     ), intent(in) :: spectral_desc
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""
  function copy_and_alloc!(this::ty_optical_props_2str, ncol, nlay, spectral_desc, name)
    err_message = ""
    is_initialized(this) && finalize!(this)
    err_message = init!(this, get_band_lims_wavenumber(spectral_desc),
                              get_band_lims_gpoint(spectral_desc), name)
    err_message ≠ "" && return err_message
    err_message = alloc!(this, ncol, nlay)
  end

"""
    copy_and_alloc_nstr!(...)

    class(ty_optical_props_nstr)             :: this
    integer,                      intent(in) :: nmom, ncol, nlay
    class(ty_optical_props     ), intent(in) :: spectral_desc
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""
  function copy_and_alloc!(this::ty_optical_props_nstr, nmom, ncol, nlay, spectral_desc, name)
    err_message = ""
    is_initialized(this) && finalize!(this)
    err_message = init!(this, get_band_lims_wavenumber(spectral_desc),
                              get_band_lims_gpoint(spectral_desc), name)
    err_message ≠ "" && return err_message
    err_message = alloc!(this, nmom, ncol, nlay)
  end

  # ------------------------------------------------------------------------------------------
  #
  #  Routines for array classes: delta-scaling, validation (ensuring all values can be used )
  #
  # ------------------------------------------------------------------------------------------
  # --- delta scaling
  # ------------------------------------------------------------------------------------------
"""
    delta_scale_1scl(...)
"""
  function delta_scale_1scl(this, for) result(err_message)
    class(ty_optical_props_1scl), intent(inout) :: this
    real(wp), dimension(:,:,:), optional, &
                                  intent(in   ) :: for
    character(128)                              :: err_message
    #
    # Nothing to do
    #
    err_message = ""
  end
  # ------------------------------------------------------------------------------------------
  function delta_scale_2str(this, for) result(err_message)
    class(ty_optical_props_2str), intent(inout) :: this
    real(wp), dimension(:,:,:), optional, &
                                  intent(in   ) :: for
    # Forward scattering fraction; g**2 if not provided
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt
    # --------------------------------
    ncol = this%get_ncol()
    nlay = this%get_nlay()
    ngpt = this%get_ngpt()
    err_message = ""

    if(present(for)) then
      if(any([size(for, 1), size(for, 2), size(for, 3)] /= [ncol, nlay, ngpt])) then
        err_message = "delta_scale: dimension of 'for' don't match optical properties arrays"
        return
      end
      if(any(for < 0._wp .or. for > 1._wp)) then
        err_message = "delta_scale: values of 'for' out of bounds [0,1]"
        return
      end
      call delta_scale_2str_kernel(ncol, nlay, ngpt, this%tau, this%ssa, this%g, for)
    else
      call delta_scale_2str_kernel(ncol, nlay, ngpt, this%tau, this%ssa, this%g)
    end

  end
  # ------------------------------------------------------------------------------------------
  function delta_scale_nstr(this, for) result(err_message)
    class(ty_optical_props_nstr), intent(inout) :: this
    real(wp), dimension(:,:,:), optional, &
                                 intent(in   ) :: for
    character(128)                             :: err_message

    err_message = 'delta_scale_nstr: Not yet implemented'
  end
  # ------------------------------------------------------------------------------------------
  #
  # --- Validation
  #
  # ------------------------------------------------------------------------------------------
  function validate_1scalar(this) result(err_message)
    class(ty_optical_props_1scl), intent(in) :: this
    character(len=128)                       :: err_message

    err_message = ''
    if(.not. allocated(this%tau)) then
      err_message = "validate: tau not allocated/initialized"
      return
    end
    if(any_vals_less_than(this%tau, 0._wp)) &
      err_message = "validate: tau values out of range"
    if(len_trim(err_message) > 0 .and. len_trim(this%get_name()) > 0) &
      err_message = trim(this%get_name()) // ': ' // trim(err_message)

  end
  # ------------------------------------------------------------------------------------------
  function validate_2stream(this) result(err_message)
    class(ty_optical_props_2str), intent(in) :: this
    character(len=128)                       :: err_message

    integer :: varSizes(3)

    err_message = ''
    #
    # Array allocation status, sizing
    #
    if(.not. all([allocated(this%tau), allocated(this%ssa), allocated(this%g)])) then
      err_message = "validate: arrays not allocated/initialized"
      return
    end
    varSizes =   [size(this%tau, 1), size(this%tau, 2), size(this%tau, 3)]
    if(.not. all([size(this%ssa, 1), size(this%ssa, 2), size(this%ssa, 3)] == varSizes) .or. &
       .not. all([size(this%g,   1), size(this%g,   2), size(this%g,   3)] == varSizes))     &
    err_message = "validate: arrays not sized consistently"
    #
    # Valid values
    #
    if(any_vals_less_than(this%tau,  0._wp)) &
      err_message = "validate: tau values out of range"
    if(any_vals_outside  (this%ssa,  0._wp, 1._wp)) &
      err_message = "validate: ssa values out of range"
    if(any_vals_outside  (this%g  , -1._wp, 1._wp)) &
      err_message = "validate: g values out of range"

    if(len_trim(err_message) > 0 .and. len_trim(this%get_name()) > 0) &
      err_message = trim(this%get_name()) // ': ' // trim(err_message)

  end

  # ------------------------------------------------------------------------------------------
  function validate_nstream(this) result(err_message)
    class(ty_optical_props_nstr), intent(in) :: this
    character(len=128)                       :: err_message

    integer :: varSizes(3)

    err_message = ''
    #
    # Array allocation status, sizing
    #
    if(.not. all([allocated(this%tau), allocated(this%ssa), allocated(this%p)])) then
      err_message = "validate: arrays not allocated/initialized"
      return
    end
    varSizes =   [size(this%tau, 1), size(this%tau, 2), size(this%tau, 3)]
    if(.not. all([size(this%ssa, 1), size(this%ssa, 2), size(this%ssa, 3)] == varSizes) .or. &
       .not. all([size(this%p,   2), size(this%p,   3), size(this%p,   4)] == varSizes))     &
    err_message = "validate: arrays not sized consistently"
    #
    # Valid values
    #
    if(any_vals_less_than(this%tau,  0._wp)) &
      err_message = "validate: tau values out of range"
    if(any_vals_outside  (this%ssa,  0._wp, 1._wp)) &
      err_message = "validate: ssa values out of range"
    if(any_vals_outside  (this%p(1,:,:,:),  &
                                                                           -1._wp, 1._wp)) &
      err_message = "validate: p(1,:,:,:)  = g values out of range"

    if(len_trim(err_message) > 0 .and. len_trim(this%get_name()) > 0) &
        err_message = trim(this%get_name()) // ': ' // trim(err_message)
  end

  # ------------------------------------------------------------------------------------------
  #
  #  Routines for array classes: subsetting of optical properties arrays along x (col) direction
  #
  # Allocate class, then arrays; copy. Could probably be more efficient if
  #   classes used pointers internally.
  #
  # This set takes start position and number as scalars
  #
  # ------------------------------------------------------------------------------------------

  function subset_1scl_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_1scl), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom

    err_message = ""
    if(.not. full%is_initialized()) then
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return
    end
    ncol = full%get_ncol()
    nlay = full%get_nlay()
    ngpt = full%get_ngpt()
    if(start < 1 .or. start + n-1 > full%get_ncol()) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    if(subset%is_initialized()) call subset%finalize()
    err_message = subset%init(full)
    # Seems like the deallocation statements should be needed under Fortran 2003
    #   but Intel compiler doesn't run without them
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props_1scl)
        err_message = subset%alloc_1scl(n, nlay)
        if(err_message /= "") return
      class is (ty_optical_props_2str)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%alloc_2str(n, nlay)
        if(err_message /= "") return
        subset%ssa(1:n,:,:) = 0._wp
        subset%g  (1:n,:,:) = 0._wp
      class is (ty_optical_props_nstr)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) then
          nmom = subset%get_nmom()
          deallocate(subset%p  )
        else
          nmom = 1
        end
        err_message = subset%alloc_nstr(nmom, n, nlay)
        if(err_message /= "") return
        subset%ssa(1:n,:,:) = 0._wp
        subset%p(:,1:n,:,:) = 0._wp
    end select
    call extract_subset(ncol, nlay, ngpt, full%tau, start, start+n-1, subset%tau)

  end
  # ------------------------------------------------------------------------------------------
  function subset_2str_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_2str), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom

    err_message = ""
    if(.not. full%is_initialized()) then
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return
    end
    ncol = full%get_ncol()
    nlay = full%get_nlay()
    ngpt = full%get_ngpt()
    if(start < 1 .or. start + n-1 > full%get_ncol()) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    if(subset%is_initialized()) call subset%finalize()
    err_message = subset%init(full)
    select type (subset)
      class is (ty_optical_props_1scl)
        err_message = subset%alloc_1scl(n, nlay)
        if(err_message /= "") return
        call extract_subset(ncol, nlay, ngpt, full%tau, full%ssa, start, start+n-1, subset%tau)
      class is (ty_optical_props_2str)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%alloc_2str(n, nlay)
        if(err_message /= "") return
        call extract_subset(ncol, nlay, ngpt, full%tau, start, start+n-1, subset%tau)
        call extract_subset(ncol, nlay, ngpt, full%ssa, start, start+n-1, subset%ssa)
        call extract_subset(ncol, nlay, ngpt, full%g  , start, start+n-1, subset%g  )
      class is (ty_optical_props_nstr)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) then
          nmom = subset%get_nmom()
          deallocate(subset%p  )
        else
          nmom = 1
        end
        err_message = subset%alloc_nstr(nmom, n, nlay)
        if(err_message /= "") return
        call extract_subset(ncol, nlay, ngpt, full%tau, start, start+n-1, subset%tau)
        call extract_subset(ncol, nlay, ngpt, full%ssa, start, start+n-1, subset%ssa)
        subset%p(1,1:n,:,:) = full%g  (start:start+n-1,:,:)
        subset%p(2:,:, :,:) = 0._wp
    end select

  end
  # ------------------------------------------------------------------------------------------
  function subset_nstr_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_nstr), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom

    err_message = ""
    if(.not. full%is_initialized()) then
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return
    end
    ncol = full%get_ncol()
    nlay = full%get_nlay()
    ngpt = full%get_ngpt()
    if(start < 1 .or. start + n-1 > full%get_ncol()) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    if(subset%is_initialized()) call subset%finalize()
    err_message = subset%init(full)
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props_1scl)
        err_message = subset%alloc_1scl(n, nlay)
        if(err_message /= "") return
        call extract_subset(ncol, nlay, ngpt, full%tau, full%ssa, start, start+n-1, subset%tau)
      class is (ty_optical_props_2str)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%alloc_2str(n, nlay)
        if(err_message /= "") return
        call extract_subset(ncol, nlay, ngpt, full%tau, start, start+n-1, subset%tau)
        call extract_subset(ncol, nlay, ngpt, full%ssa, start, start+n-1, subset%ssa)
        subset%g  (1:n,:,:) = full%p(1,start:start+n-1,:,:)
      class is (ty_optical_props_nstr)
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) deallocate(subset%p  )
        err_message = subset%alloc_nstr(nmom, n, nlay)
        if(err_message /= "") return
        call extract_subset(      ncol, nlay, ngpt, full%tau, start, start+n-1, subset%tau)
        call extract_subset(      ncol, nlay, ngpt, full%ssa, start, start+n-1, subset%ssa)
        call extract_subset(nmom, ncol, nlay, ngpt, full%p  , start, start+n-1, subset%p  )
    end select
  end

  # ------------------------------------------------------------------------------------------
  #
  #  Routines for array classes: incrementing
  #   a%increment(b) adds the values of a to b, changing b and leaving a untouched
  #
  # -----------------------------------------------------------------------------------------
  function increment(op_in, op_io) result(err_message)
    class(ty_optical_props_arry), intent(in   ) :: op_in
    class(ty_optical_props_arry), intent(inout) :: op_io
    character(128)                              :: err_message
    # -----
    integer :: ncol, nlay, ngpt, nmom
    # -----
    err_message = ""
    if(.not. op_in%bands_are_equal(op_io)) then
      err_message = "ty_optical_props%increment: optical properties objects have different band structures"
      return
    end
    ncol = op_io%get_ncol()
    nlay = op_io%get_nlay()
    ngpt = op_io%get_ngpt()
    if(op_in%gpoints_are_equal(op_io)) then
      #
      # Increment by gpoint
      #   (or by band if both op_in and op_io are defined that way)
      #
      select type (op_io)
      class is (ty_optical_props_1scl)
          select type (op_in)
           class is (ty_optical_props_1scl)
             call increment_1scalar_by_1scalar(ncol, nlay, ngpt, &
                                               op_io%tau,          &
                                               op_in%tau)
           class is (ty_optical_props_2str)
             call increment_1scalar_by_2stream(ncol, nlay, ngpt, &
                                               op_io%tau,          &
                                               op_in%tau, op_in%ssa)

           class is (ty_optical_props_nstr)
             call increment_1scalar_by_nstream(ncol, nlay, ngpt, &
                                               op_io%tau,          &
                                               op_in%tau, op_in%ssa)
          end select
      class is (ty_optical_props_2str)
        select type (op_in)
          class is (ty_optical_props_1scl)
            call increment_2stream_by_1scalar(ncol, nlay, ngpt,   &
                                              op_io%tau, op_io%ssa,&
                                              op_in%tau)
          class is (ty_optical_props_2str)
            call increment_2stream_by_2stream(ncol, nlay, ngpt,        &
                                              op_io%tau, op_io%ssa, op_io%g, &
                                              op_in%tau, op_in%ssa, op_in%g)
          class is (ty_optical_props_nstr)
            call increment_2stream_by_nstream(ncol, nlay, ngpt, op_in%get_nmom(), &
                                              op_io%tau, op_io%ssa, op_io%g, &
                                              op_in%tau, op_in%ssa, op_in%p)
        end select

      class is (ty_optical_props_nstr)
        select type (op_in)
          class is (ty_optical_props_1scl)
            call increment_nstream_by_1scalar(ncol, nlay, ngpt, &
                                              op_io%tau, op_io%ssa, &
                                              op_in%tau)
          class is (ty_optical_props_2str)
            call increment_nstream_by_2stream(ncol, nlay, ngpt, op_io%get_nmom(), &
                                              op_io%tau, op_io%ssa, op_io%p, &
                                              op_in%tau, op_in%ssa, op_in%g)
          class is (ty_optical_props_nstr)
            call increment_nstream_by_nstream(ncol, nlay, ngpt, op_io%get_nmom(), op_in%get_nmom(), &
                                              op_io%tau, op_io%ssa, op_io%p, &
                                              op_in%tau, op_in%ssa, op_in%p)
        end select
      end select
    else
      #
      # Values defined by-band will have ngpt() = nband()
      # We can use values by band in op_in to increment op_io
      #   Anything else is an error
      #
      if(op_in%get_ngpt() /= op_io%get_nband()) then
        err_message = "ty_optical_props%increment: optical properties objects have incompatible g-point structures"
        return
      end
      #
      # Increment by band
      #
      select type (op_io)
        class is (ty_optical_props_1scl)
          select type (op_in)
          class is (ty_optical_props_1scl)
              call inc_1scalar_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                                op_io%tau,          &
                                                op_in%tau,          &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
            class is (ty_optical_props_2str)
              call inc_1scalar_by_2stream_bybnd(ncol, nlay, ngpt, &
                                                op_io%tau,          &
                                                op_in%tau, op_in%ssa, &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
            class is (ty_optical_props_nstr)
              call inc_1scalar_by_nstream_bybnd(ncol, nlay, ngpt, &
                                                op_io%tau,          &
                                                op_in%tau, op_in%ssa, &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
          end select

        class is (ty_optical_props_2str)
          select type (op_in)
            class is (ty_optical_props_1scl)
              call inc_2stream_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                                op_io%tau, op_io%ssa, &
                                                op_in%tau,          &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
            class is (ty_optical_props_2str)
              call inc_2stream_by_2stream_bybnd(ncol, nlay, ngpt,        &
                                                op_io%tau, op_io%ssa, op_io%g, &
                                                op_in%tau, op_in%ssa, op_in%g, &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
            class is (ty_optical_props_nstr)
              call inc_2stream_by_nstream_bybnd(ncol, nlay, ngpt, op_in%get_nmom(), &
                                                op_io%tau, op_io%ssa, op_io%g, &
                                                op_in%tau, op_in%ssa, op_in%p, &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
          end select

        class is (ty_optical_props_nstr)
          select type (op_in)
            class is (ty_optical_props_1scl)
              call inc_nstream_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                                op_io%tau, op_io%ssa, &
                                                op_in%tau,          &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
            class is (ty_optical_props_2str)
              call inc_nstream_by_2stream_bybnd(ncol, nlay, ngpt, op_io%get_nmom(), &
                                                op_io%tau, op_io%ssa, op_io%p, &
                                                op_in%tau, op_in%ssa, op_in%g, &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
            class is (ty_optical_props_nstr)
              call inc_nstream_by_nstream_bybnd(ncol, nlay, ngpt, op_io%get_nmom(), op_in%get_nmom(), &
                                                op_io%tau, op_io%ssa, op_io%p, &
                                                op_in%tau, op_in%ssa, op_in%p, &
                                                op_io%get_nband(), op_io%get_band_lims_gpoint())
          end select
      end select
    end
  end
  # -----------------------------------------------------------------------------------------------
  #
  #  Routines for array classes: problem sizes
  #
  # -----------------------------------------------------------------------------------------------
  function get_arry_extent(this, dim)
    class(ty_optical_props_arry), intent(in   ) :: this
    integer,                      intent(in   ) :: dim
    integer                                     :: get_arry_extent

    if(allocated(this%tau)) then
      get_arry_extent = size(this%tau, dim)
    else
      get_arry_extent = 0
    end
  end
  # ------------------------------------------------------------------------------------------
  function get_ncol(this)
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_ncol

    get_ncol = get_arry_extent(this, 1)
  end
  # ------------------------------------------------------------------------------------------
  function get_nlay(this)
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_nlay

    get_nlay = get_arry_extent(this, 2)
  end
  # ------------------------------------------------------------------------------------------
  function get_nmom(this)
    class(ty_optical_props_nstr), intent(in   ) :: this
    integer                                     :: get_nmom

    if(allocated(this%p)) then
      get_nmom = size(this%p, 1)
    else
      get_nmom = 0
    end
  end
  # -----------------------------------------------------------------------------------------------
  #
  #  Routines for base class: spectral discretization
  #
  # -----------------------------------------------------------------------------------------------
  #
  # Number of bands
  #
  function get_nband(this)
    class(ty_optical_props), intent(in) :: this
    integer                             :: get_nband

    if(this%is_initialized()) then
      get_nband = size(this%band2gpt,dim=2)
    else
      get_nband = 0
    end
  end
  # -----------------------------------------------------------------------------------------------
  #
  # Number of g-points
  #
  function get_ngpt(this)
    class(ty_optical_props), intent(in) :: this
    integer                             :: get_ngpt

    if(this%is_initialized()) then
      get_ngpt = maxval(this%band2gpt)
    else
      get_ngpt = 0
    end
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # The first and last g-point of all bands at once
  # dimension (2, nbands)
  #
  function get_band_lims_gpoint(this)
    class(ty_optical_props), intent(in) :: this
    integer, dimension(size(this%band2gpt,dim=1), size(this%band2gpt,dim=2)) &
                                        :: get_band_lims_gpoint

    get_band_lims_gpoint = this%band2gpt
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # First and last g-point of a specific band
  #
  function convert_band2gpt(this, band)
    class(ty_optical_props), intent(in) :: this
    integer,                 intent(in) :: band
    integer, dimension(2)               :: convert_band2gpt

    if(this%is_initialized()) then
      convert_band2gpt(:) = this%band2gpt(:,band)
    else
      convert_band2gpt(:) = 0
    end
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Lower and upper wavenumber of all bands
  # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  #
  function get_band_lims_wavenumber(this)
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2)) &
                                        :: get_band_lims_wavenumber

    if(this%is_initialized()) then
      get_band_lims_wavenumber(:,:) = this%band_lims_wvn(:,:)
    else
      get_band_lims_wavenumber(:,:) = 0._wp
    end
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Lower and upper wavelength of all bands
  #
  function get_band_lims_wavelength(this)
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2)) &
                                        :: get_band_lims_wavelength

    if(this%is_initialized()) then
      get_band_lims_wavelength(:,:) = 1._wp/this%band_lims_wvn(:,:)
    else
      get_band_lims_wavelength(:,:) = 0._wp
    end
  end
  #--------------------------------------------------------------------------------------------------------------------
  # Bands for all the g-points at once
  # dimension (ngpt)
  #
  function get_gpoint_bands(this)
    class(ty_optical_props), intent(in) :: this
    integer, dimension(size(this%gpt2band,dim=1)) &
                                        :: get_gpoint_bands

    if(this%is_initialized()) then
      get_gpoint_bands(:) = this%gpt2band(:)
    else
      get_gpoint_bands(:) = 0
    end
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Band associated with a specific g-point
  #
  function convert_gpt2band(this, gpt)
    class(ty_optical_props), intent(in) :: this
    integer,                            intent(in) :: gpt
    integer                             :: convert_gpt2band

    if(this%is_initialized()) then
      convert_gpt2band = this%gpt2band(gpt)
    else
      convert_gpt2band = 0
    end
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Expand an array of dimension arr_in(nband) to dimension arr_out(ngpt)
  #
  function expand(this, arr_in) result(arr_out)
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(:),  intent(in) :: arr_in # (nband)
    real(wp), dimension(size(this%gpt2band)) :: arr_out

    integer :: iband

    do iband=1,this%get_nband()
      arr_out(this%band2gpt(1,iband):this%band2gpt(2,iband)) = arr_in(iband)
    end
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Are the bands of two objects the same? (same number, same wavelength limits)
  #
  function bands_are_equal(this, that)
    class(ty_optical_props), intent(in) :: this, that
    logical                             :: bands_are_equal

    bands_are_equal = this%get_nband() == that%get_nband() .and. &
                      this%get_nband() > 0
    if(.not. bands_are_equal) return
    bands_are_equal = &
      all(abs(this%get_band_lims_wavenumber() - that%get_band_lims_wavenumber()) < &
          5._wp * spacing(this%get_band_lims_wavenumber()))
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Is the g-point structure of two objects the same?
  #   (same bands, same number of g-points, same mapping between bands and g-points)
  #
  function gpoints_are_equal(this, that)
    class(ty_optical_props), intent(in) :: this, that
    logical                             :: gpoints_are_equal

    gpoints_are_equal = this%bands_are_equal(that) .and. &
                        this%get_ngpt() == that%get_ngpt()
    if(.not. gpoints_are_equal) return
    gpoints_are_equal = &
      all(this%get_gpoint_bands() == that%get_gpoint_bands())
  end
  # -----------------------------------------------------------------------------------------------
  #
  # --- Setting/getting the name
  #
  # -----------------------------------------------------------------------------------------------
  function set_name!(this, name)
    class(ty_optical_props),  intent(inout) :: this
    character(len=*),         intent(in   ) :: name

    this%name = trim(name)
  end
  # --------------------------------------------------------
  function get_name(this)
    class(ty_optical_props),  intent(in   ) :: this
    character(len=name_len)                 :: get_name

      get_name = trim(this%name)
  end
  # ------------------------------------------------------------------------------------------

end # module
