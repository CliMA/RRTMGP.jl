"""
This code is part of Radiative Transfer for Energetics (RTE)

Contacts: Robert Pincus and Eli Mlawer
email:  rrtmgp@aer.com

Copyright 2015-2018,  Atmospheric and Environmental Research and
Regents of the University of Colorado.  All right reserved.

Use and duplication is permitted under the terms of the
  BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
-------------------------------------------------------------------------------------------------

Encapsulate optical properties defined on a spectral grid of N bands.
 The bands are described by their limiting wavenumbers. They need not be contiguous or complete.
 A band may contain more than one spectral sub-point (g-point) in which case a mapping must be supplied.
 A name may be provided and will be prepended to error messages.
 The base class (ty_optical_props) encapsulates only this spectral discretization and must be initialized
    with the spectral information before use.

 Optical properties may be represented as arrays with dimensions ncol, nlay, ngpt
 (abstract class ty_optical_props_arry).
 The type holds arrays depending on how much information is needed
 There are three possibilites
    ty_optical_props_1scl holds absorption optical depth tau, used in calculations accounting for extinction and emission
    ty_optical_props_2str holds extincion optical depth tau, single-scattering albedo (ssa), and asymmetry parameter g.
    ty_optical_props_nstr holds extincion optical depth tau, ssa, and phase function moments p with leading dimension nmom.
 These classes must be allocated before use. Initialization and allocation can be combined.
 The classes have a validate() function that checks all arrays for valid values (e.g. tau > 0.)

Optical properties can be delta-scaled (though this is currently implemented only for two-stream arrays)

Optical properties can increment or "add themselves to" a set of properties represented with arrays
 as long as both sets have the same underlying band structure. Properties defined by band
 may be added to properties defined by g-point; the same value is assumed for all g-points with each band.

Subsets of optical properties held as arrays may be extracted along the column dimension.

-------------------------------------------------------------------------------------------------
"""
module mo_optical_props

using ..fortran_intrinsics
using ..mo_util_array

export init!,
       init_base_from_copy!,
       finalize!,
       alloc!,
       init_and_alloc!,
       copy_and_alloc!,
       delta_scale!,
       validate!,
       subset_range!,
       increment!,
       convert_band2gpt,
       get_band_lims_wavenumber,
       get_band_lims_wavelength,
       get_gpoint_bands,
       expand,
       bands_are_equal,
       convert_gpt2band,
       gpoints_are_equal

export ty_optical_props,
       ty_optical_props_1scl,
       ty_optical_props_2str,
       ty_optical_props_nstr

any_vals_less_than(a, t) = a .< t
any_vals_outside(a,t1,t2) = !(t1<a<t2)

  # -------------------------------------------------------------------------------------------------
  #
  # Base class for optical properties
  #   Describes the spectral discretization including the wavenumber limits
  #   of each band (spectral region) and the mapping between g-points and bands
  #
  # -------------------------------------------------------------------------------------------------

export ty_optical_props
abstract type ty_optical_props{T,I} end
export ty_optical_props_arry
abstract type ty_optical_props_arry{T,I} <: ty_optical_props{T,I} end

struct ty_optical_props_1scl{T,I} <: ty_optical_props{T,I}
  band2gpt#::Array{T,2}        # (begin g-point, end g-point) = band2gpt(2,band)
  gpt2band#::Array{I,1}        # band = gpt2band(g-point)
  band_lims_wvn#::Array{T,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  name#::String
  tau#::Array{T,3}
end

struct ty_optical_props_2str{T,I} <: ty_optical_props{T,I}
  band2gpt#::Array{T,2}        # (begin g-point, end g-point) = band2gpt(2,band)
  gpt2band#::Array{I,1}        # band = gpt2band(g-point)
  band_lims_wvn#::Array{T,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  name#::String
  tau#::Array{T,3}
  ssa#::Array{T,3}
  g#::Array{T,3}
end

struct ty_optical_props_nstr{T,I} <: ty_optical_props{T,I}
  band2gpt#::Array{T,2}        # (begin g-point, end g-point) = band2gpt(2,band)
  gpt2band#::Array{I,1}        # band = gpt2band(g-point)
  band_lims_wvn#::Array{T,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  name#::String
  tau#::Array{T,3}
  ssa#::Array{T,3}
  p#::Array{T,4}
end


# ty_optical_props_1scl() = ty_optical_props_1scl(nothing,nothing,nothing,nothing,nothing)
# ty_optical_props_2str() = ty_optical_props_2str(nothing,nothing,nothing,nothing,nothing,nothing,nothing)
# ty_optical_props_nstr() = ty_optical_props_nstr(nothing,nothing,nothing,nothing,nothing,nothing,nothing)

# band_lims_wvn, band_lims_gpt
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
    integer,  dimension(:,:),
                      optional, intent(in   ) :: band_lims_gpt
    character(len=*), optional, intent(in   ) :: name
    character(len = 128)                      :: err_message

    integer :: iband
    integer, dimension(2, size(band_lims_wvn, 2)) :: band_lims_gpt_lcl
"""
  function ty_optical_props_1scl(name, band_lims_wvn, band_lims_gpt=nothing)
    DT = eltype(band_lims_wvn)

    @assert size(band_lims_wvn,1) == 2
    @show typeof(band_lims_wvn)
    @assert !any(band_lims_wvn.<DT(0))
    @show band_lims_gpt
    band_lims_gpt_lcl = Array{Int}(undef, 2, size(band_lims_wvn, 2))
    if band_lims_gpt ≠ nothing
      @assert size(band_lims_gpt,1) == 2
      @assert size(band_lims_gpt,2) == size(band_lims_wvn,2)
      @assert !any(band_lims_gpt .< 1)
      band_lims_gpt_lcl[:,:] .= band_lims_gpt[:,:]
    else
      for iband in 1:size(band_lims_wvn, 2)
        band_lims_gpt_lcl[1:2,iband] = iband
      end
    end

    band2gpt = band_lims_gpt_lcl

    # Make a map between g-points and bands
    #   Efficient only when g-point indexes start at 1 and are contiguous.

    gpt2band = Array{Int}(undef, max(band_lims_gpt_lcl...))
    for iband in 1:size(band_lims_gpt_lcl, 2)
      gpt2band[band_lims_gpt_lcl[1,iband]:band_lims_gpt_lcl[2,iband]] .= iband
    end
    tau = nothing
    return ty_optical_props_1scl{DT, Int}(band2gpt,gpt2band,band_lims_wvn,name,tau)
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
  is_initialized(this::ty_optical_props) = allocated(this.band2gpt)
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
    allocated(this.band2gpt) && deallocate!(this.band2gpt)
    allocated(this.gpt2band) && deallocate!(this.gpt2band)
    allocated(this.band_lims_wvn) && deallocate!(this.band_lims_wvn)
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
      allocated(this.tau) && deallocate!(this.tau)
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
      allocated(this.tau) && deallocate!(this.tau)
      this.tau = Array(undef, ncol,nlay,get_ngpt(this))
    end
    allocated(this.ssa) && deallocate!(this.ssa)
    this.ssa = Array(undef, ncol,nlay,get_ngpt(this))
    allocated(this.g) && deallocate!(this.g)
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
      allocated(this.tau) && deallocate!(this.tau)
      this.tau = Array(undef, ncol,nlay,get_ngpt(this))
    end
    allocated(this.ssa) && deallocate!(this.ssa)
    this.ssa = Array(undef, ncol,nlay,get_ngpt(this))
    allocated(this.p) && deallocate!(this.p)
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
    integer,  dimension(:,:),
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
    integer,  dimension(:,:),
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
    integer,  dimension(:,:),
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
    delta_scale!(...)

    class(ty_optical_props_1scl), intent(inout) :: this
    real(wp), dimension(:,:,:), optional,
                                  intent(in   ) :: for_
    character(128)                              :: err_message
"""
  delta_scale!(this::ty_optical_props_1scl, for_) = ""


"""
    delta_scale!()

    class(ty_optical_props_2str), intent(inout) :: this
    real(wp), dimension(:,:,:), optional,
                                  intent(in   ) :: for_
    # Forward scattering fraction; g**2 if not provided
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt
    # --------------------------------
"""
  function delta_scale!(this::ty_optical_props_2str{DT}, for_ = nothing) where DT
    ncol = get_ncol(this)
    nlay = get_nlay(this)
    ngpt = get_ngpt(this)
    err_message = ""

    if for_ ≠ nothing
      if any([size(for_, 1), size(for_, 2), size(for_, 3)] ≠ [ncol, nlay, ngpt])
        err_message = "delta_scale: dimension of 'for_' don't match optical properties arrays"
        return err_message
      end
      if any(for_ < DT(0) || for_ > DT(1))
        err_message = "delta_scale: values of 'for_' out of bounds [0,1]"
        return err_message
      end
      delta_scale_2str_kernel!(ncol, nlay, ngpt, this.tau, this.ssa, this.g, for_)
    else
      delta_scale_2str_kernel!(ncol, nlay, ngpt, this.tau, this.ssa, this.g)
    end
    return err_message
  end
  # ------------------------------------------------------------------------------------------

"""
    delta_scale!(...)

    class(ty_optical_props_nstr), intent(inout) :: this
    real(wp), dimension(:,:,:), optional,
                                 intent(in   ) :: for_
    character(128)                             :: err_message
"""
  delta_scale!(this::ty_optical_props_nstr, for_) = "delta_scale_nstr: Not yet implemented"

  # ------------------------------------------------------------------------------------------
  #
  # --- Validation
  #
  # ------------------------------------------------------------------------------------------
"""
    validate!(...)

    class(ty_optical_props_1scl), intent(in) :: this
    character(len=128)                       :: err_message
"""

  function validate!(this::ty_optical_props_1scl{DT}) where DT
    err_message = ""
    if !allocated(this)
      err_message = "validate: tau not allocated/initialized"
      return err_message
    end
    any_vals_less_than(this.tau, DT(0)) && (err_message = "validate: tau values out of range")
    if len_trim(err_message) > 0 && len_trim(this.name) > 0
      err_message = strip(this.name) * ": " * strip(err_message)
    end
    return err_message
  end

"""

    class(ty_optical_props_2str), intent(in) :: this
    character(len=128)                       :: err_message

    integer :: varSizes(3)
"""
  function validate!(this::ty_optical_props_2str{DT}) where DT

    err_message = ""
    #
    # Array allocation status, sizing
    #
    if !all([allocated(this.tau), allocated(this.ssa), allocated(this.g)])
      err_message = "validate: arrays not allocated/initialized"
      return err_message
    end
    varSizes =   [size(this.tau, 1), size(this.tau, 2), size(this.tau, 3)]
    if !all([size(this.ssa, 1), size(this.ssa, 2), size(this.ssa, 3)] == varSizes) ||
       !all([size(this.g,   1), size(this.g,   2), size(this.g,   3)] == varSizes)
    err_message = "validate: arrays not sized consistently"
    end
    #
    # Valid values
    #
    any_vals_less_than(this.tau,  DT(0)) && (err_message = "validate: tau values out of range")
    any_vals_outside(this.ssa,  DT(0), DT(1)) && (err_message = "validate: ssa values out of range")
    any_vals_outside(this.g  , DT(-1), DT(1)) && (err_message = "validate: g values out of range")

    if(len_trim(err_message) > 0 && len_trim(this.name) > 0)
      err_message = trim(this.name) * ": " * trim(err_message)
    end

  end

  # ------------------------------------------------------------------------------------------
"""

    class(ty_optical_props_nstr), intent(in) :: this
    character(len=128)                       :: err_message

    integer :: varSizes(3)
"""
  function validate!(this::ty_optical_props_nstr{DT}) where DT

    err_message = ""
    #
    # Array allocation status, sizing
    #
    if !all([allocated(this.tau), allocated(this.ssa), allocated(this.p)])
      err_message = "validate: arrays not allocated/initialized"
      return err_message
    end
    varSizes =   [size(this.tau, 1), size(this.tau, 2), size(this.tau, 3)]
    if      !all([size(this.ssa, 1), size(this.ssa, 2), size(this.ssa, 3)] == varSizes) ||
            !all([size(this.p,   2), size(this.p,   3), size(this.p,   4)] == varSizes)
      err_message = "validate: arrays not sized consistently"
    end
    #
    # Valid values
    #
    any_vals_less_than(this.tau,  DT(0)) && (err_message = "validate: tau values out of range")
    any_vals_outside(this.ssa,  DT(0), DT(1)) && (err_message = "validate: ssa values out of range")
    any_vals_outside(this.p[1,:,:,:], DT(-1), DT(1)) && (err_message = "validate: p(1,:,:,:)  = g values out of range")

    if len_trim(err_message) > 0 && len_trim(this.name) > 0
      err_message = trim(this.name) * ": " * trim(err_message)
    end
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

"""
    class(ty_optical_props_1scl), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom
"""
  function subset_range!(full::ty_optical_props_1scl{DT}, start, n, subset) where DT
    err_message = ""
    if !is_initialized(full)
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return err_message
    end
    ncol = get_ncol(full)
    nlay = get_nlay(full)
    ngpt = get_ngpt(full)
    if start < 1 || start + n-1 > get_ncol(full)
       err_message = "optical_props%subset: Asking for columns outside range"
    end
    err_message ≠ "" && return

    is_initialized(subset) && finalize!(subset)
    err_message = init!(subset, full)
    # Seems like the deallocation statements should be needed under Fortran 2003
    #   but Intel compiler doesn't run without them

    allocated(subset.tau) && deallocate!(subset.tau)
    if subset isa ty_optical_props_1scl
        err_message = alloc!(subset, n, nlay)
        err_message ≠ "" && return
    elseif subset isa ty_optical_props_2str
        allocated(subset.ssa) && deallocate!(subset.ssa)
        allocated(subset.g  ) && deallocate!(subset.g  )
        err_message = alloc!(subset, n, nlay)
        err_message ≠ "" && return
        subset.ssa[1:n,:,:] .= DT(0)
        subset.g[1:n,:,:] .= DT(0)
    elseif subset isa ty_optical_props_nstr
        allocated(subset.ssa) && deallocate!(subset.ssa)
        if allocated(subset.p)
          nmom = get_nmom(subset)
          allocated(subset.p  ) && deallocate!(subset.p  )
        else
          nmom = 1
        end
        err_message = alloc!(subset, nmom, n, nlay)
        err_message ≠ "" && return
        subset.ssa[1:n,:,:] .= DT(0)
        subset.p[:,1:n,:,:] .= DT(0)
    else
      error("Uncaught case in subset_range!(full::ty_optical_props_1scl{DT}")
    end
    extract_subset!(ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)

  end
  # ------------------------------------------------------------------------------------------

"""
    subset_range!(...)

    class(ty_optical_props_2str), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom
"""
  function subset_range!(full::ty_optical_props_2str{DT}, start, n, subset) where DT

    err_message = ""
    if !is_initialized(full)
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return err_message
    end
    ncol = get_ncol(full)
    nlay = get_nlay(full)
    ngpt = get_ngpt(full)
    if start < 1 || start + n-1 > get_ncol(full)
       err_message = "optical_props%subset: Asking for columns outside range"
    end
    err_message ≠ "" && return

    is_initialized(subset) && finalize!(subset)
    err_message = init(subset, full)

    if subset isa ty_optical_props_1scl # TODO: check logic
      err_message = alloc!(subset, n, nlay)
      err_message ≠ "" && return err_message
      extract_subset!(ncol, nlay, ngpt, full.tau, full.ssa, start, start+n-1, subset.tau)
    elseif subset isa ty_optical_props_2str
      allocated(subset.ssa) && deallocate!(subset.ssa)
      allocated(subset.g  ) && deallocate!(subset.g  )
      err_message = alloc!(subset, n, nlay)
      err_message ≠ "" && return
      extract_subset!(ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)
      extract_subset!(ncol, nlay, ngpt, full.ssa, start, start+n-1, subset.ssa)
      extract_subset!(ncol, nlay, ngpt, full.g  , start, start+n-1, subset.g  )
    elseif subset isa ty_optical_props_nstr
      allocated(subset.ssa) && deallocate!(subset.ssa)
      if allocated(subset.p)
        nmom = get_nmom(subset)
        allocated(subset.p  ) && deallocate!(subset.p  )
      else
        nmom = 1
      end
      err_message = alloc!(subset, nmom, n, nlay)
      err_message ≠ "" && return
      extract_subset!(ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)
      extract_subset!(ncol, nlay, ngpt, full.ssa, start, start+n-1, subset.ssa)
      subset.p[1,1:n,:,:] .= full.g[start:start+n-1,:,:]
      subset.p[2:end,:, :,:] .= DT(0) # TODO Verify this line
    else
      error("Uncaught case in subset_range!(full::ty_optical_props_2str{DT}")
    end
  end

  # ------------------------------------------------------------------------------------------
"""
    subset_range!(full::ty_optical_props_nstr, start, n, subset)

    class(ty_optical_props_nstr), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom
"""
  function subset_range!(full::ty_optical_props_nstr, start, n, subset)

    err_message = ""
    if !is_initialized(full)
      err_message = "optical_props%subset: Asking for a subset of uninitialized data"
      return err_message
    end
    ncol = get_ncol(full)
    nlay = get_nlay(full)
    ngpt = get_ngpt(full)
    if(start < 1 || start + n-1 > get_ncol(full))
       err_message = "optical_props%subset: Asking for columns outside range"
    end
    err_message ≠ "" && return

    is_initialized(subset) && finalize!(subset)
    err_message = init!(subset, full)

    allocated(subset.tau) && deallocate!(subset.tau)
    if subset isa ty_optical_props_1scl # TODO: check logic
      err_message = alloc!(subset, n, nlay)
      err_message ≠ "" && return
      extract_subset!(ncol, nlay, ngpt, full.tau, full.ssa, start, start+n-1, subset.tau)
    elseif subset isa ty_optical_props_2str
      allocated(subset.ssa) && deallocate!(subset.ssa)
      allocated(subset.g  ) && deallocate!(subset.g  )
      err_message = alloc!(subset, n, nlay)
      err_message ≠ "" && return
      extract_subset!(ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)
      extract_subset!(ncol, nlay, ngpt, full.ssa, start, start+n-1, subset.ssa)
      subset.g[1:n,:,:] .= full.p[1,start:start+n-1,:,:]
    elseif subset isa ty_optical_props_nstr
      allocated(subset.ssa) && deallocate!(subset.ssa)
      allocated(subset.p  ) && deallocate!(subset.p  )
      err_message = alloc!(subset, nmom, n, nlay)
      err_message ≠ "" && return
      extract_subset!(      ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)
      extract_subset!(      ncol, nlay, ngpt, full.ssa, start, start+n-1, subset.ssa)
      extract_subset!(nmom, ncol, nlay, ngpt, full.p  , start, start+n-1, subset.p  )
    else
      error("Uncaught case in subset_range!(full::ty_optical_props_nstr, start, n, subset)")
    end
  end

  # ------------------------------------------------------------------------------------------
  #
  #  Routines for array classes: incrementing
  #   a%increment(b) adds the values of a to b, changing b and leaving a untouched
  #
  # -----------------------------------------------------------------------------------------
"""
    increment(op_in, op_io)

    class(ty_optical_props_arry), intent(in   ) :: op_in
    class(ty_optical_props_arry), intent(inout) :: op_io
    character(128)                              :: err_message
    # -----
    integer :: ncol, nlay, ngpt, nmom
    # -----
"""
  function increment!(op_in, op_io)
    err_message = ""
    if !bands_are_equal(op_in, op_io)
      err_message = "ty_optical_props%increment: optical properties objects have different band structures"
      return
    end
    ncol = get_ncol(op_io)
    nlay = get_nlay(op_io)
    ngpt = get_ngpt(op_io)
    if gpoints_are_equal(op_in, op_io)
      #
      # Increment by gpoint
      #   (or by band if both op_in and op_io are defined that way)
      #
      increment_by_gpoint!(op_io, op_in)

    else

      # Values defined by-band will have ngpt() = nband()
      # We can use values by band in op_in to increment op_io
      #   Anything else is an error

      if get_ngpt(op_in) ≠ get_nband(op_io)
        err_message = "ty_optical_props%increment: optical properties objects have incompatible g-point structures"
        return err_message
      end
      #
      # Increment by band
      #
      increment_bybnd!(op_io, op_in)
    end
  end

  # -----------------------------------------------------------------------------------------------
  #
  #  Routines for array classes: problem sizes
  #
  # -----------------------------------------------------------------------------------------------

"""
    class(ty_optical_props_arry), intent(in   ) :: this
    integer,                      intent(in   ) :: dim
    integer                                     :: get_arry_extent
"""
  get_arry_extent(this::ty_optical_props, dim) = allocated(this.tau) ? size(this.tau, dim) : 0
  # ------------------------------------------------------------------------------------------

"""
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_ncol
"""
  get_ncol(this::ty_optical_props) = get_arry_extent(this, 1)

  # ------------------------------------------------------------------------------------------

"""
    class(ty_optical_props_arry), intent(in   ) :: this
    integer                                     :: get_nlay
"""
  get_nlay(this::ty_optical_props) = get_arry_extent(this, 2)

  # ------------------------------------------------------------------------------------------

"""
    class(ty_optical_props_nstr), intent(in   ) :: this
    integer                                     :: get_nmom
"""
  get_nmom(this::ty_optical_props) = allocated(this.p) ? size(this.p, 1) : 0

  # -----------------------------------------------------------------------------------------------
  #
  #  Routines for base class: spectral discretization
  #
  # -----------------------------------------------------------------------------------------------
  #
  # Number of bands
  #

"""
    class(ty_optical_props), intent(in) :: this
    integer                             :: get_nband
"""
  get_nband(this::ty_optical_props) = is_initialized(this) ? size(this.band2gpt,dim=2) : 0

  # -----------------------------------------------------------------------------------------------
  #
  # Number of g-points
  #
"""
    class(ty_optical_props), intent(in) :: this
    integer                             :: get_ngpt
"""
  get_ngpt(this::ty_optical_props) = is_initialized(this) ? maxval(this.band2gpt) : 0

  #--------------------------------------------------------------------------------------------------------------------
  #
  # The first and last g-point of all bands at once
  # dimension (2, nbands)
  #

"""
    class(ty_optical_props), intent(in) :: this
    integer, dimension(size(this%band2gpt,dim=1), size(this%band2gpt,dim=2))
                                        :: get_band_lims_gpoint
"""
  get_band_lims_gpoint(this::ty_optical_props) = this.band2gpt

  #--------------------------------------------------------------------------------------------------------------------
  #
  # First and last g-point of a specific band
  #
"""
    class(ty_optical_props), intent(in) :: this
    integer,                 intent(in) :: band
    integer, dimension(2)               :: convert_band2gpt
"""
  function convert_band2gpt(this::ty_optical_props{DT}, band) where DT
    is_initialized(this) ? this.band2gpt[:,band] : zeros(DT,length(this.band2gpt[:,band]))
  end

  #--------------------------------------------------------------------------------------------------------------------
  #
  # Lower and upper wavenumber of all bands
  # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  #
"""
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2))
                                        :: get_band_lims_wavenumber
"""
  function get_band_lims_wavenumber(this::ty_optical_props{DT}) where DT
    is_initialized(this) ? this.band_lims_wvn[:,:] : zeros(DT, size(this.band_lims_wvn))
  end

  #--------------------------------------------------------------------------------------------------------------------
  #
  # Lower and upper wavelength of all bands
  #
"""
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2))
                                        :: get_band_lims_wavelength
"""
  function get_band_lims_wavelength(this::ty_optical_props{DT}) where DT
    is_initialized(this) ? 1 ./ this.band_lims_wvn[:,:] : zeros(DT, size(this.band_lims_wvn))
  end

  #--------------------------------------------------------------------------------------------------------------------
  # Bands for all the g-points at once
  # dimension (ngpt)
  #
"""
    class(ty_optical_props), intent(in) :: this
    integer, dimension(size(this%gpt2band,dim=1))
                                        :: get_gpoint_bands
"""
  function get_gpoint_bands(this::ty_optical_props)
    is_initialized(this) ? this.gpt2band[:] : zeros(Int,length(this.gpt2band))
  end

  #--------------------------------------------------------------------------------------------------------------------
  #
  # Band associated with a specific g-point
  #

"""
    class(ty_optical_props), intent(in) :: this
    integer,                            intent(in) :: gpt
    integer                             :: convert_gpt2band
"""
  convert_gpt2band(this::ty_optical_props, gpt) = is_initialized(this) ? this.gpt2band[gpt] : 0

  #--------------------------------------------------------------------------------------------------------------------
  #
  # Expand an array of dimension arr_in(nband) to dimension arr_out(ngpt)
  #
"""
    class(ty_optical_props), intent(in) :: this
    real(wp), dimension(:),  intent(in) :: arr_in # (nband)
    real(wp), dimension(size(this%gpt2band)) :: arr_out

    integer :: iband
"""
  function expand(this::ty_optical_props{DT}, arr_in) where DT
    arr_out = Vector{DT}(undef, size(this.gpt2band))
    for iband in 1:get_nband(this)
      arr_out[this.band2gpt[1,iband]:this.band2gpt[2,iband]] = arr_in[iband]
    end
    return arr_out
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Are the bands of two objects the same? (same number, same wavelength limits)
  #
"""
    class(ty_optical_props), intent(in) :: this, that
    logical                             :: bands_are_equal
"""
  function bands_are_equal(this::ty_optical_props{DT}, that::ty_optical_props{DT}) where DT
    bands_are_equal = get_nband(this) == get_nband(that) && get_nband(this) > 0
    if !bands_are_equal
      return bands_are_equal
    end
    bands_are_equal = all(abs(get_band_lims_wavenumber(this) - get_band_lims_wavenumber(that)) < DT(5) * spacing(get_band_lims_wavenumber(this)))
    return bands_are_equal
  end
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Is the g-point structure of two objects the same?
  #   (same bands, same number of g-points, same mapping between bands and g-points)
  #
"""
    class(ty_optical_props), intent(in) :: this, that
    logical                             :: gpoints_are_equal
"""
  function gpoints_are_equal(this::T1, that::T2) where {T1,T2}

    gpoints_are_equal = bands_are_equal(this, that) && get_ngpt(this) == get_ngpt(that)
    if !gpoints_are_equal
      return gpoints_are_equal
    end
    gpoints_are_equal = all(get_gpoint_bands(this) == get_gpoint_bands(that))
  end
  # -----------------------------------------------------------------------------------------------
  #
  # --- Setting/getting the name
  #
  # -----------------------------------------------------------------------------------------------
"""
    class(ty_optical_props),  intent(inout) :: this
    character(len=*),         intent(in   ) :: name
"""
  set_name!(this, name) = (this.name = strip(name))

  # --------------------------------------------------------
"""
    class(ty_optical_props),  intent(in   ) :: this
    character(len=name_len)                 :: get_name
"""
  get_name(this) = trim(this.name)
  # ------------------------------------------------------------------------------------------

include(joinpath("kernels","mo_optical_props_kernels.jl"))

end # module
