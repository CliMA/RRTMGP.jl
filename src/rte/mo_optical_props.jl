"""
    mo_optical_props

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
    ty_optical_props_2str holds extinction optical depth tau, single-scattering albedo (ssa), and asymmetry parameter g.
    ty_optical_props_nstr holds extinction optical depth tau, ssa, and phase function moments p with leading dimension nmom.
 These classes must be allocated before use. Initialization and allocation can be combined.
 The classes have a validate() function that checks all arrays for valid values (e.g. tau > 0.)

Optical properties can be delta-scaled (though this is currently implemented only for two-stream arrays)

Optical properties can increment or "add themselves to" a set of properties represented with arrays
 as long as both sets have the same underlying band structure. Properties defined by band
 may be added to properties defined by g-point; the same value is assumed for all g-points with each band.

Subsets of optical properties held as arrays may be extracted along the column dimension.

"""
module mo_optical_props

using DocStringExtensions
using ..fortran_intrinsics
import ..fortran_intrinsics: is_initialized
using ..mo_util_array

export init_base_from_copy!,
       alloc!,
       init_and_alloc!,
       copy_and_alloc!,
       delta_scale!,
       validate!,
       increment!,
       convert_band2gpt,
       get_band_lims_wavenumber,
       get_band_lims_wavelength,
       get_gpoint_bands,
       expand,
       bands_are_equal,
       convert_gpt2band,
       get_band_lims_gpoint,
       get_ncol,
       get_nlay,
       gpoints_are_equal,
       is_initialized

export ty_optical_props,
       ty_optical_props_1scl,
       ty_optical_props_2str,
       ty_optical_props_nstr,
       init!,
       ty_optical_props_base

export get_nband, get_ngpt

export ty_optical_props
abstract type ty_optical_props{FT,I} end
export ty_optical_props_arry
abstract type ty_optical_props_arry{FT,I} <: ty_optical_props{FT,I} end

"""
    ty_optical_props_base{FT,I} <: ty_optical_props{FT,I}

Base class for optical properties. Describes the spectral
discretization including the wavenumber limits of each band
(spectral region) and the mapping between g-points and bands.

 - (begin g-point, end g-point) = band2gpt(2,band)
 - band = gpt2band(g-point)
 - (upper and lower wavenumber by band) = band_lims_wvn(2,band)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ty_optical_props_base{FT,I} <: ty_optical_props{FT,I}
  band2gpt::Array{I,2}
  gpt2band::Array{I,1}
  band_lims_wvn::Array{FT,2}
  name::AbstractString
end

function ty_optical_props_base(name, band_lims_wvn::Array{FT}, band_lims_gpt=nothing) where FT
  @assert size(band_lims_wvn,1) == 2
  @assert !any(band_lims_wvn.<FT(0))
  band_lims_gpt_lcl = Array{Int}(undef, 2, size(band_lims_wvn, 2))
  if band_lims_gpt ≠ nothing
    @assert size(band_lims_gpt,1) == 2
    @assert size(band_lims_gpt,2) == size(band_lims_wvn,2)
    @assert !any(band_lims_gpt .< 1)
    band_lims_gpt_lcl .= band_lims_gpt
  else
    for iband in 1:size(band_lims_wvn, 2)
      band_lims_gpt_lcl[1:2,iband] .= iband
    end
  end

  band2gpt = band_lims_gpt_lcl

  # Make a map between g-points and bands
  gpt2band = Array{Int}(undef, max(band_lims_gpt_lcl...))
  for iband in 1:size(band_lims_gpt_lcl, 2)
    gpt2band[band_lims_gpt_lcl[1,iband]:band_lims_gpt_lcl[2,iband]] .= iband
  end
  return ty_optical_props_base{FT,Int}(band2gpt,gpt2band,band_lims_wvn,name)
end

mutable struct ty_optical_props_1scl{FT,I} <: ty_optical_props_arry{FT,I}
  base
  tau#::Array{FT,3}
end

mutable struct ty_optical_props_2str{FT,I} <: ty_optical_props_arry{FT,I}
  base
  tau#::Array{FT,3}
  ssa#::Array{FT,3}
  g#::Array{FT,3}
end

mutable struct ty_optical_props_nstr{FT,I} <: ty_optical_props_arry{FT,I}
  base
  tau#::Array{FT,3}
  ssa#::Array{FT,3}
  p#::Array{FT,4}
end

ty_optical_props_1scl(FT,I) = ty_optical_props_1scl{FT,I}(ntuple(i->nothing, 2)...)
ty_optical_props_2str(FT,I) = ty_optical_props_2str{FT,I}(ntuple(i->nothing, 4)...)

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
  init_base_from_copy!(...)

  class(ty_optical_props),    intent(inout) :: this
  class(ty_optical_props),    intent(in   ) :: spectral_desc
  character(len = 128)                        :: err_message
"""
init_base_from_copy(spectral_desc::ty_optical_props) =
 ty_optical_props_base("called_from_init_base_from_copy" ,get_band_lims_wavenumber(spectral_desc), get_band_lims_gpoint(spectral_desc))

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
is_initialized(this::ty_optical_props_base) = allocated(this.band2gpt)
is_initialized(this::ty_optical_props) = is_initialized(this.base)

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
function alloc!(this::ty_optical_props_1scl{FT}, ncol, nlay) where FT
  @assert !any([ncol, nlay] .<= 0)
  allocated(this.tau) && deallocate!(this.tau)
  this.tau = Array{FT}(undef, ncol, nlay, get_ngpt(this))
end

# --- 2 stream ------------------------------------------------------------------------
"""
  alloc_only_2str!(...)

  class(ty_optical_props_2str)    :: this
  integer,             intent(in) :: ncol, nlay
  character(len=128)              :: err_message
"""

function alloc!(this::ty_optical_props_2str{FT}, ncol, nlay) where FT
  @assert !any([ncol, nlay] .<= 0)

  allocated(this.tau) && deallocate!(this.tau)
  this.tau = Array{FT}(undef, ncol,nlay,get_ngpt(this))

  allocated(this.ssa) && deallocate!(this.ssa)
  this.ssa = Array{FT}(undef, ncol,nlay,get_ngpt(this))
  allocated(this.g) && deallocate!(this.g)
  this.g = Array{FT}(undef, ncol,nlay,get_ngpt(this))
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
    real(FT), dimension(:,:),     intent(in) :: band_lims_wvn
    integer,  dimension(:,:),
                      optional,   intent(in) :: band_lims_gpt
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""
function init_and_alloc!(this::ty_optical_props_1scl, ncol, nlay, band_lims_wvn, band_lims_gpt=nothing, name=nothing)
  if band_lims_gpt==nothing
    this.base = ty_optical_props_base(name, band_lims_wvn)
  else
    this.base = ty_optical_props_base(name, band_lims_wvn, band_lims_gpt)
  end
  alloc!(this, ncol, nlay)
end

#-------------------------------------------------------------------------------------------------
#
# Initialization from an existing spectral discretization/ty_optical_props
#
#-------------------------------------------------------------------------------------------------
"""
    copy_and_alloc!(...)

    class(ty_optical_props_2str)             :: this
    integer,                      intent(in) :: ncol, nlay
    class(ty_optical_props     ), intent(in) :: spectral_desc
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""
function copy_and_alloc!(this::ty_optical_props, ncol, nlay, spectral_desc::ty_optical_props, name="optical props")
  this.base = ty_optical_props_base(name, get_band_lims_wavenumber(spectral_desc),
                                     get_band_lims_gpoint(spectral_desc))
  alloc!(this, ncol, nlay)
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
    real(FT), dimension(:,:,:), optional,
                                  intent(in   ) :: for_
    character(128)                              :: err_message
"""
delta_scale!(this::ty_optical_props_1scl, for_) = ""


"""
    delta_scale!()

    class(ty_optical_props_2str), intent(inout) :: this
    real(FT), dimension(:,:,:), optional,
                                  intent(in   ) :: for_
    # Forward scattering fraction; g**2 if not provided
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt
"""
function delta_scale!(this::ty_optical_props_2str{FT}, for_ = nothing) where FT
  ncol = get_ncol(this)
  nlay = get_nlay(this)
  ngpt = get_ngpt(this)

  if for_ ≠ nothing
    @assert all(size(for_) .== (ncol, nlay, ngpt))
    @assert !any(for_ < FT(0) || for_ > FT(1))
    # delta_scale_2str_kernel!(ncol, nlay, ngpt, this.tau, this.ssa, this.g, for_)
    delta_scale_2str_kernel!(this, for_)
  else
    # delta_scale_2str_kernel!(ncol, nlay, ngpt, this.tau, this.ssa, this.g)
    delta_scale_2str_kernel!(this)
  end
end


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
function validate!(this::ty_optical_props_1scl{FT}) where FT
  @assert allocated(this.tau)
  # Valid values
  @assert !any_vals_less_than(this.tau, FT(0))
end

"""

  class(ty_optical_props_2str), intent(in) :: this
  character(len=128)                       :: err_message

  integer :: varSizes(3)
"""
function validate!(this::ty_optical_props_2str{FT}) where FT

  #
  # Array allocation status, sizing
  #
  @assert allocated(this.tau)
  @assert allocated(this.ssa)
  @assert allocated(this.g)
  @assert all(size(this.ssa) == size(this.tau))
  @assert all(size(this.g) == size(this.tau))
  #
  # Valid values
  #
  @assert !any_vals_less_than(this.tau,  FT(0))
  @assert !any_vals_outside(this.ssa,  FT(0), FT(1))
  @assert !any_vals_outside(this.g  , FT(-1), FT(1))
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
  @assert bands_are_equal(op_in, op_io)

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

    @assert get_ngpt(op_in) == get_nband(op_io)
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
get_nband(this::ty_optical_props) = get_nband(this.base)
get_nband(this::ty_optical_props_base) = is_initialized(this) ? size(this.band2gpt,2) : 0

# -----------------------------------------------------------------------------------------------
#
# Number of g-points
#
"""
  class(ty_optical_props), intent(in) :: this
  integer                             :: get_ngpt
"""
get_ngpt(this::ty_optical_props) = get_ngpt(this.base)
get_ngpt(this::ty_optical_props_base) = is_initialized(this) ? max(this.band2gpt...) : 0

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
get_band_lims_gpoint(this::ty_optical_props) = get_band_lims_gpoint(this.base)
get_band_lims_gpoint(this::ty_optical_props_base) = this.band2gpt


# First and last g-point of a specific band
"""
  class(ty_optical_props), intent(in) :: this
  integer,                 intent(in) :: band
  integer, dimension(2)               :: convert_band2gpt
"""
convert_band2gpt(this::ty_optical_props, band) = convert_band2gpt(this.base)
convert_band2gpt(this::ty_optical_props_base{FT}, band) where FT =
  is_initialized(this) ? this.band2gpt[:,band] : zeros(FT,length(this.band2gpt[:,band]))


# Lower and upper wavenumber of all bands
# (upper and lower wavenumber by band) = band_lims_wvn(2,band)
"""
  class(ty_optical_props), intent(in) :: this
  real(FT), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2))
                                      :: get_band_lims_wavenumber
"""
get_band_lims_wavenumber(this::ty_optical_props) = get_band_lims_wavenumber(this.base)
get_band_lims_wavenumber(this::ty_optical_props_base{FT}) where FT =
  is_initialized(this) ? this.band_lims_wvn[:,:] : zeros(FT, size(this.band_lims_wvn))


# Lower and upper wavelength of all bands
"""
  class(ty_optical_props), intent(in) :: this
  real(FT), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2))
                                      :: get_band_lims_wavelength
"""
get_band_lims_wavelength(this::ty_optical_props) = get_band_lims_wavelength(this.base)
get_band_lims_wavelength(this::ty_optical_props_base{FT}) where FT =
  is_initialized(this) ? 1 ./ this.band_lims_wvn[:,:] : zeros(FT, size(this.band_lims_wvn))


# Bands for all the g-points at once
# dimension (ngpt)
"""
  class(ty_optical_props), intent(in) :: this
  integer, dimension(size(this%gpt2band,dim=1))
                                      :: get_gpoint_bands
"""
get_gpoint_bands(this::ty_optical_props) = get_gpoint_bands(this.base)
get_gpoint_bands(this::ty_optical_props_base) =
  is_initialized(this) ? this.gpt2band[:] : zeros(Int,length(this.gpt2band))

#
# Band associated with a specific g-point

"""
    class(ty_optical_props), intent(in) :: this
    integer,                            intent(in) :: gpt
    integer                             :: convert_gpt2band
"""
convert_gpt2band(this::ty_optical_props, gpt) = convert_gpt2band(this.base)
convert_gpt2band(this::ty_optical_props_base, gpt) = is_initialized(this) ? this.gpt2band[gpt] : 0

#
# Expand an array of dimension arr_in(nband) to dimension arr_out(ngpt)
#
"""
    class(ty_optical_props), intent(in) :: this
    real(FT), dimension(:),  intent(in) :: arr_in # (nband)
    real(FT), dimension(size(this%gpt2band)) :: arr_out

    integer :: iband
"""
function expand(this::ty_optical_props{FT}, arr_in) where FT
  arr_out = Vector{FT}(undef, size(this.gpt2band))
  for iband in 1:get_nband(this)
    arr_out[this.band2gpt[1,iband]:this.band2gpt[2,iband]] = arr_in[iband]
  end
  return arr_out
end

#
# Are the bands of two objects the same? (same number, same wavelength limits)

"""
    class(ty_optical_props), intent(in) :: this, that
    logical                             :: bands_are_equal
"""
function bands_are_equal(this::ty_optical_props{FT}, that::ty_optical_props{FT}) where FT
  bands_equal = get_nband(this) == get_nband(that) && get_nband(this) > 0
  if !bands_equal
    return bands_equal
  end
  bands_equal = all(abs.(get_band_lims_wavenumber(this) .- get_band_lims_wavenumber(that)) .< FT(5) * spacing.(get_band_lims_wavenumber(this)))
  return bands_equal
end

#
# Is the g-point structure of two objects the same?
#   (same bands, same number of g-points, same mapping between bands and g-points)

"""
    class(ty_optical_props), intent(in) :: this, that
    logical                             :: gpoints_are_equal
"""
function gpoints_are_equal(this::T1, that::T2) where {T1,T2}

  gpoints_equal = bands_are_equal(this, that) && get_ngpt(this) == get_ngpt(that)
  if !gpoints_equal
    return gpoints_equal
  end
  gpoints_equal = all(get_gpoint_bands(this) == get_gpoint_bands(that))
end

"""
    class(ty_optical_props),  intent(inout) :: this
    character(len=*),         intent(in   ) :: name
"""
set_name!(this, name) = (this.name = name)

"""
    class(ty_optical_props),  intent(in   ) :: this
    character(len=name_len)                 :: get_name
"""
get_name(this) = this.name


include(joinpath("kernels","mo_optical_props_kernels.jl"))

end # module
