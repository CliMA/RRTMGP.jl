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
 There are three possibilities
    `ty_optical_props_1scl` holds absorption optical depth `tau`, used in calculations accounting for extinction and emission
    `ty_optical_props_2str` holds extinction optical depth `tau`, single-scattering albedo (`ssa`), and asymmetry parameter `g`.
 These classes must be allocated before use. Initialization and allocation can be combined.
 The classes have a validate() function that checks all arrays for valid values (e.g. `tau` > 0.)

Optical properties can be delta-scaled (though this is currently implemented only for two-stream arrays)

Optical properties can increment or "add themselves to" a set of properties represented with arrays
 as long as both sets have the same underlying band structure. Properties defined by band
 may be added to properties defined by g-point; the same value is assumed for all g-points with each band.

Subsets of optical properties held as arrays may be extracted along the column dimension.

"""
module mo_optical_props

using DocStringExtensions
using ..fortran_intrinsics
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
       ProblemSize,
       bands_are_equal,
       convert_gpt2band,
       get_band_lims_gpoint,
       get_ncol,
       get_nlay,
       gpoints_are_equal

export ty_optical_props,
       ty_optical_props_1scl,
       ty_optical_props_2str,
       ty_optical_props_base

export get_nband, get_ngpt

export ty_optical_props
abstract type ty_optical_props{FT,I} end
export ty_optical_props_arry
abstract type ty_optical_props_arry{FT,I} <: ty_optical_props{FT,I} end

struct ProblemSize{I}
  ncol::I
  nlay::I
  ngpt::I
  s::NTuple{3,I}
  ProblemSize(ncol::I,nlay::I,ngpt::I) where {I<:Int} =
    new{I}(ncol,nlay,ngpt, (ncol,nlay,ngpt))
end

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
    return new{FT,Int}(band2gpt,gpt2band,band_lims_wvn,name)
  end
end


"""
    ty_optical_props_1scl{FT,I} <: ty_optical_props_arry{FT,I}

Holds absorption optical depth `tau`, used in calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ty_optical_props_1scl{FT,I} <: ty_optical_props_arry{FT,I}
  base
  tau#::Array{FT,3}
end

ty_optical_props_1scl(base::ty_optical_props_base{FT}, ps::ProblemSize{I}) where {FT<:AbstractFloat,I<:Int} =
  ty_optical_props_1scl{FT,I}(base, Array{FT}(undef, ps.s...))

"""
    ty_optical_props_2str{FT,I} <: ty_optical_props_arry{FT,I}

Holds extinction optical depth `tau`, single-scattering albedo (`ssa`), and asymmetry parameter `g`.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ty_optical_props_2str{FT,I} <: ty_optical_props_arry{FT,I}
  base
  tau#::Array{FT,3}
  ssa#::Array{FT,3}
  g#::Array{FT,3}
end

ty_optical_props_2str(base::ty_optical_props_base{FT}, ps::ProblemSize{I}) where {FT<:AbstractFloat,I<:Int} =
  ty_optical_props_2str{FT,I}(base, ntuple(i->Array{FT}(undef, ps.s...),3)... )

#####
#####  Delta-scaling
#####

"""
    delta_scale!(...)

    class(ty_optical_props_1scl), intent(inout) :: this
    real(FT), dimension(:,:,:), optional, intent(in   ) :: for_
"""
delta_scale!(this::ty_optical_props_1scl, for_) = ""

"""
    delta_scale!()

    class(ty_optical_props_2str), intent(inout) :: this
    real(FT), dimension(:,:,:), optional, intent(in   ) :: for_
    # Forward scattering fraction; g**2 if not provided

    integer :: ncol, nlay, ngpt
"""
function delta_scale!(this::ty_optical_props_2str{FT}, for_ = nothing) where FT
  ncol = get_ncol(this)
  nlay = get_nlay(this)
  ngpt = get_ngpt(this)

  if for_ ≠ nothing
    @assert all(size(for_) .== (ncol, nlay, ngpt))
    @assert !any(for_ < FT(0) || for_ > FT(1))
    delta_scale_2str_kernel!(this, for_)
  else
    delta_scale_2str_kernel!(this)
  end
end

#####
##### Validation
#####

"""
  validate!(this::ty_optical_props_1scl{FT}) where FT

"""
function validate!(this::ty_optical_props_1scl{FT}) where FT
  # Validate sizes
  @assert !any_vals_less_than(this.tau, FT(0))
end

"""
    validate!(this::ty_optical_props_2str{FT}) where FT

"""
function validate!(this::ty_optical_props_2str{FT}) where FT
  # Validate sizes
  @assert all(size(this.ssa) == size(this.tau))
  @assert all(size(this.g) == size(this.tau))
  # Valid values
  @assert !any_vals_less_than(this.tau,  FT(0))
  @assert !any_vals_outside(this.ssa,  FT(0), FT(1))
  @assert !any_vals_outside(this.g  , FT(-1), FT(1))
end

#####
#####  Routines for array classes: incrementing
#####

"""
  increment!(op_in::ty_optical_props_arry, op_io::ty_optical_props_arry)

"""
function increment!(op_in::ty_optical_props_arry, op_io::ty_optical_props_arry)
  @assert bands_are_equal(op_in, op_io)

  if gpoints_are_equal(op_in, op_io)

    # Increment by gpoint
    #   (or by band if both op_in and op_io are defined that way)
    increment_by_gpoint!(op_io, op_in)

  else

    # Values defined by-band will have ngpt = nband
    # We can use values by band in op_in to increment op_io
    #   Anything else is an error

    @assert get_ngpt(op_in) == get_nband(op_io)
    # Increment by band
    increment_bybnd!(op_io, op_in)
  end
end

#####
#####  Routines for array classes: problem sizes
#####

"""
  class(ty_optical_props_arry), intent(in   ) :: this
  integer,                      intent(in   ) :: dim
  integer                                     :: get_arry_extent
"""
get_arry_extent(this::ty_optical_props, dim) = size(this.tau, dim)

"""
  class(ty_optical_props_arry), intent(in   ) :: this
  integer                                     :: get_ncol
"""
get_ncol(this::ty_optical_props) = get_arry_extent(this, 1)

"""
  class(ty_optical_props_arry), intent(in   ) :: this
  integer                                     :: get_nlay
"""
get_nlay(this::ty_optical_props) = get_arry_extent(this, 2)

#####
#####  Routines for base class: spectral discretization
#####

"""
    get_nband(this::ty_optical_props)

Number of bands
"""
get_nband(this::ty_optical_props) = get_nband(this.base)
get_nband(this::ty_optical_props_base) = size(this.band2gpt,2)

"""
    get_ngpt(this::ty_optical_props)

Number of g-points
  class(ty_optical_props), intent(in) :: this
  integer                             :: get_ngpt
"""
get_ngpt(this::ty_optical_props) = get_ngpt(this.base)
get_ngpt(this::ty_optical_props_base) = max(this.band2gpt...)


"""
    get_band_lims_gpoint(this::ty_optical_props)

The first and last g-point of all bands at once
dimension (2, nbands)

  class(ty_optical_props), intent(in) :: this
  integer, dimension(size(this%band2gpt,dim=1), size(this%band2gpt,dim=2))
                                      :: get_band_lims_gpoint
"""
get_band_lims_gpoint(this::ty_optical_props) = get_band_lims_gpoint(this.base)
get_band_lims_gpoint(this::ty_optical_props_base) = this.band2gpt

"""
    convert_band2gpt(this::ty_optical_props, band)

First and last g-point of a specific band

  class(ty_optical_props), intent(in) :: this
  integer,                 intent(in) :: band
  integer, dimension(2)               :: convert_band2gpt
"""
convert_band2gpt(this::ty_optical_props, band) = convert_band2gpt(this.base)
convert_band2gpt(this::ty_optical_props_base, band) = this.band2gpt[:,band]

"""
    get_band_lims_wavenumber(this::ty_optical_props)

Lower and upper wavenumber of all bands
(upper and lower wavenumber by band) = band_lims_wvn(2,band)

  class(ty_optical_props), intent(in) :: this
  real(FT), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2))
                                      :: get_band_lims_wavenumber
"""
get_band_lims_wavenumber(this::ty_optical_props) = get_band_lims_wavenumber(this.base)
get_band_lims_wavenumber(this::ty_optical_props_base) = this.band_lims_wvn

"""
    get_band_lims_wavelength(this::ty_optical_props)

Lower and upper wavelength of all bands

class(ty_optical_props), intent(in) :: this
real(FT), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2))
                                    :: get_band_lims_wavelength
"""
get_band_lims_wavelength(this::ty_optical_props) = get_band_lims_wavelength(this.base)
get_band_lims_wavelength(this::ty_optical_props_base) = 1 ./ this.band_lims_wvn

"""
    get_gpoint_bands(this::ty_optical_props)

Bands for all the g-points at once
dimension (ngpt)

  class(ty_optical_props), intent(in) :: this
  integer, dimension(size(this%gpt2band,dim=1)) :: get_gpoint_bands
"""
get_gpoint_bands(this::ty_optical_props) = get_gpoint_bands(this.base)
get_gpoint_bands(this::ty_optical_props_base) = this.gpt2band

"""
    convert_gpt2band(this::ty_optical_props, gpt)

Band associated with a specific g-point

    class(ty_optical_props), intent(in) :: this
    integer,                            intent(in) :: gpt
    integer                             :: convert_gpt2band
"""
convert_gpt2band(this::ty_optical_props, gpt) = convert_gpt2band(this.base)
convert_gpt2band(this::ty_optical_props_base, gpt) = this.gpt2band[gpt]

"""
    expand(this::ty_optical_props{FT}, arr_in) where FT

Expand an array of dimension arr_in(nband) to dimension arr_out(ngpt)

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

"""
    bands_are_equal(this::ty_optical_props{FT}, that::ty_optical_props{FT}) where FT

Are the bands of two objects the same? (same number, same wavelength limits)

    class(ty_optical_props), intent(in) :: this, that
    logical                             :: bands_are_equal
"""
function bands_are_equal(this::ty_optical_props{FT}, that::ty_optical_props{FT}) where FT
  if get_nband(this) == get_nband(that) && get_nband(this) > 0
    return all(abs.(get_band_lims_wavenumber(this) .-
                    get_band_lims_wavenumber(that)) .<
    FT(5) * spacing.(get_band_lims_wavenumber(this)))
  else
    return false
  end
end

"""
    gpoints_are_equal(this::ty_optical_props{FT}, that::ty_optical_props{FT}) where {FT}

Is the g-point structure of two objects the same?
  (same bands, same number of g-points, same mapping between bands and g-points)

    class(ty_optical_props), intent(in) :: this, that
    logical                             :: gpoints_are_equal
"""
function gpoints_are_equal(this::ty_optical_props{FT}, that::ty_optical_props{FT}) where {FT}
  if bands_are_equal(this, that) && get_ngpt(this) == get_ngpt(that)
    return all(get_gpoint_bands(this) == get_gpoint_bands(that))
  else
    return false
  end
end

include(joinpath("kernels","mo_optical_props_kernels.jl"))

end # module
