"""
    OpticalProps

Encapsulate optical properties defined on a spectral grid of N bands.
 The bands are described by their limiting wavenumbers. They need not be contiguous or complete.
 A band may contain more than one spectral sub-point (g-point) in which case a mapping must be supplied.
 A name may be provided and will be prepended to error messages.
 The base class (AbstractOpticalProps) encapsulates only this spectral discretization and must be initialized
    with the spectral information before use.

 Optical properties may be represented as arrays with dimensions ncol, nlay, ngpt
 (abstract class AbstractOpticalPropsArry).
 The type holds arrays depending on how much information is needed
 There are three possibilities
    `OneScalar` holds absorption optical depth `τ`, used in calculations accounting for extinction and emission
    `TwoStream` holds extinction optical depth `τ`, single-scattering albedo (`ssa`), and asymmetry parameter `g`.
 These classes must be allocated before use. Initialization and allocation can be combined.
 The classes have a validate() function that checks all arrays for valid values (e.g. `τ` > 0.)

Optical properties can be delta-scaled (though this is currently implemented only for two-stream arrays)

Optical properties can increment or "add themselves to" a set of properties represented with arrays
 as long as both sets have the same underlying band structure. Properties defined by band
 may be added to properties defined by g-point; the same value is assumed for all g-points with each band.

Subsets of optical properties held as arrays may be extracted along the column dimension.

"""
module OpticalProps

using DocStringExtensions
using ..FortranIntrinsics
using ..Utilities

export delta_scale!,
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

export AbstractOpticalProps,
       OneScalar,
       TwoStream,
       OpticalPropsBase

export get_nband, get_ngpt

export AbstractOpticalProps
abstract type AbstractOpticalProps{FT,I} end
export AbstractOpticalPropsArry
abstract type AbstractOpticalPropsArry{FT,I} <: AbstractOpticalProps{FT,I} end

struct ProblemSize{I}
  ncol::I
  nlay::I
  ngpt::I
  s::NTuple{3,I}
  ProblemSize(ncol::I,nlay::I,ngpt::I) where {I<:Int} =
    new{I}(ncol,nlay,ngpt, (ncol,nlay,ngpt))
end

"""
    OpticalPropsBase{FT,I} <: AbstractOpticalProps{FT,I}

Base class for optical properties. Describes the spectral
discretization including the wavenumber limits of each band
(spectral region) and the mapping between g-points and bands.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OpticalPropsBase{FT,I} <: AbstractOpticalProps{FT,I}
  "begin g-point, end g-point"
  band2gpt::Array{I,2}
  "band = gpt2band[g-point]"
  gpt2band::Array{I,1}
  "(upper and lower wavenumber by band) = band_lims_wvn[2,band]"
  band_lims_wvn::Array{FT,2}
  "name of particular optical properties"
  name::AbstractString
  function OpticalPropsBase(name::AbstractString, band_lims_wvn::Array{FT}, band_lims_gpt=nothing) where FT
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
    OneScalar{FT,I} <: AbstractOpticalPropsArry{FT,I}

Holds absorption optical depth `τ`, used in calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OneScalar{FT,I} <: AbstractOpticalPropsArry{FT,I}
  base::OpticalPropsBase{FT}
  τ::Array{FT,3}
end

OneScalar(base::OpticalPropsBase{FT}, ps::ProblemSize{I}) where {FT<:AbstractFloat,I<:Int} =
  OneScalar{FT,I}(base, Array{FT}(undef, ps.s...))

"""
    TwoStream{FT,I} <: AbstractOpticalPropsArry{FT,I}

Holds extinction optical depth `τ`, single-scattering albedo (`ssa`), and asymmetry parameter `g`.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStream{FT,I} <: AbstractOpticalPropsArry{FT,I}
  base::OpticalPropsBase{FT}
  τ::Array{FT,3}
  ssa::Array{FT,3}
  g::Array{FT,3}
end

TwoStream(base::OpticalPropsBase{FT}, ps::ProblemSize{I}) where {FT<:AbstractFloat,I<:Int} =
  TwoStream{FT,I}(base, ntuple(i->Array{FT}(undef, ps.s...),3)... )

#####
#####  Delta-scaling
#####

"""
    delta_scale!(...)

    class(OneScalar), intent(inout) :: this
    real(FT), dimension(:,:,:), optional, intent(in   ) :: for_
"""
delta_scale!(this::OneScalar, for_) = ""

"""
    delta_scale!()

    class(TwoStream), intent(inout) :: this
    real(FT), dimension(:,:,:), optional, intent(in   ) :: for_
    # Forward scattering fraction; g**2 if not provided

    integer :: ncol, nlay, ngpt
"""
function delta_scale!(this::TwoStream{FT}, for_ = nothing) where FT
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
  validate!(this::OneScalar{FT}) where FT

"""
function validate!(this::OneScalar{FT}) where FT
  # Validate sizes
  @assert !any_vals_less_than(this.τ, FT(0))
end

"""
    validate!(this::TwoStream{FT}) where FT

"""
function validate!(this::TwoStream{FT}) where FT
  # Validate sizes
  @assert all(size(this.ssa) == size(this.τ))
  @assert all(size(this.g) == size(this.τ))
  # Valid values
  @assert !any_vals_less_than(this.τ,  FT(0))
  @assert !any_vals_outside(this.ssa,  FT(0), FT(1))
  @assert !any_vals_outside(this.g  , FT(-1), FT(1))
end

#####
#####  Routines for array classes: incrementing
#####

"""
  increment!(op_in::AbstractOpticalPropsArry, op_io::AbstractOpticalPropsArry)

"""
function increment!(op_in::AbstractOpticalPropsArry, op_io::AbstractOpticalPropsArry)
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
    get_ncol(this::AbstractOpticalProps)

Number of columns
"""
get_ncol(this::AbstractOpticalProps) = size(this.τ, 1)

"""
    get_ncol(this::AbstractOpticalProps)

Number of layers
"""
get_nlay(this::AbstractOpticalProps) = size(this.τ, 2)

#####
#####  Routines for base class: spectral discretization
#####

"""
    get_nband(this::AbstractOpticalProps)

Number of bands
"""
get_nband(this::AbstractOpticalProps) = get_nband(this.base)
get_nband(this::OpticalPropsBase) = size(this.band2gpt,2)

"""
    get_ngpt(this::AbstractOpticalProps)

Number of g-points
  class(AbstractOpticalProps), intent(in) :: this
  integer                             :: get_ngpt
"""
get_ngpt(this::AbstractOpticalProps) = get_ngpt(this.base)
get_ngpt(this::OpticalPropsBase) = max(this.band2gpt...)


"""
    get_band_lims_gpoint(this::AbstractOpticalProps)

The first and last g-point of all bands at once
dimension (2, nbands)

  class(AbstractOpticalProps), intent(in) :: this
  integer, dimension(size(this%band2gpt,dim=1), size(this%band2gpt,dim=2))
                                      :: get_band_lims_gpoint
"""
get_band_lims_gpoint(this::AbstractOpticalProps) = get_band_lims_gpoint(this.base)
get_band_lims_gpoint(this::OpticalPropsBase) = this.band2gpt

"""
    convert_band2gpt(this::AbstractOpticalProps, band)

First and last g-point of a specific band

  class(AbstractOpticalProps), intent(in) :: this
  integer,                 intent(in) :: band
  integer, dimension(2)               :: convert_band2gpt
"""
convert_band2gpt(this::AbstractOpticalProps, band) = convert_band2gpt(this.base)
convert_band2gpt(this::OpticalPropsBase, band) = this.band2gpt[:,band]

"""
    get_band_lims_wavenumber(this::AbstractOpticalProps)

Lower and upper wavenumber of all bands
(upper and lower wavenumber by band) = band_lims_wvn(2,band)

  class(AbstractOpticalProps), intent(in) :: this
  real(FT), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2))
                                      :: get_band_lims_wavenumber
"""
get_band_lims_wavenumber(this::AbstractOpticalProps) = get_band_lims_wavenumber(this.base)
get_band_lims_wavenumber(this::OpticalPropsBase) = this.band_lims_wvn

"""
    get_band_lims_wavelength(this::AbstractOpticalProps)

Lower and upper wavelength of all bands

class(AbstractOpticalProps), intent(in) :: this
real(FT), dimension(size(this%band_lims_wvn,1), size(this%band_lims_wvn,2))
                                    :: get_band_lims_wavelength
"""
get_band_lims_wavelength(this::AbstractOpticalProps) = get_band_lims_wavelength(this.base)
get_band_lims_wavelength(this::OpticalPropsBase) = 1 ./ this.band_lims_wvn

"""
    get_gpoint_bands(this::AbstractOpticalProps)

Bands for all the g-points at once
dimension (ngpt)

  class(AbstractOpticalProps), intent(in) :: this
  integer, dimension(size(this%gpt2band,dim=1)) :: get_gpoint_bands
"""
get_gpoint_bands(this::AbstractOpticalProps) = get_gpoint_bands(this.base)
get_gpoint_bands(this::OpticalPropsBase) = this.gpt2band

"""
    convert_gpt2band(this::AbstractOpticalProps, gpt)

Band associated with a specific g-point

    class(AbstractOpticalProps), intent(in) :: this
    integer,                            intent(in) :: gpt
    integer                             :: convert_gpt2band
"""
convert_gpt2band(this::AbstractOpticalProps, gpt) = convert_gpt2band(this.base)
convert_gpt2band(this::OpticalPropsBase, gpt) = this.gpt2band[gpt]

"""
    bands_are_equal(this::AbstractOpticalProps{FT}, that::AbstractOpticalProps{FT}) where FT

Are the bands of two objects the same? (same number, same wavelength limits)

    class(AbstractOpticalProps), intent(in) :: this, that
    logical                             :: bands_are_equal
"""
function bands_are_equal(this::AbstractOpticalProps{FT}, that::AbstractOpticalProps{FT}) where FT
  if get_nband(this) == get_nband(that) && get_nband(this) > 0
    return all(abs.(get_band_lims_wavenumber(this) .-
                    get_band_lims_wavenumber(that)) .<
    FT(5) * spacing.(get_band_lims_wavenumber(this)))
  else
    return false
  end
end

"""
    gpoints_are_equal(this::AbstractOpticalProps{FT}, that::AbstractOpticalProps{FT}) where {FT}

Is the g-point structure of two objects the same?
  (same bands, same number of g-points, same mapping between bands and g-points)

    class(AbstractOpticalProps), intent(in) :: this, that
    logical                             :: gpoints_are_equal
"""
function gpoints_are_equal(this::AbstractOpticalProps{FT}, that::AbstractOpticalProps{FT}) where {FT}
  if bands_are_equal(this, that) && get_ngpt(this) == get_ngpt(that)
    return all(get_gpoint_bands(this) == get_gpoint_bands(that))
  else
    return false
  end
end

include(joinpath("kernels","OpticalProps_kernels.jl"))

end # module
