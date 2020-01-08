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
       get_band_lims_wavenumber,
       get_gpoint_bands,
       bands_are_equal,
       gpt_range,
       get_τ′_size,
       get_band_lims_gpoint

export get_nband,
       get_ngpt,
       get_ncol,
       get_nlay

export AbstractOpticalProps,
       AbstractOpticalPropsArry,
       AbstractOpticalPropsPGP,
       OneScalar,
       TwoStream,
       OneScalarPGP,
       TwoStreamPGP,
       OpticalPropsBase

"""
    AbstractOpticalProps{FT<:AbstractFloat,I<:Int}

Abstract optical properties
"""
abstract type AbstractOpticalProps{FT<:AbstractFloat,I<:Int} end

"""
    AbstractOpticalPropsArry{FT,I} <: AbstractOpticalProps{FT,I}

Abstract optical properties array. Includes arrays for
optical depth and, depending on the approximation used,
other variables.
"""
abstract type AbstractOpticalPropsArry{FT,I} <: AbstractOpticalProps{FT,I} end

"""
    AbstractOpticalPropsPGP{FT,I} <: AbstractOpticalProps{FT,I}

Abstract optical properties array per-grid-point. Includes
arrays for optical depth and, depending on the approximation
used, other variables.
"""
abstract type AbstractOpticalPropsPGP{FT,I} <: AbstractOpticalProps{FT,I} end

"""
    OpticalPropsBase{FT,I} <: AbstractOpticalProps{FT,I}

Base class for optical properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OpticalPropsBase{FT,I} <: AbstractOpticalProps{FT,I}
  "Map from band to g-point. `(g_point_begin, g_point_end) = band2gpt[2,band]`"
  band2gpt::Array{I,2}
  "Map from g-point to band. `band = gpt2band(gpt)`"
  gpt2band::Array{I,1}
  "Wavenumber band limits. `(upper_band_wvn, lower_band_wvn) = band_lims_wvn[2,band]`"
  band_lims_wvn::Array{FT,2}
  "Name of optical properties"
  name::AbstractString
  function OpticalPropsBase(name, band_lims_wvn::Array{FT}, band_lims_gpt=nothing) where FT
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

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OneScalar{FT,I} <: AbstractOpticalPropsArry{FT,I}
  "Base optical properties, see [`OpticalPropsBase`](@ref)"
  base::Union{OpticalPropsBase{FT,I}, Nothing}
  "Optical depth"
  τ::Array{FT,3}
end

OneScalar(base::OpticalPropsBase{FT,I}, ncol::I, nlay::I, ngpt::I) where {FT<:AbstractFloat,I<:Int} =
  OneScalar{FT,I}(base, Array{FT}(undef, ncol, nlay, ngpt))

"""
    OneScalarPGP{FT,I} <: AbstractOpticalPropsArry{FT,I}

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission
per-grid-point

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OneScalarPGP{FT,I} <: AbstractOpticalPropsPGP{FT,I}
  "Base optical properties, see [`OpticalPropsBase`](@ref)"
  base::Union{OpticalPropsBase{FT,I}, Nothing}
  "Optical depth"
  τ::Array{FT,1}
end

function Base.convert(::Type{OneScalar}, data::Array{OneScalarPGP{FT,I}}) where {FT,I}
  s = size(data)
  ngpt = length(first(data).τ)
  return OneScalar{FT,I}(first(data).base, Array{FT}([data[i,j].τ[k] for i in 1:s[1], j in 1:s[2], k in 1:ngpt]))
end

Base.convert(::Type{Array{OneScalarPGP}}, data::OneScalar{FT,I}) where {FT,I} =
  [OneScalarPGP{FT,I}(data.base, data.τ[i,j,:]) for i in 1:size(data.τ,1), j in 1:size(data.τ,2)]

"""
    TwoStream{FT,I} <: AbstractOpticalPropsArry{FT,I}

Two stream approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStream{FT,I} <: AbstractOpticalPropsArry{FT,I}
  "Base optical properties, see [`OpticalPropsBase`](@ref)"
  base::Union{OpticalPropsBase{FT,I}, Nothing}
  "Optical depth"
  τ::Array{FT,3}
  "Single-scattering albedo"
  ssa::Array{FT,3}
  "Asymmetry parameter"
  g::Array{FT,3}
end

TwoStream(base::OpticalPropsBase{FT,I}, ncol::I, nlay::I, ngpt::I) where {FT<:AbstractFloat,I<:Int} =
  TwoStream{FT,I}(base, ntuple(i->Array{FT}(undef, ncol, nlay, ngpt),3)... )
TwoStream(::Type{FT}, ncol::I, nlay::I, ngpt::I) where {FT<:AbstractFloat,I<:Int} =
  TwoStream{FT,I}(nothing, ntuple(i->Array{FT}(undef, ncol, nlay, ngpt),3)... )

"""
    TwoStreamPGP{FT,I} <: AbstractOpticalPropsArry{FT,I}

Two stream approximation for optical depth, used in
calculations accounting for extinction and emission
per-grid-point

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStreamPGP{FT,I} <: AbstractOpticalPropsPGP{FT,I}
  "Base optical properties, see [`OpticalPropsBase`](@ref)"
  base::Union{OpticalPropsBase{FT,I}, Nothing}
  "Optical depth"
  τ::Array{FT,1}
  "Single-scattering albedo"
  ssa::Array{FT,1}
  "Asymmetry parameter"
  g::Array{FT,1}
end

TwoStreamPGP(::Type{FT}, ngpt::I) where {FT<:AbstractFloat,I<:Int} =
  TwoStreamPGP{FT,I}(nothing, ntuple(i->Array{FT}(undef, ngpt),3)... )

function Base.convert(::Type{TwoStream}, data::Array{TwoStreamPGP{FT,I}}) where {FT,I}
  s = size(data)
  ngpt = length(first(data).τ)
  return TwoStream{FT,I}(first(data).base, Array{FT}([data[i,j].τ[k]   for i in 1:s[1], j in 1:s[2], k in 1:ngpt]),
                                           Array{FT}([data[i,j].ssa[k] for i in 1:s[1], j in 1:s[2], k in 1:ngpt]),
                                           Array{FT}([data[i,j].g[k]   for i in 1:s[1], j in 1:s[2], k in 1:ngpt]))
end

Base.convert(::Type{Array{TwoStreamPGP}}, data::TwoStream{FT,I}) where {FT,I} =
  [TwoStreamPGP{FT,I}(data.base, data.τ[i,j,:], data.ssa[i,j,:],data.g[i,j,:]) for i in 1:size(data.τ,1), j in 1:size(data.τ,2)]

#####
#####  Delta-scaling
#####

"""
    delta_scale!(this::OneScalar{FT}, fwd_scat_frac)

 - `this` optical properties, see [`OneScalar`](@ref)
 - `fwd_scat_frac` forward scattering fraction
"""
delta_scale!(this::OneScalar{FT}, fwd_scat_frac) where {FT} = nothing

"""
    delta_scale!(this::TwoStream{FT}, for_ = nothing) where {FT<:AbstractFloat}

 - `this` optical properties, see [`TwoStream`](@ref)
 - `fwd_scat_frac` forward scattering fraction
"""
function delta_scale!(this::TwoStream{FT}, fwd_scat_frac = nothing) where {FT<:AbstractFloat}
  ncol = get_ncol(this)
  nlay = get_nlay(this)
  ngpt = get_ngpt(this)

  if fwd_scat_frac ≠ nothing
    @assert all(size(fwd_scat_frac) .== (ncol, nlay, ngpt))
    @assert !any(fwd_scat_frac < FT(0) || fwd_scat_frac > FT(1))
    delta_scale_kernel!(this, fwd_scat_frac)
  else
    delta_scale_kernel!(this)
  end
  return nothing
end

#####
##### Validation
#####

"""
  validate!(this::OneScalar{FT}) where FT

Validate values of optical properties
"""
function validate!(this::OneScalar{FT}) where {FT<:AbstractFloat}
  # Validate sizes
  @assert !any_vals_less_than(this.τ, FT(0))
end

"""
    validate!(this::TwoStream{FT}) where FT

Validate values and sizes of optical properties
"""
function validate!(this::TwoStream{FT}) where {FT<:AbstractFloat}
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

Increment (combine) two sets of optical properties
"""
function increment!(op_in::AbstractOpticalPropsArry, op_io::AbstractOpticalPropsArry)
  @assert bands_are_equal(op_in, op_io)
  if gpoints_are_equal(op_in, op_io)
    # Increment by gpoint or by band if both op_in and op_io are defined that way
    increment_by_gpoint!(op_io, op_in)
  else
    # Values defined by-band will have ngpt = nband
    # We can use values by band in op_in to increment op_io
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

Number of columns, given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
get_ncol(this::AbstractOpticalProps) = size(this.τ, 1)

"""
    get_ncol(this::AbstractOpticalProps)

Number of layers, given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
get_nlay(this::AbstractOpticalProps) = size(this.τ, 2)

#####
#####  Routines for base class: spectral discretization
#####

"""
    get_nband(this::AbstractOpticalProps)

Number of bands, given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
get_nband(this::AbstractOpticalProps) = get_nband(this.base)
get_nband(this::OpticalPropsBase) = size(this.band2gpt,2)

"""
    get_ngpt(this::AbstractOpticalProps)

Number of g-points, given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
get_ngpt(this::AbstractOpticalProps) = get_ngpt(this.base)
get_ngpt(this::OpticalPropsBase) = max(this.band2gpt...)

"""
    get_τ′_size(this::AbstractOpticalProps)

Size of τ′

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
get_τ′_size(this::AbstractOpticalPropsArry) = (size(this.τ,3),size(this.τ,2),size(this.τ,1))
get_τ′_size(this::AbstractOpticalPropsPGP) = length(this.τ)


"""
    get_band_lims_gpoint(this::AbstractOpticalProps)

The first and last g-point of all bands at once, given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
get_band_lims_gpoint(this::AbstractOpticalProps) = get_band_lims_gpoint(this.base)
get_band_lims_gpoint(this::OpticalPropsBase) = this.band2gpt

"""
    get_band_lims_wavenumber(this::AbstractOpticalProps)

Lower and upper wavenumber of all bands, given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
get_band_lims_wavenumber(this::AbstractOpticalProps) = get_band_lims_wavenumber(this.base)
get_band_lims_wavenumber(this::OpticalPropsBase) = this.band_lims_wvn

"""
    get_gpoint_bands(this::AbstractOpticalProps)

Bands for all the g-points at once, given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
get_gpoint_bands(this::AbstractOpticalProps) = get_gpoint_bands(this.base)
get_gpoint_bands(this::OpticalPropsBase) = this.gpt2band

"""
    gpt_range(this::AbstractOpticalProps)

A range of g-point bands, given the i-th band

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
"""
gpt_range(this::AbstractOpticalProps{FT,I}, ibnd::I) where {FT,I} = gpt_range(this.base, ibnd)
gpt_range(this::OpticalPropsBase{FT,I}, ibnd::I) where {FT,I<:Int} = this.band2gpt[1,ibnd]:this.band2gpt[2,ibnd]

"""
    bands_are_equal(this::AbstractOpticalProps{FT}, that::AbstractOpticalProps{FT}) where FT

Boolean that indicates if the bands of two objects
the same, (same number, same wavelength limits), given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
 - `that` optical properties, see [`AbstractOpticalProps`](@ref)
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

Boolean that indicates if the g-point structure of two
objects the same, (same bands, same number of g-points,
same mapping between bands and g-points), given

 - `this` optical properties, see [`AbstractOpticalProps`](@ref)
 - `that` optical properties, see [`AbstractOpticalProps`](@ref)
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
