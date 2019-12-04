"""
    SourceFunctions

Encapsulate source function arrays for longwave/lw/internal sources
    and shortwave/sw/external source.
"""
module SourceFunctions

using DocStringExtensions
using ..OpticalProps
using ..FortranIntrinsics
import ..OpticalProps:get_ncol,get_nlay,get_ngpt

export SourceFuncLW
export SourceFuncSW
export get_ncol, get_nlay, get_ngpt

"""
    SourceFuncLW{FT, I} <: AbstractOpticalProps{FT, I}

Longwave sources: computed at layer center, at layer edges using
   spectral mapping in each direction separately, and at the surface

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceFuncLW{FT, I} <: AbstractOpticalProps{FT, I}
  "optical properties, see [`AbstractOpticalProps`](@ref)"
  optical_props::AbstractOpticalProps{FT,I}
  "Planck source at layer average temperature [W/m2] (ncol, nlay, ngpt)"
  lay_source::Array{FT,3}
  "Planck source at layer edge in increasing ilay direction [W/m2] (ncol, nlay+1, ngpt), includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
  lev_source_inc::Array{FT,3}
  "Planck source at layer edge in decreasing ilay direction [W/m2] (ncol, nlay+1, ngpt), includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
  lev_source_dec::Array{FT,3}
  "Surface source"
  sfc_source::Array{FT,2}
  function SourceFuncLW(ncol::I, nlay::I, optical_props::AbstractOpticalProps{FT}) where {FT<:AbstractFloat,I<:Int}
    ngpt = get_ngpt(optical_props)
    sfc_source = Array{FT}(undef, ncol,ngpt)
    lay_source = Array{FT}(undef, ncol,nlay,ngpt)
    lev_source_inc = Array{FT}(undef, ncol,nlay,ngpt)
    lev_source_dec = Array{FT}(undef, ncol,nlay,ngpt)
    opt_props = deepcopy(optical_props)
    return new{FT,I}(opt_props,
                     lay_source,
                     lev_source_inc,
                     lev_source_dec,
                     sfc_source)
  end
end

"""
    SourceFuncSW{FT, I} <: AbstractOpticalProps{FT, I}

Shortwave sources

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceFuncSW{FT, I} <: AbstractOpticalProps{FT, I}
  band2gpt::Array{FT,2}        # (begin g-point, end g-point) = band2gpt(2,band)
  gpt2band::Array{I,1}         # band = gpt2band(g-point)
  band_lims_wvn::Array{FT,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  name::String
  τ::Array{FT,3}
  toa_source
  lev_source_inc
  lev_source_dec
end

get_ncol(this::SourceFuncLW) = size(this.lay_source,1)

get_nlay(this::SourceFuncLW) = size(this.lay_source,2)

get_ngpt(this::SourceFuncLW) = size(this.lay_source,3)

end # module