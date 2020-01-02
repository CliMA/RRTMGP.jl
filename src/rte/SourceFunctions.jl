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

export SourceFuncLongWave, SourceFuncLongWavePGP
export SourceFuncShortWave
export get_ncol, get_nlay, get_ngpt

abstract type AbstractSourceFunc{FT,I} <: AbstractOpticalProps{FT, I} end

"""
    SourceFuncLongWave{FT, I} <: AbstractOpticalProps{FT, I}

Longwave sources: computed at layer center, at layer edges using
spectral mapping in each direction separately, and at the surface

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceFuncLongWave{FT, I} <: AbstractSourceFunc{FT, I}
  "optical properties, see [`OpticalPropsBase`](@ref)"
  optical_props::OpticalPropsBase{FT,I}
  "Planck source at layer average temperature [W/m2] (ncol, nlay, ngpt)"
  lay_source::Array{FT,3}
  "Planck source at layer edge in increasing ilay direction [W/m2] (ncol, nlay+1, ngpt), includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
  lev_source_inc::Array{FT,3}
  "Planck source at layer edge in decreasing ilay direction [W/m2] (ncol, nlay+1, ngpt), includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
  lev_source_dec::Array{FT,3}
  "Surface source"
  sfc_source::Array{FT,2}
  "Temporary array for Planck fraction"
  p_frac::Array{FT,3}
end

function SourceFuncLongWave(ncol::I, nlay::I, optical_props::OpticalPropsBase{FT}) where {FT<:AbstractFloat,I<:Int}
  ngpt = get_ngpt(optical_props)
  op = deepcopy(optical_props)
  p_frac         = zeros(FT, ncol,nlay,ngpt)
  lay_source     = zeros(FT, ncol,nlay,ngpt)
  lev_source_inc = zeros(FT, ncol,nlay,ngpt)
  lev_source_dec = zeros(FT, ncol,nlay,ngpt)
  sfc_source     = zeros(FT, ncol,ngpt)
  return SourceFuncLongWave{FT,I}(op,
                                  lay_source,
                                  lev_source_inc,
                                  lev_source_dec,
                                  sfc_source,
                                  p_frac)
end

"""
    SourceFuncLongWavePGP{FT, I} <: OpticalPropsBase{FT, I}

Longwave sources: computed at layer center, at layer edges using
spectral mapping in each direction separately, and at the surface
per grid point

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceFuncLongWavePGP{FT, I} <: AbstractSourceFunc{FT, I}
  "optical properties, see [`OpticalPropsBase`](@ref)"
  optical_props::OpticalPropsBase{FT,I}
  "Planck source at layer average temperature [W/m2] (ngpt)"
  lay_source::Array{FT,1}
  "Planck source at layer edge in increasing ilay direction [W/m2] (ngpt), includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
  lev_source_inc::Array{FT,1}
  "Planck source at layer edge in decreasing ilay direction [W/m2] (ngpt), includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
  lev_source_dec::Array{FT,1}
  "Surface source"
  sfc_source::Array{FT,1}
  "Temporary array for Planck fraction"
  p_frac::Array{FT,1}
end

function Base.convert(::Type{SourceFuncLongWave}, data::Array{SourceFuncLongWavePGP{FT,I}}, sfc_lay::I) where {FT,I}
  s = size(data)
  op = first(data).optical_props
  ngpt = get_ngpt(op)
  return SourceFuncLongWave{FT,I}(
    op,
    Array{FT}([data[i,j].lay_source[k]       for i in 1:s[1], j in 1:s[2], k in 1:ngpt]),
    Array{FT}([data[i,j].lev_source_inc[k]   for i in 1:s[1], j in 1:s[2], k in 1:ngpt]),
    Array{FT}([data[i,j].lev_source_dec[k]   for i in 1:s[1], j in 1:s[2], k in 1:ngpt]),
    Array{FT}([data[i,sfc_lay].sfc_source[k] for i in 1:s[1], k in 1:ngpt]),
    Array{FT}([data[i,j].p_frac[k]           for i in 1:s[1], j in 1:s[2], k in 1:ngpt])
    )
end

Base.convert(::Type{Array{SourceFuncLongWavePGP}}, data::SourceFuncLongWave{FT,I}) where {FT,I} =
  [SourceFuncLongWavePGP{FT,I}(
    data.optical_props,
    data.lay_source[i,j,:],
    data.lev_source_inc[i,j,:],
    data.lev_source_dec[i,j,:],
    data.sfc_source[i,:],
    data.p_frac[i,j,:]
    ) for i in 1:size(data.lay_source,1),
          j in 1:size(data.lay_source,2)]

"""
    SourceFuncShortWave{FT, I} <: AbstractOpticalProps{FT, I}

Shortwave sources

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceFuncShortWave{FT, I} <: AbstractSourceFunc{FT, I}
  "optical properties, see [`OpticalPropsBase`](@ref)"
  optical_props::OpticalPropsBase{FT,I}
  "top of atmosphere source"
  toa_source
  "Source at layer edge in increasing ilay direction"
  lev_source_inc
  "Source at layer edge in decreasing ilay direction"
  lev_source_dec
end

get_ncol(this::SourceFuncLongWave) = size(this.lay_source,1)
get_nlay(this::SourceFuncLongWave) = size(this.lay_source,2)
get_ngpt(this::SourceFuncLongWave) = size(this.lay_source,3)
get_ngpt(this::SourceFuncLongWavePGP) = size(this.lay_source,1)

end # module