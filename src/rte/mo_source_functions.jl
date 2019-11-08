"""
    mo_source_functions

Encapsulate source function arrays for longwave/lw/internal sources
    and shortwave/sw/external source.
"""
module mo_source_functions

using DocStringExtensions
using ..mo_optical_props
using ..fortran_intrinsics
import ..mo_optical_props:get_ncol,get_nlay,get_ngpt

export ty_source_func_lw
export get_ncol, get_nlay, get_ngpt
export ty_source_func_sw

"""
    ty_source_func_lw{FT, I} <: ty_optical_props{FT, I}

Longwave sources: computed at layer center, at layer edges using
   spectral mapping in each direction separately, and at the surface

# Fields

$(DocStringExtensions.FIELDS)
"""
struct ty_source_func_lw{FT, I} <: ty_optical_props{FT, I}
  optical_props::ty_optical_props{FT,I}
  lay_source     # Planck source at layer average temperature [W/m2] (ncol, nlay, ngpt)
  lev_source_inc # Planck source at layer edge in increasing ilay direction [W/m2] (ncol, nlay+1, ngpt)
  lev_source_dec # Planck source at layer edge in decreasing ilay direction [W/m2] (ncol, nlay+1, ngpt)
  sfc_source
  function ty_source_func_lw(ncol::I, nlay::I, optical_props::ty_optical_props{FT}) where {FT,I}
    ngpt = get_ngpt(optical_props)
    sfc_source = Array{FT}(undef, ncol,ngpt)
    lay_source = Array{FT}(undef, ncol,nlay,ngpt)
    lev_source_inc = Array{FT}(undef, ncol,nlay,ngpt)
    lev_source_dec = Array{FT}(undef, ncol,nlay,ngpt)
    opt_props = deepcopy(optical_props)
    new{FT,I}(opt_props,
              lay_source,
              lev_source_inc,
              lev_source_dec,
              sfc_source)
  end
end

"""
    ty_source_func_sw{FT, I} <: ty_optical_props{FT, I}

Shortwave sources

# Fields

$(DocStringExtensions.FIELDS)
"""
struct ty_source_func_sw{FT, I} <: ty_optical_props{FT, I}
  band2gpt::Array{FT,2}        # (begin g-point, end g-point) = band2gpt(2,band)
  gpt2band::Array{I,1}         # band = gpt2band(g-point)
  band_lims_wvn::Array{FT,2}   # (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  name::String
  tau::Array{FT,3}
  toa_source
  lev_source_inc
  lev_source_dec
end

is_allocated(this::ty_source_func_lw) = is_initialized(this.optical_props) && allocated(this.sfc_source)

get_ncol(this::ty_source_func_lw) = is_allocated(this) ? size(this.lay_source,1) : 0

get_nlay(this::ty_source_func_lw) = is_allocated(this) ? size(this.lay_source,2) : 0

get_ngpt(this::ty_source_func_lw) = is_allocated(this) ? size(this.lay_source,3) : 0

end # module