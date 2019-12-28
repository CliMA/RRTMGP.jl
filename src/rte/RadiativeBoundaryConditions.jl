"""
    RadiativeBoundaryConditions

Atmospheric conditions used as inputs to RRTMGP.
"""
module RadiativeBoundaryConditions

using DocStringExtensions
using ..Utilities

export ShortwaveBCs
export LongwaveBCs

abstract type AbstractRadiativeBoundaryConditions{FT<:AbstractFloat} end

"""
    ShortwaveBCs{FT}

Shortwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ShortwaveBCs{FT} <: AbstractRadiativeBoundaryConditions{FT}
  "top of atmosphere flux"
  toa_flux::Array{FT}
  "surface albedo for specular (direct) radiation"
  sfc_alb_direct::Array{FT}
  "surface albedo for diffuse radiation"
  sfc_alb_diffuse::Array{FT}
  "incident diffuse flux at top of domain [W/m2] (ncol, ngpt)"
  inc_flux_diffuse::Union{Array{FT},Nothing}
  function ShortwaveBCs(toa_flux::Array{FT},
                        sfc_alb_direct::Array{FT},
                        sfc_alb_diffuse::Array{FT},
                        inc_flux_diffuse::Union{Array{FT},Nothing}=nothing) where {FT<:AbstractFloat}
    if inc_flux_diffuse ≠ nothing
      # @assert all(size(inc_flux_diffuse) .== (ncol, ngpt)) # TODO: Check sizes somehow
      @assert !any_vals_less_than(inc_flux_diffuse, FT(0))
    end
    # @assert all(size(toa_flux) .== (ncol, ngpt)) # TODO: Check sizes somehow
    @assert !any_vals_less_than(toa_flux, FT(0))
    # @assert all(size(sfc_alb_dir) .== (nband, ncol)) # TODO: Check sizes somehow
    @assert !any_vals_outside(sfc_alb_direct,  FT(0), FT(1))
    # @assert all(size(sfc_alb_dif) .== (nband, ncol)) # TODO: Check sizes somehow
    @assert !any_vals_outside(sfc_alb_diffuse, FT(0), FT(1))
    return new{FT}(toa_flux, sfc_alb_direct, sfc_alb_diffuse, inc_flux_diffuse)
  end
end



"""
    LongwaveBCs{FT}

Longwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LongwaveBCs{FT} <: AbstractRadiativeBoundaryConditions{FT}
  "spectrally-resolved emissivity (nbands, block_size)"
  sfc_emis::Array{FT,2}
  "incident flux at domain top [W/m2] (ncol, ngpts)"
  inc_flux::Union{Array{FT,2},Nothing}
  function LongwaveBCs(sfc_emis::Array{FT},
                       inc_flux::Union{Array{FT},Nothing}=nothing) where {FT<:AbstractFloat}
    # if inc_flux is present it has the right dimensions, is positive definite
    # @assert all(size(sfc_emis) .== (nband, ncol)) # TODO: check sizes somehow
    @assert !any_vals_outside(sfc_emis, FT(0), FT(1))
    if inc_flux ≠ nothing
      # @assert all(size(inc_flux) .== (ncol, ngpt)) # TODO: check sizes somehow
      @assert !any_vals_less_than(inc_flux, FT(0))
    end
    return new{eltype(sfc_emis)}(sfc_emis,inc_flux)
  end
end


end #module
