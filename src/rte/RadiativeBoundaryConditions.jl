"""
    RadiativeBoundaryConditions

Atmospheric conditions used as inputs to RRTMGP.
"""
module RadiativeBoundaryConditions

using DocStringExtensions

abstract type AbstractRadiativeBoundaryConditions{AbstractFloat} end

  if is_sw
    # toa_flux is thread private
    toa_flux = zeros(FT, ncol, ngpt)
    #
    sfc_alb_dir = zeros(FT, nbnd, ncol)
    sfc_alb_dif = zeros(FT, nbnd, ncol)
    mu0 = zeros(FT, ncol)
    # Ocean-ish values for no particular reason
    sfc_alb_dir .= FT(0.06)
    sfc_alb_dif .= FT(0.06)
    mu0 .= FT(.86)
  else
    # lw_sorces is thread private
    lw_sources = SourceFuncLW(ncol, nlay, k_dist.optical_props)

    t_sfc = zeros(FT, ncol)
    emis_sfc = zeros(FT, nbnd, ncol)
    # Surface temperature
    t_sfc .= atmos_state.t_lev[1, fmerge(nlay+1, 1, top_at_1)]
    emis_sfc .= FT(0.98)
  end


"""
    ShortwaveBCs{FT}

Shortwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ShortwaveBCs{FT} <: AbstractRadiativeBoundaryConditions{FT}
  "top of atmosphere flux"
  toa_flux
  "surface albedo for direct radiation"
  sfc_alb_dir
  "surface albedo for diffuse radiation"
  sfc_alb_dif
  inc_solar_irradiance
  flux
  inc_flux
  inc_flux_dif
end
fluxes

"""
    LongwaveBCs{FT}

Longwave boundary conditions

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LongwaveBCs{FT} <: AbstractRadiativeBoundaryConditions{FT}
  "surface skin temperatures [K]"
  t_sfc
  "spectrally-resolved emissivity (nbands, block_size)"
  emis_sfc
  "spectrally-resolved emissivity (nbands, block_size)"
  sfc_emis_spec
  function LongwaveBCs(t_sfc, emis_sfc, nbnd, block_size)
    FT = eltype(t_sfc)
    sfc_emis_spec = Array{FT}(undef, nbnd,block_size)
    for icol = 1:block_size
      for ibnd = 1:nbnd
        sfc_emis_spec[ibnd,icol] = emis_sfc[icol,b]
      end
    end
    return new{FT}(t_sfc, emis_sfc, sfc_emis_spec)
  end
end


end #module
