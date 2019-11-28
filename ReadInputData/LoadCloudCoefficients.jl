using RRTMGP.OpticalProps
using RRTMGP.CloudOptics

"""
    load_cld_lutcoeff(cloud_spec, ds_cld_coeff)

read cloud optical property LUT coefficients from NetCDF file
"""
function load_cld_lutcoeff(::Type{FT}, ds_cld_coeff, icergh) where FT
  # Open cloud optical property coefficient file
  # Read LUT coefficient dimensions
  band_lims_wvn = Array{FT}(ds_cld_coeff["bnd_limits_wavenumber"][:])

  # Read LUT constants
  liq = LookUpTable(ds_cld_coeff["radliq_lwr"][:],
                    ds_cld_coeff["radliq_upr"][:],
                    ds_cld_coeff["radliq_fac"][:],
                    Array{FT}(ds_cld_coeff["lut_extliq"][:]),
                    Array{FT}(ds_cld_coeff["lut_ssaliq"][:]),
                    Array{FT}(ds_cld_coeff["lut_asyliq"][:]))

  ice = LookUpTable(ds_cld_coeff["radice_lwr"][:],
                    ds_cld_coeff["radice_upr"][:],
                    ds_cld_coeff["radice_fac"][:],
                    Array{FT}(ds_cld_coeff["lut_extice"][:]),
                    Array{FT}(ds_cld_coeff["lut_ssaice"][:]),
                    Array{FT}(ds_cld_coeff["lut_asyice"][:]))

  cloud_spec = CloudOpticsLUT{FT,Int}(
    OpticalPropsBase("RRTMGP cloud optics LUT", band_lims_wvn),
    icergh,
    liq,
    ice)
  return cloud_spec
end

"""
    load_cld_padecoeff

read cloud optical property Pade coefficients from NetCDF file
"""
function load_cld_padecoeff(::Type{FT}, ds_cld_coeff, icergh) where FT
  # # Spectral discretization
  # Read Pade coefficient dimensions
  band_lims_wvn = Array{FT}(ds_cld_coeff["bnd_limits_wavenumber"][:])

  # Allocate cloud property Pade coefficient input arrays

  liq = PadeMethod(Array{FT}(ds_cld_coeff["pade_extliq"][:]),
                   Array{FT}(ds_cld_coeff["pade_ssaliq"][:]),
                   Array{FT}(ds_cld_coeff["pade_asyliq"][:]),
                   Array{FT}(ds_cld_coeff["pade_sizreg_extliq"][:]),
                   Array{FT}(ds_cld_coeff["pade_sizreg_ssaliq"][:]),
                   Array{FT}(ds_cld_coeff["pade_sizreg_asyliq"][:]))

  ice = PadeMethod(Array{FT}(ds_cld_coeff["pade_extice"][:]),
                   Array{FT}(ds_cld_coeff["pade_ssaice"][:]),
                   Array{FT}(ds_cld_coeff["pade_asyice"][:]),
                   Array{FT}(ds_cld_coeff["pade_sizreg_extice"][:]),
                   Array{FT}(ds_cld_coeff["pade_sizreg_ssaice"][:]),
                   Array{FT}(ds_cld_coeff["pade_sizreg_asyice"][:]))

  cloud_spec = CloudOpticsPade{FT,Int}(
    OpticalPropsBase("RRTMGP cloud optics Pade", band_lims_wvn),
    icergh,
    liq,
    ice)
  return cloud_spec
end
