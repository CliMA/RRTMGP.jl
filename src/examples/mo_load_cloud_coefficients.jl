module mo_load_cloud_coefficients

using ..mo_optical_props
using ..mo_cloud_optics
using ..mo_simple_netcdf

export load_cld_lutcoeff, load_cld_padecoeff

"""
    load_cld_lutcoeff(cloud_spec, ds_cld_coeff)

read cloud optical property LUT coefficients from NetCDF file
"""
function load_cld_lutcoeff(::Type{FT}, ds_cld_coeff, icergh) where FT
  # class(ty_cloud_optics),     intent(inout) :: cloud_spec
  # Open cloud optical property coefficient file
  # Read LUT coefficient dimensions
  band_lims_wvn = Array{FT}(read_field(ds_cld_coeff, "bnd_limits_wavenumber"))

  # Read LUT constants
  liq = LookUpTable(read_field(ds_cld_coeff, "radliq_lwr"),
                    read_field(ds_cld_coeff, "radliq_upr"),
                    read_field(ds_cld_coeff, "radliq_fac"),
                    Array{FT}(read_field(ds_cld_coeff, "lut_extliq")),
                    Array{FT}(read_field(ds_cld_coeff, "lut_ssaliq")),
                    Array{FT}(read_field(ds_cld_coeff, "lut_asyliq")))

  ice = LookUpTable(read_field(ds_cld_coeff, "radice_lwr"),
                    read_field(ds_cld_coeff, "radice_upr"),
                    read_field(ds_cld_coeff, "radice_fac"),
                    Array{FT}(read_field(ds_cld_coeff, "lut_extice")),
                    Array{FT}(read_field(ds_cld_coeff, "lut_ssaice")),
                    Array{FT}(read_field(ds_cld_coeff, "lut_asyice")))

  cloud_spec = ty_cloud_optics_lut{FT,Int}(
    ty_optical_props_base("RRTMGP cloud optics LUT", band_lims_wvn),
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
  # class(ty_cloud_optics),       intent(inout) :: cloud_spec
  # # Spectral discretization
  # Read Pade coefficient dimensions
  band_lims_wvn = Array{FT}(read_field(ds_cld_coeff, "bnd_limits_wavenumber"))

  # Allocate cloud property Pade coefficient input arrays

  liq = PadeMethod(Array{FT}(read_field(ds_cld_coeff, "pade_extliq")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_ssaliq")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_asyliq")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_sizreg_extliq")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_sizreg_ssaliq")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_sizreg_asyliq")))

  ice = PadeMethod(Array{FT}(read_field(ds_cld_coeff, "pade_extice")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_ssaice")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_asyice")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_sizreg_extice")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_sizreg_ssaice")),
                   Array{FT}(read_field(ds_cld_coeff, "pade_sizreg_asyice")))

  cloud_spec = ty_cloud_optics_pade{FT,Int}(
    ty_optical_props_base("RRTMGP cloud optics Pade", band_lims_wvn),
    icergh,
    liq,
    ice)
  return cloud_spec
end

end # module