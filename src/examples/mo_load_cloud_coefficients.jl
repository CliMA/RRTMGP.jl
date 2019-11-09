module mo_load_cloud_coefficients

using ..mo_optical_props
using ..mo_cloud_optics
using ..mo_simple_netcdf

export load_cld_lutcoeff!, load_cld_padecoeff!

"""
    load_cld_lutcoeff(cloud_spec, ds_cld_coeff)

read cloud optical property LUT coefficients from NetCDF file
"""
function load_cld_lutcoeff!(cloud_spec::ty_cloud_optics{FT}, ds_cld_coeff) where FT
  # class(ty_cloud_optics),     intent(inout) :: cloud_spec
  # character(len=*),           intent(in   ) :: ds_cld_coeff
  # # -----------------
  # # Local variables
  # integer :: ds_cld_coeff, nband, nrghice, nsize_liq, nsize_ice

  # real(wp), dimension(:,:), allocatable                :: band_lims_wvn
  # # Lookup table interpolation constants
  # real(wp) :: radliq_lwr          # liquid particle size lower bound for interpolation
  # real(wp) :: radliq_upr          # liquid particle size upper bound for interpolation
  # real(wp) :: radliq_fac          # constant for calculating interpolation indices for liquid
  # real(wp) :: radice_lwr          # ice particle size lower bound for interpolation
  # real(wp) :: radice_upr          # ice particle size upper bound for interpolation
  # real(wp) :: radice_fac          # constant for calculating interpolation indices for ice
  # # LUT coefficients
  # real(wp), dimension(:,:),   allocatable :: lut_extliq   # extinction: liquid
  # real(wp), dimension(:,:),   allocatable :: lut_ssaliq   # single scattering albedo: liquid
  # real(wp), dimension(:,:),   allocatable :: lut_asyliq   # asymmetry parameter: liquid
  # real(wp), dimension(:,:,:), allocatable :: lut_extice   # extinction: ice
  # real(wp), dimension(:,:,:), allocatable :: lut_ssaice   # single scattering albedo: ice
  # real(wp), dimension(:,:,:), allocatable :: lut_asyice   # asymmetry parameter: ice
  # -----------------
  # Open cloud optical property coefficient file

  # Read LUT coefficient dimensions
  nband     = ds_cld_coeff.dim["nband"]
  nrghice   = ds_cld_coeff.dim["nrghice"]
  nsize_liq = ds_cld_coeff.dim["nsize_liq"]
  nsize_ice = ds_cld_coeff.dim["nsize_ice"]

  band_lims_wvn = Array{FT}(undef, 2, nband)
  band_lims_wvn .= read_field(ds_cld_coeff, "bnd_limits_wavenumber", 2, nband)

  # Read LUT constants
  radliq_lwr = read_field(ds_cld_coeff, "radliq_lwr")
  radliq_upr = read_field(ds_cld_coeff, "radliq_upr")
  radliq_fac = read_field(ds_cld_coeff, "radliq_fac")
  radice_lwr = read_field(ds_cld_coeff, "radice_lwr")
  radice_upr = read_field(ds_cld_coeff, "radice_upr")
  radice_fac = read_field(ds_cld_coeff, "radice_fac")

  # Allocate cloud property lookup table input arrays
  lut_extliq = Array{FT}(undef, nsize_liq, nband)
  lut_ssaliq = Array{FT}(undef, nsize_liq, nband)
  lut_asyliq = Array{FT}(undef, nsize_liq, nband)
  lut_extice = Array{FT}(undef, nsize_ice, nband, nrghice)
  lut_ssaice = Array{FT}(undef, nsize_ice, nband, nrghice)
  lut_asyice = Array{FT}(undef, nsize_ice, nband, nrghice)

  # Read LUT coefficients
  lut_extliq .= read_field(ds_cld_coeff, "lut_extliq",  nsize_liq, nband)
  lut_ssaliq .= read_field(ds_cld_coeff, "lut_ssaliq",  nsize_liq, nband)
  lut_asyliq .= read_field(ds_cld_coeff, "lut_asyliq",  nsize_liq, nband)
  lut_extice .= read_field(ds_cld_coeff, "lut_extice",  nsize_ice, nband, nrghice)
  lut_ssaice .= read_field(ds_cld_coeff, "lut_ssaice",  nsize_ice, nband, nrghice)
  lut_asyice .= read_field(ds_cld_coeff, "lut_asyice",  nsize_ice, nband, nrghice)

  load_lut!(cloud_spec, band_lims_wvn,
        radliq_lwr, radliq_upr, radliq_fac,
        radice_lwr, radice_upr, radice_fac,
        lut_extliq, lut_ssaliq, lut_asyliq,
        lut_extice, lut_ssaice, lut_asyice)
end

"""
    load_cld_padecoeff

read cloud optical property Pade coefficients from NetCDF file
"""
function load_cld_padecoeff!(cloud_spec::ty_cloud_optics{FT}, ds_cld_coeff) where FT
  # class(ty_cloud_optics),       intent(inout) :: cloud_spec
  # character(len=*),             intent(in   ) :: ds_cld_coeff
  # # -----------------
  # # Local variables
  # integer :: ds_cld_coeff, nband, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g, nbound

  # # Spectral discretization
  # real(wp), dimension(:,:), allocatable :: band_lims_wvn

  # # Pade coefficients
  # real(wp), dimension(:,:,:),   allocatable :: pade_extliq   # extinction: liquid
  # real(wp), dimension(:,:,:),   allocatable :: pade_ssaliq   # single scattering albedo: liquid
  # real(wp), dimension(:,:,:),   allocatable :: pade_asyliq   # asymmetry parameter: liquid
  # real(wp), dimension(:,:,:,:), allocatable :: pade_extice   # extinction: ice
  # real(wp), dimension(:,:,:,:), allocatable :: pade_ssaice   # single scattering albedo: ice
  # real(wp), dimension(:,:,:,:), allocatable :: pade_asyice   # asymmetry parameter: ice

  # # Pade particle size regime boundaries
  # real(wp),  dimension(:),       allocatable :: pade_sizreg_extliq
  # real(wp),  dimension(:),       allocatable :: pade_sizreg_ssaliq
  # real(wp),  dimension(:),       allocatable :: pade_sizreg_asyliq
  # real(wp),  dimension(:),       allocatable :: pade_sizreg_extice
  # real(wp),  dimension(:),       allocatable :: pade_sizreg_ssaice
  # real(wp),  dimension(:),       allocatable :: pade_sizreg_asyice
  # # -----------------
  # Open cloud optical property coefficient file

  # Read Pade coefficient dimensions
  nband        = ds_cld_coeff.dim["nband"]
  nrghice      = ds_cld_coeff.dim["nrghice"]
  nsizereg     = ds_cld_coeff.dim["nsizereg"]
  ncoeff_ext   = ds_cld_coeff.dim["ncoeff_ext"]
  ncoeff_ssa_g = ds_cld_coeff.dim["ncoeff_ssa_g"]
  nbound       = ds_cld_coeff.dim["nbound"]

  #
  band_lims_wvn = Array{FT}(undef, 2, nband)
  band_lims_wvn .= read_field(ds_cld_coeff, "bnd_limits_wavenumber", 2, nband)

  # Allocate cloud property Pade coefficient input arrays
  pade_extliq = Array{FT}(undef, nband, nsizereg, ncoeff_ext)
  pade_ssaliq = Array{FT}(undef, nband, nsizereg, ncoeff_ssa_g)
  pade_asyliq = Array{FT}(undef, nband, nsizereg, ncoeff_ssa_g)
  pade_extice = Array{FT}(undef, nband, nsizereg, ncoeff_ext,   nrghice)
  pade_ssaice = Array{FT}(undef, nband, nsizereg, ncoeff_ssa_g, nrghice)
  pade_asyice = Array{FT}(undef, nband, nsizereg, ncoeff_ssa_g, nrghice)

  pade_extliq .= read_field(ds_cld_coeff, "pade_extliq", nband, nsizereg, ncoeff_ext)
  pade_ssaliq .= read_field(ds_cld_coeff, "pade_ssaliq", nband, nsizereg, ncoeff_ssa_g)
  pade_asyliq .= read_field(ds_cld_coeff, "pade_asyliq", nband, nsizereg, ncoeff_ssa_g)
  pade_extice .= read_field(ds_cld_coeff, "pade_extice", nband, nsizereg, ncoeff_ext, nrghice)
  pade_ssaice .= read_field(ds_cld_coeff, "pade_ssaice", nband, nsizereg, ncoeff_ssa_g, nrghice)
  pade_asyice .= read_field(ds_cld_coeff, "pade_asyice", nband, nsizereg, ncoeff_ssa_g, nrghice)

  # Allocate cloud property Pade coefficient particle size boundary input arrays
  pade_sizreg_extliq = Array{FT}(undef, nbound)
  pade_sizreg_ssaliq = Array{FT}(undef, nbound)
  pade_sizreg_asyliq = Array{FT}(undef, nbound)
  pade_sizreg_extice = Array{FT}(undef, nbound)
  pade_sizreg_ssaice = Array{FT}(undef, nbound)
  pade_sizreg_asyice = Array{FT}(undef, nbound)

  pade_sizreg_extliq .= read_field(ds_cld_coeff, "pade_sizreg_extliq", nbound)
  pade_sizreg_ssaliq .= read_field(ds_cld_coeff, "pade_sizreg_ssaliq", nbound)
  pade_sizreg_asyliq .= read_field(ds_cld_coeff, "pade_sizreg_asyliq", nbound)
  pade_sizreg_extice .= read_field(ds_cld_coeff, "pade_sizreg_extice", nbound)
  pade_sizreg_ssaice .= read_field(ds_cld_coeff, "pade_sizreg_ssaice", nbound)
  pade_sizreg_asyice .= read_field(ds_cld_coeff, "pade_sizreg_asyice", nbound)

  load_Pade!(cloud_spec, band_lims_wvn,
        pade_extliq, pade_ssaliq, pade_asyliq,
        pade_extice, pade_ssaice, pade_asyice,
        pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq,
        pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice)
end

end # module