module mo_load_coefficients

  using ..mo_gas_concentrations
  using ..mo_gas_optics_rrtmgp
  using ..fortran_intrinsics
  using NCDatasets

  export load_and_init

  function read_char_vec(ds, var_name)
    tmp = ds[var_name][:]
    return Array{String}([trim(join(tmp[:,i])) for i = 1:size(tmp, 2)])
  end

#--------------------------
  function load_and_init(ds, available_gases::ty_gas_concs{FT}) where FT

  # Reading the properties from the NetCDF file

  band2gpt                        = Array{Int}(ds["bnd_limits_gpt"][:])
  minor_limits_gpt_lower          = Array{Int}(ds["minor_limits_gpt_lower"][:])
  minor_limits_gpt_upper          = Array{Int}(ds["minor_limits_gpt_upper"][:])
  key_species                     = Array{Int}(ds["key_species"][:])
  kminor_start_lower              = Array{Int}(ds["kminor_start_lower"][:])
  kminor_start_upper              = Array{Int}(ds["kminor_start_upper"][:])
  minor_scales_with_density_lower = Array{Bool}(ds["minor_scales_with_density_lower"][:])
  minor_scales_with_density_upper = Array{Bool}(ds["minor_scales_with_density_upper"][:])
  scale_by_complement_lower       = Array{Bool}(ds["scale_by_complement_lower"][:])
  scale_by_complement_upper       = Array{Bool}(ds["scale_by_complement_upper"][:])
  band_lims                       = Array{FT}(ds["bnd_limits_wavenumber"][:])
  press_ref                       = Array{FT}(ds["press_ref"][:])
  temp_ref                        = Array{FT}(ds["temp_ref"][:])
  kminor_lower                    = Array{FT}(ds["kminor_lower"][:])
  kminor_upper                    = Array{FT}(ds["kminor_upper"][:])
  vmr_ref                         = Array{FT}(ds["vmr_ref"][:])
  kmajor                          = Array{FT}(ds["kmajor"][:])
  temp_ref_p                      = FT(ds["absorption_coefficient_ref_P"][:])
  temp_ref_t                      = FT(ds["absorption_coefficient_ref_T"][:])
  press_ref_trop                  = FT(ds["press_ref_trop"][:])

  gas_names                       = read_char_vec(ds, "gas_names")
  gas_minor                       = read_char_vec(ds, "gas_minor")
  identifier_minor                = read_char_vec(ds, "identifier_minor")
  minor_gases_lower               = read_char_vec(ds, "minor_gases_lower")
  minor_gases_upper               = read_char_vec(ds, "minor_gases_upper")
  scaling_gas_lower               = read_char_vec(ds, "scaling_gas_lower")
  scaling_gas_upper               = read_char_vec(ds, "scaling_gas_upper")

  rayl_lower = haskey(ds,"rayl_lower") ? Array{FT}(ds["rayl_lower"][:]) : nothing
  rayl_upper = haskey(ds,"rayl_upper") ? Array{FT}(ds["rayl_upper"][:]) : nothing

  args = (available_gases,
          gas_names,
          key_species,
          band2gpt,
          band_lims,
          press_ref,
          press_ref_trop,
          temp_ref,
          temp_ref_p,
          temp_ref_t,
          vmr_ref,
          kmajor,
          kminor_lower,
          kminor_upper,
          gas_minor,
          identifier_minor,
          minor_gases_lower,
          minor_gases_upper,
          minor_limits_gpt_lower,
          minor_limits_gpt_upper,
          minor_scales_with_density_lower,
          minor_scales_with_density_upper,
          scaling_gas_lower,
          scaling_gas_upper,
          scale_by_complement_lower,
          scale_by_complement_upper,
          kminor_start_lower,
          kminor_start_upper,)

 #   ! Initialize the gas optics class with data. The calls look slightly different depending
 #   !   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
 #   ! gas_optics%load() returns a string; a non-empty string indicates an error.

  if haskey(ds,"totplnk")
    totplnk     =      ds["totplnk"][:]
    planck_frac =      ds["plank_fraction"][:]
    kdist = load_totplnk(totplnk, planck_frac, rayl_lower, rayl_upper, args...)
  else
    solar_src = ds["solar_source"][:]
    kdist = load_solar_source(solar_src, rayl_lower, rayl_upper, args...)
  end
  return kdist

  #-------------------------------

  end
#--------------------------


end
