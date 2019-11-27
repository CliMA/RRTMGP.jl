module LoadCoefficients

using ..GasConcentrations
using ..GasOptics
using ..FortranIntrinsics
using NCDatasets

export load_and_init

function read_char_vec(ds, var_name)
  tmp = ds[var_name][:]
  return Array{String}([trim(join(tmp[:,i])) for i = 1:size(tmp, 2)])
end

"""
    load_and_init(ds, available_gases::GasConcs{FT}) where FT

Initialize the gas optics class with data. The calls look slightly different depending
  on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
"""
function load_and_init(ds, ::Type{FT}, available_gases::Vector{String}) where {FT<:AbstractFloat}

  # Reading the properties from the NetCDF file

  I = Int
  kmajor                          = Array{FT}(ds["kmajor"][:])
  key_species                     = Array{I}(ds["key_species"][:])
  band2gpt                        = Array{I}(ds["bnd_limits_gpt"][:])
  band_lims                       = Array{FT}(ds["bnd_limits_wavenumber"][:])
  gas_names                       = read_char_vec(ds, "gas_names")
  gas_minor                       = read_char_vec(ds, "gas_minor")
  identifier_minor                = read_char_vec(ds, "identifier_minor")

  rayl_lower = haskey(ds,"rayl_lower") ? Array{FT}(ds["rayl_lower"][:]) : nothing
  rayl_upper = haskey(ds,"rayl_upper") ? Array{FT}(ds["rayl_upper"][:]) : nothing
  @assert haskey(ds,"rayl_lower") == haskey(ds,"rayl_upper")

  lower = GasOpticsVars{FT,I}(Array{I}(ds["minor_limits_gpt_lower"][:]),
                              Array{Bool}(ds["minor_scales_with_density_lower"][:]),
                              Array{Bool}(ds["scale_by_complement_lower"][:]),
                              Array{I}(ds["kminor_start_lower"][:]),
                              Array{FT}(ds["kminor_lower"][:]),
                              read_char_vec(ds, "scaling_gas_lower"),
                              read_char_vec(ds, "minor_gases_lower"))
  upper = GasOpticsVars{FT,I}(Array{I}(ds["minor_limits_gpt_upper"][:]),
                              Array{Bool}(ds["minor_scales_with_density_upper"][:]),
                              Array{Bool}(ds["scale_by_complement_upper"][:]),
                              Array{I}(ds["kminor_start_upper"][:]),
                              Array{FT}(ds["kminor_upper"][:]),
                              read_char_vec(ds, "scaling_gas_upper"),
                              read_char_vec(ds, "minor_gases_upper"))

  ref = Reference(Array{FT}(ds["press_ref"][:]),
                  Array{FT}(ds["temp_ref"][:]),
                  FT(ds["press_ref_trop"][:]),
                  Array{FT}(ds["vmr_ref"][:]),
                  available_gases,
                  gas_names
                  )

  args = (available_gases,
          gas_names,
          key_species,
          band2gpt,
          band_lims,
          kmajor,
          gas_minor,
          identifier_minor,
          lower,
          upper)

  if haskey(ds,"totplnk")
    totplnk     = Array{FT}(ds["totplnk"][:])
    planck_frac = Array{FT}(ds["plank_fraction"][:])
    kdist = load_totplnk(totplnk, planck_frac, rayl_lower, rayl_upper, ref, args...)
  else
    solar_src = ds["solar_source"][:]
    kdist = load_solar_source(solar_src, rayl_lower, rayl_upper, ref, args...)
  end
  return kdist

end


end #module