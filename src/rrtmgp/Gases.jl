"""
    Gases

Class for specifying gases present in the atmosphere
"""
module Gases

using InteractiveUtils
export AbstractGas
export chem_name, rfmip_name
export DryGas, h2o, o3, no2, co2, CCl4, CH3Br, CH3Cl, n2o, co, ch4, o2, n2, cfc22
export ccl4, cfc11, cfc12, hfc143a, hfc125, hfc23, hfc32, hfc134a, cf4, ch3br, ch3cl # not tested
export h2o_frgn, h2o_self # TODO: figure out what these gases are
export UncaughtGas

"""
    AbstractGas

Abstract gas class, for specifying present set of gases
"""
abstract type AbstractGas end

struct DryGas      <: AbstractGas end # Only dry component
struct h2o         <: AbstractGas end # Included by default
struct o3          <: AbstractGas end # Included by default
struct no2         <: AbstractGas end # Included by default
struct co2         <: AbstractGas end # Chosen by forcing_index
struct CCl4        <: AbstractGas end # Chosen by forcing_index
struct CH3Br       <: AbstractGas end # Chosen by forcing_index
struct CH3Cl       <: AbstractGas end # Chosen by forcing_index
struct n2o         <: AbstractGas end # Chosen by forcing_index
struct co          <: AbstractGas end # Chosen by forcing_index
struct ch4         <: AbstractGas end # Chosen by forcing_index
struct o2          <: AbstractGas end # Chosen by forcing_index
struct n2          <: AbstractGas end # Chosen by forcing_index
struct cfc22       <: AbstractGas end # Chosen by forcing_index
struct ccl4        <: AbstractGas end # Not tested
struct cfc11       <: AbstractGas end # Not tested
struct cfc12       <: AbstractGas end # Not tested
struct hfc143a     <: AbstractGas end # Not tested
struct hfc125      <: AbstractGas end # Not tested
struct hfc23       <: AbstractGas end # Not tested
struct hfc32       <: AbstractGas end # Not tested
struct hfc134a     <: AbstractGas end # Not tested
struct cf4         <: AbstractGas end # Not tested
struct ch3br       <: AbstractGas end # Not tested
struct ch3cl       <: AbstractGas end # Not tested
struct h2o_frgn    <: AbstractGas end # TODO: figure out what these gases are
struct h2o_self    <: AbstractGas end # TODO: figure out what these gases are
struct UncaughtGas <: AbstractGas end

chem_name(g::AbstractGas) = first(split(string(g), "("))

# TODO: This is not efficient and may need to be improved
function Base.convert(::Type{gas}, s::String) where {gas<:AbstractGas}
  for g in InteractiveUtils.subtypes(AbstractGas)
    if chem_name(g())==s
      return g()
    end
  end
  # throw(ArgumentError("convert: No chem_name exists for gas $(s)"))
  return UncaughtGas()
end

# Corresponding names in the RFMIP file:
rfmip_name(::DryGas)      = "dry_component"
rfmip_name(::h2o)         = "water_vapor"
rfmip_name(::o3)          = "ozone"
rfmip_name(::no2)         = "no2"
rfmip_name(::co2)         = "carbon_dioxide"
rfmip_name(::CCl4)        = "carbon_tetrachloride"
rfmip_name(::CH3Br)       = "methyl_bromide"
rfmip_name(::CH3Cl)       = "methyl_chloride"
rfmip_name(::n2o)         = "nitrous_oxide"
rfmip_name(::co)          = "carbon_monoxide"
rfmip_name(::ch4)         = "methane"
rfmip_name(::o2)          = "oxygen"
rfmip_name(::n2)          = "nitrogen"
rfmip_name(::cfc22)       = "hcfc22"
rfmip_name(::ccl4)        = "carbon_tetrachloride"
rfmip_name(::cfc11)       = "cfc11"
rfmip_name(::cfc12)       = "cfc12"
rfmip_name(::hfc143a)     = "hfc143a"
rfmip_name(::hfc125)      = "hfc125"
rfmip_name(::hfc23)       = "hfc23"
rfmip_name(::hfc32)       = "hfc32"
rfmip_name(::hfc134a)     = "hfc134a"
rfmip_name(::cf4)         = "cf4"
rfmip_name(::ch3br)       = "methyl_bromide"
rfmip_name(::ch3cl)       = "methyl_chloride"
rfmip_name(::h2o_frgn)    = "h2o_frgn"
rfmip_name(::h2o_self)    = "h2o_self"
rfmip_name(::UncaughtGas) = "UncaughtGas"

end # module