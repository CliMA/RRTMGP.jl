"""
    AtmosphericStates

Atmospheric conditions used as inputs to RRTMGP.
"""
module AtmosphericStates

using DocStringExtensions

using ..Utilities
using ..FortranIntrinsics
using ..PhysicalConstants
using ..Gases
using ..GasConcentrations
using ..ReferenceStates

export AtmosphericState, AtmosphericStatePGP

abstract type AbstractAtmosphericState{AbstractFloat,Integer} end

include("Interpolation.jl")

"""
    SimpleGrid{FT}

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SimpleGrid{FT}
  z
end

"""
    get_col_dry(vmr_h2o, plev, tlay, latitude=nothing)

Utility function, provided for user convenience
computes column amounts of dry air using hydrostatic equation

# input
 - `vmr_h2o` volume mixing ratio of water vapor to dry air; (ncol,nlay)
 - `plev` Layer boundary pressures [Pa] (ncol,nlay+1)
 - `tlay` Layer temperatures [K] (ncol,nlay)
 - `latitude` Latitude [degrees] (ncol)
"""
function get_col_dry(vmr_h2o, plev, tlay, latitude=nothing)
  FT = eltype(plev)

  # first and second term of Helmert formula
  helmert1 = FT(9.80665)
  helmert2 = FT(0.02586)
  # local variables
  g0         = zeros(FT, size(tlay,1)             ) # (ncol)
  delta_plev = zeros(FT, size(tlay,1),size(tlay,2)) # (ncol,nlay)
  m_air      = zeros(FT, size(tlay,1),size(tlay,2)) # average mass of air; (ncol,nlay)

  nlay = size(tlay, 2)
  nlev = size(plev, 2)

  if latitude ≠ nothing
    g0 .= helmert1 - helmert2 * cos(FT(2) * π * latitude / FT(180)) # acceleration due to gravity [m/s^2]
  else
    g0 .= grav(FT)
  end
  delta_plev .= abs.(plev[:,1:nlev-1] .- plev[:,2:nlev])

  # Get average mass of moist air per mole of moist air
  m_air .= (m_dry(FT) .+ m_h2o(FT) .* vmr_h2o) ./ (1 .+ vmr_h2o)

  # Hydrostatic equation
  col_dry = zeros(FT, size(tlay,1),size(tlay,2))
  col_dry .= FT(10) .* delta_plev .* avogad(FT) ./ (FT(1000)*m_air .* FT(100) .* spread(g0, 2, nlay))
  col_dry .= col_dry ./ (FT(1) .+ vmr_h2o)
  return col_dry
end

"""
    AtmosphericState{FT}

as = AtmosphericState(gas_conc, p_lay, p_lev, t_lay, t_lev)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AtmosphericState{FT,I} <: AbstractAtmosphericState{FT,I}
  "Gas concentrations, in the form of volume mixing ratios"
  gas_conc::GasConcs{FT,I}
  "layer pressures [Pa, mb]; (ncol,nlay)"
  p_lay::Array{FT,2}
  "level pressures [Pa, mb]; (ncol,nlay+1)"
  p_lev::Array{FT,2}
  "layer temperatures [K]; (ncol,nlay)"
  t_lay::Array{FT,2}
  "level temperatures [K]; (ncol,nlay+1)"
  t_lev::Array{FT,2}
  "surface temperatures [K]; (ncol)"
  t_sfc::Vector{FT}
  "column amounts for each gas, plus col_dry. gas amounts [molec/cm^2]"
  col_gas::Array{FT,3}
  "Number of molecules per cm-2 of dry air"
  col_dry::Array{FT,2}
  "Indicates whether arrays are ordered in the vertical with 1 at the top or the bottom of the domain."
  top_at_1::Bool
  "Number of layers."
  nlay::I
  "Number of columns."
  ncol::I
end
function AtmosphericState(gas_conc::GasConcs{FT,I},
                          p_lay::Array{FT,2},
                          p_lev::Array{FT,2},
                          t_lay::Array{FT,2},
                          t_lev::Union{Array{FT,2},Nothing},
                          ref::ReferenceState,
                          col_dry::Union{Array{FT,2},Nothing}=nothing,
                          t_sfc::Union{Vector{FT},Nothing}=nothing) where {I<:Int,FT<:AbstractFloat}
  nlay = size(p_lay, 2)
  ncol = size(p_lay, 1)
  if t_lev==nothing
    t_lev = interpolate_var(p_lay, p_lev, t_lay, ncol, nlay)
  end

  check_extent(p_lay, (ncol, nlay  ), "p_lay")
  check_extent(p_lev, (ncol, nlay+1), "p_lev")
  check_extent(t_lay, (ncol, nlay  ), "t_lay")
  check_extent(t_lev, (ncol, nlay+1), "t_lev")
  check_range(p_lay, ref.press_min, ref.press_max, "p_lay")
  check_range(p_lev, ref.press_min, ref.press_max, "p_lev")
  check_range(t_lay, ref.temp_min,  ref.temp_max,  "t_lay")
  check_range(t_lev, ref.temp_min,  ref.temp_max,  "t_lev")

  ngas    = length(gas_conc.gas_names)
  vmr     = Array{FT}(undef, ncol,nlay,ngas) # volume mixing ratios
  col_gas = Array{FT}(undef, ncol,nlay,ngas+1) # column amounts for each gas, plus col_dry

  # Fill out the array of volume mixing ratios
  for gas in gas_conc.gas_names
    i_gas = loc_in_array(gas, gas_conc.gas_names)
    vmr[:,:,i_gas] .= gas_conc.concs[i_gas,:,:]
  end

  # Compute dry air column amounts (number of molecule per cm^2) if user hasn't provided them
  idx_h2o = loc_in_array(h2o(), gas_conc.gas_names)
  if col_dry == nothing
    col_dry = get_col_dry(vmr[:,:,idx_h2o], p_lev, t_lay) # dry air column amounts computation
  end

  check_extent(col_dry, (ncol, nlay), "col_dry")
  check_range(col_dry, FT(0), floatmax(FT), "col_dry")

  # compute column gas amounts [molec/cm^2]
  col_gas[:,:,1] .= col_dry
  for igas = 1:ngas
    col_gas[:,:,igas+1] .= vmr[:,:,igas] .* col_dry
  end
  # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
  top_at_1 = p_lay[1, 1] < p_lay[1, nlay]

  if t_sfc == nothing
    t_sfc = Array{FT}(undef, ncol)
    t_sfc .= t_lev[1, fmerge(nlay+1, 1, top_at_1)]
  end
  check_extent(t_sfc, ncol, "t_sfc")
  check_range(t_sfc, ref.temp_min,  ref.temp_max,  "t_sfc")

  return AtmosphericState{FT,I}(gas_conc,
                                p_lay,
                                p_lev,
                                t_lay,
                                t_lev,
                                t_sfc,
                                col_gas,
                                col_dry,
                                top_at_1,
                                nlay,
                                ncol)
end

"""
    AtmosphericStatePGP{FT}

as = AtmosphericStatePGP(gas_conc, p_lay, p_lev, t_lay, t_lev)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct AtmosphericStatePGP{FT,I} <: AbstractAtmosphericState{FT,I}
  "Gas concentrations, in the form of volume mixing ratios"
  gas_conc::GasConcsPGP{FT}
  "layer pressures [Pa, mb]"
  p_lay::FT
  "level pressures [Pa, mb]"
  p_lev::Vector{FT}
  "layer temperatures [K]"
  t_lay::FT
  "level temperatures [K]"
  t_lev::Vector{FT}
  "surface temperatures [K]"
  t_sfc::FT
  "column amounts for each gas, plus col_dry. gas amounts [molec/cm^2]"
  col_gas::Array{FT,1}
  "Number of molecules per cm-2 of dry air"
  col_dry::FT
  "Indicates whether arrays are ordered in the vertical with 1 at the top or the bottom of the domain."
  top_at_1::Bool
end

function AtmosphericStatePGP(gas_conc::GasConcsPGP{FT},
                             p_lay::Array{FT},
                             p_lev::FT,
                             t_lay::Array{FT},
                             t_lev::Union{FT,Nothing},
                             ref::ReferenceState,
                             col_dry::FT,
                             t_sfc::FT,
                             top_at_1::Bool) where {FT<:AbstractFloat}
  if t_lev==nothing
    t_lev = interpolate_var(p_lay, p_lev, t_lay)
  end

  check_range(p_lay, ref.press_min, ref.press_max, "p_lay")
  check_range(p_lev, ref.press_min, ref.press_max, "p_lev")
  check_range(t_lay, ref.temp_min,  ref.temp_max,  "t_lay")
  check_range(t_lev, ref.temp_min,  ref.temp_max,  "t_lev")

  ngas    = length(gas_conc.gas_names)
  I = eltype(ngas)
  vmr     = Array{FT}(undef, ngas) # volume mixing ratios
  col_gas = Array{FT}(undef, ngas+1) # column amounts for each gas, plus col_dry

  # Fill out the array of volume mixing ratios
  for gas in gas_conc.gas_names
    i_gas = loc_in_array(gas, gas_conc.gas_names)
    vmr[i_gas] .= gas_conc.concs[i_gas]
  end

  check_range(col_dry, FT(0), floatmax(FT), "col_dry")

  # compute column gas amounts [molec/cm^2]
  col_gas[1] .= col_dry
  for igas = 1:ngas
    col_gas[igas+1] .= vmr[igas] .* col_dry
  end
  # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?

  check_range(t_sfc, ref.temp_min,  ref.temp_max,  "t_sfc")

  return AtmosphericStatePGP{FT,I}(gas_conc,
                                   p_lay,
                                   p_lev,
                                   t_lay,
                                   t_lev,
                                   t_sfc,
                                   col_gas,
                                   col_dry,
                                   top_at_1)
end

function Base.convert(::Type{AtmosphericState}, data::Array{AtmosphericStatePGP{FT,I}}) where {FT,I}
  s = size(data)
  ncol,nlay = s
  gas_conc = convert(GasConcs, [data[i,j].gas_conc for i in 1:s[1], j in 1:s[2]])
  top_at_1 = first(data).top_at_1
  ngas    = length(gas_conc.gas_names)

  col_gas = Array{FT}(undef, ncol,nlay,ngas+1)
  for i in 1:s[1]
    for j in 1:s[2]
      for igas in 1:ngas+1
        col_gas[i,j,igas] = data[i,j].col_gas[igas]
      end
    end
  end

  col_dry = [data[i,j].col_dry for i in 1:s[1], j in 1:s[2]]
  p_lay = [data[i,j].p_lay for i in 1:s[1], j in 1:s[2]]
  t_lay = [data[i,j].t_lay for i in 1:s[1], j in 1:s[2]]
  t_sfc = [data[i].t_sfc for i in 1:s[1]]

  p_lev,t_lev = ntuple(i->Array{FT}(undef, ncol,nlay+1),2)
  for i in 1:s[1]
    for j in 1:s[2]
      p_lev[i,j] = data[i,j].p_lev[1]
      t_lev[i,j] = data[i,j].t_lev[1]
      if j==s[2]
        p_lev[i,j+1] = data[i,j].p_lev[2]
        t_lev[i,j+1] = data[i,j].t_lev[2]
      end
    end
  end

  return AtmosphericState{FT,I}(gas_conc,
                                p_lay,
                                p_lev,
                                t_lay,
                                t_lev,
                                t_sfc,
                                col_gas,
                                col_dry,
                                top_at_1,
                                nlay, ncol)
end

function Base.convert(::Type{Array{AtmosphericStatePGP}}, data::AtmosphericState{FT,I}) where {FT,I}
  s = size(data.p_lay)
  nlay, ncol = size(data.p_lay)
  top_at_1 = data.top_at_1
  gas_conc = convert(Array{GasConcsPGP}, data.gas_conc)
  [AtmosphericStatePGP{FT,I}(gas_conc[i,j],
                             data.p_lay[i,j],
                             [data.p_lev[i,j], data.p_lev[i,j+1]],
                             data.t_lay[i,j],
                             [data.t_lev[i,j], data.t_lev[i,j+1]],
                             data.t_sfc[i],
                             data.col_gas[i,j,:],
                             data.col_dry[i,j],
                             data.top_at_1) for i in 1:s[1], j in 1:s[2]]

end

end #module
