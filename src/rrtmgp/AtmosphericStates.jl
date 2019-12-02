"""
    AtmosphericStates

Atmospheric conditions used as inputs to RRTMGP.
"""
module AtmosphericStates

using DocStringExtensions

using ..Utilities
using ..GasConcentrations
using ..ReferenceStates

export AtmosphericState

abstract type AbstractAtmosphericState{AbstractFloat,Integer} end

"""
    SimpleGrid{FT}

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SimpleGrid{FT}
  z
end

function interpolate_temperature(p_lay::Array{FT},
                                 p_lev::Array{FT},
                                 t_lay::Array{FT},
                                 ncol::I,
                                 nlay::I) where {FT<:AbstractFloat,I<:Int}
  t_lev = zeros(FT, ncol,nlay+1)
  # Interpolate temperature to levels if not provided
  #   Interpolation and extrapolation at boundaries is weighted by pressure
  #
  for icol = 1:ncol
     t_lev[icol,1] = t_lay[icol,1] +
                    (p_lev[icol,1]-p_lay[icol,1])*
                    (t_lay[icol,2]-t_lay[icol,1])/
                    (p_lay[icol,2]-p_lay[icol,1])
  end
  for ilay in 2:nlay
    for icol in 1:ncol
       t_lev[icol,ilay] = (p_lay[icol,ilay-1]*t_lay[icol,ilay-1]*
                          (p_lev[icol,ilay]-p_lay[icol,ilay]) +
                           p_lay[icol,ilay]*t_lay[icol,ilay]*
                          (p_lay[icol,ilay-1]-p_lev[icol,ilay]))/
                          (p_lev[icol,ilay]*(p_lay[icol,ilay-1] - p_lay[icol,ilay]))
    end
  end
  for icol = 1:ncol
     t_lev[icol,nlay+1] = t_lay[icol,nlay] +
                         (p_lev[icol,nlay+1]-p_lay[icol,nlay])*
                         (t_lay[icol,nlay]-t_lay[icol,nlay-1])/
                         (p_lay[icol,nlay]-p_lay[icol,nlay-1])
  end
  return t_lev
end

"""
    AtmosphericState{FT}

as = AtmosphericState(gas_conc,p_lay, p_lev, t_lay, t_lev)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AtmosphericState{FT,I} <: AbstractAtmosphericState{FT,I}
  "Grid over which vertical integration is performed"
  grid::SimpleGrid{FT}
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
  "Indicates whether arrays are ordered in the vertical with 1 at the top or the bottom of the domain."
  top_at_1::Bool
  nlay::I
  ncol::I
  function AtmosphericState(gas_conc::GasConcs{FT,I},
                            p_lay::Array{FT,2},
                            p_lev::Array{FT,2},
                            t_lay::Array{FT,2},
                            t_lev::Union{Array{FT,2},Nothing},
                            ref::ReferenceState) where {I<:Int,FT<:AbstractFloat}
    nlay = size(p_lay, 2)
    ncol = size(p_lay, 1)
    if t_lev==nothing
      t_lev = interpolate_temperature(p_lay, p_lev, t_lay, ncol, nlay)
    end

    check_extent(p_lay, (ncol, nlay  ), "p_lay")
    check_extent(p_lev, (ncol, nlay+1), "p_lev")
    check_extent(t_lay, (ncol, nlay  ), "t_lay")
    check_extent(t_lev, (ncol, nlay+1), "t_lev")
    check_range(p_lay, ref.press_min, ref.press_max, "p_lay")
    check_range(p_lev, ref.press_min, ref.press_max, "p_lev")
    check_range(t_lay, ref.temp_min,  ref.temp_max,  "t_lay")
    check_range(t_lev, ref.temp_min,  ref.temp_max,  "t_lev")

    # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
    top_at_1 = p_lay[1, 1] < p_lay[1, nlay]
    grid = SimpleGrid{FT}(p_lay)
    return new{FT,I}(grid,gas_conc,p_lay,p_lev,t_lay,t_lev,top_at_1,nlay,ncol)
  end
end

end #module
