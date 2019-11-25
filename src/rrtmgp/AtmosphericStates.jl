"""
    AtmosphericStates

Atmospheric conditions used as inputs to RRTMGP.
"""
module AtmosphericStates

using DocStringExtensions

using ..GasConcentrations

export AtmosphericState, update_view!

"""
    SimpleGrid{FT}

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SimpleGrid{FT}
  z
end


"""
    AtmosphericVars{FT}

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AtmosphericVars{FT}
  "Gas concentrations, in the form of volume mixing ratios"
  gas_concs
  "layer pressures [Pa, mb]; (ncol,nlay)"
  p_lay
  "level pressures [Pa, mb]; (ncol,nlay+1)"
  p_lev
  "layer temperatures [K]; (ncol,nlay)"
  t_lay
  "level temperatures [K]; (ncol,nlay+1)"
  t_lev
  function AtmosphericVars(gas_concs,p_lay,p_lev,t_lay,t_lev=nothing)
    # if t_lev == nothing
    #   tlev_arr = zeros(FT, ncol,nlay+1)

    #   # Interpolate temperature to levels if not provided
    #   #   Interpolation and extrapolation at boundaries is weighted by pressure
    #   #
    #   for icol = 1:ncol
    #      tlev_arr[icol,1] = t_lay[icol,1] +
    #                        (p_lev[icol,1]-p_lay[icol,1])*
    #                        (t_lay[icol,2]-t_lay[icol,1])/
    #                        (p_lay[icol,2]-p_lay[icol,1])
    #   end
    #   for ilay in 2:nlay
    #     for icol in 1:ncol
    #        tlev_arr[icol,ilay] = (p_lay[icol,ilay-1]*t_lay[icol,ilay-1]*
    #                              (p_lev[icol,ilay]-p_lay[icol,ilay]) +
    #                               p_lay[icol,ilay]*t_lay[icol,ilay]*
    #                              (p_lay[icol,ilay-1]-p_lev[icol,ilay]))/
    #                              (p_lev[icol,ilay]*(p_lay[icol,ilay-1] - p_lay[icol,ilay]))
    #     end
    #   end
    #   for icol = 1:ncol
    #      tlev_arr[icol,nlay+1] = t_lay[icol,nlay] +
    #                             (p_lev[icol,nlay+1]-p_lay[icol,nlay])*
    #                             (t_lay[icol,nlay]-t_lay[icol,nlay-1])/
    #                             (p_lay[icol,nlay]-p_lay[icol,nlay-1])
    #   end
    #   t_lev = tlev_arr
    # end
    return new{eltype(p_lay)}(gas_concs,p_lay,p_lev,t_lay,t_lev)
  end
end


"""
    AtmosphericState{FT}

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct AtmosphericState{FT}
  grid
  _full::AtmosphericVars
  _view::AtmosphericVars
  "Indicates whether arrays are ordered in the vertical with 1 at the top or the bottom of the domain."
  top_at_1
  nlay
  ncol
  function AtmosphericState(gas_concs,p_lay,p_lev,t_lay,t_lev=nothing)
    _full = AtmosphericVars(gas_concs,p_lay,p_lev,t_lay,t_lev)
    _view = _full.gas_concs
    _view = _full.p_lay
    _view = _full.p_lev
    _view = _full.t_lay
    _view = _full.t_lev

    FT = eltype(p_lay)
    nlay = size(p_lay, 2)
    ncol = size(p_lay, 1)
    top_at_1 = p_lay[1, 1, 1] < p_lay[1, nlay, 1]
    grid = SimpleGrid{FT}(p_lay)
    return new{FT}(grid,_full,_view,top_at_1,nlay,ncol)
  end
end

function update_view!(atmos_state::AtmosphericState, b::Int)
    atmos_state.p_lay     = @view(atmos_state.p_lay_full[:,:,b])
    atmos_state.p_lev     = @view(atmos_state.p_lev_full[:,:,b])
    atmos_state.t_lay     = @view(atmos_state.t_lay_full[:,:,b])
    atmos_state.t_lev     = @view(atmos_state.t_lev_full[:,:,b])
    atmos_state.gas_concs = @view(atmos_state.gas_concs_full[b])
end


end #module
