"""
    ReferenceStates

Reference state variables for look-up tables / interpolation
"""
module ReferenceStates

using DocStringExtensions
using ..Gases
using ..Utilities
export ReferenceState
export get_press_min

abstract type AbstractReferenceState{FT<:AbstractFloat} end

"""
    ReferenceState{FT}

Reference state variables for look-up tables / interpolation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ReferenceState{FT} <: AbstractReferenceState{FT}
  "Pressure"
  press::Vector{FT}
  "Log of pressure"
  press_log::Vector{FT}
  "Temperature"
  temp::Vector{FT}
  "Minimum of pressure"
  press_min::FT
  "Maximum of pressure"
  press_max::FT
  "Minimum of temperature"
  temp_min::FT
  "Maximum of temperature"
  temp_max::FT
  "Maximum of temperature"
  press_log_delta::FT
  temp_delta::FT
  press_trop_log::FT
  vmr::Array{FT,3}   # vmr(lower or upper atmosphere, gas, temp)
  function ReferenceState(press::Array{FT},
                          temp::Array{FT},
                          press_ref_trop::FT,
                          vmr_ref::Array{FT},
                          gases_prescribed::Vector{AbstractGas},
                          gases_in_database::Vector{AbstractGas}) where {FT<:AbstractFloat}

    gas_is_present = map(x->x in gases_prescribed, gases_in_database)
    ngas = count(gas_is_present)

    press_log = log.(press)
    # TODO: remove assumption of ordering
    temp_min = temp[1]
    temp_max = temp[length(temp)]
    press_min = press[length(press)]
    press_max = press[1]

    press_trop_log = log(press_ref_trop)

    press_log_delta = (log(press_min)-log(press_max))/(length(press)-1)
    temp_delta      = (temp_max-temp_min)/(length(temp)-1)

    vmr_ref_red = Array{FT}(undef, size(vmr_ref, 1),ngas+1, size(vmr_ref, 3))

    # Gas 0 is used in single-key species method, set to 1.0 (col_dry)
    vmr_ref_red[:,1,:] = vmr_ref[:,1,:]
    for i = 1:ngas
      idx = loc_in_array(gases_prescribed[i], gases_in_database)
      vmr_ref_red[:,i+1,:] = vmr_ref[:,idx+1,:]
    end
    return new{FT}(press,
                   press_log,
                   temp,
                   press_min,
                   press_max,
                   temp_min,
                   temp_max,
                   press_log_delta,
                   temp_delta,
                   press_trop_log,
                   vmr_ref_red)
  end
end

"""
    get_press_min(ref::ReferenceState)

Minimum pressure on the interpolation grids
"""
get_press_min(ref::ReferenceState) = ref.press_min

end