"""
    ReferenceStates

Reference state variables for look-up tables / interpolation
"""
module ReferenceStates

using DocStringExtensions
using OffsetArrays
using ..Gases
using ..Utilities
export ReferenceStateProfile, ReferenceState
export get_press_min

"""
    ReferenceStateProfile{FT}

Reference state variables for look-up tables / interpolation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ReferenceStateProfile{FT}
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
  function ReferenceStateProfile(press::Array{FT},
                                 temp::Array{FT},
                                 press_ref_trop::FT) where {FT<:AbstractFloat}

    press_log = log.(press)
    # TODO: remove assumption of ordering
    temp_min = temp[1]
    temp_max = temp[length(temp)]
    press_min = press[length(press)]
    press_max = press[1]

    press_trop_log = log(press_ref_trop)

    press_log_delta = (log(press_min)-log(press_max))/(length(press)-1)
    temp_delta      = (temp_max-temp_min)/(length(temp)-1)

    return new{FT}(press,
                   press_log,
                   temp,
                   press_min,
                   press_max,
                   temp_min,
                   temp_max,
                   press_log_delta,
                   temp_delta,
                   press_trop_log)
  end
end

"""
    ReferenceStatePerGridPoint{FT}

Reference state variables for look-up tables / interpolation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ReferenceStatePerGridPoint{FT}
  "Pressure"
  press::FT
  "Log of pressure"
  press_log::FT
  "Temperature"
  temp::FT
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
end
ReferenceStatePerGridPoint(ref::ReferenceStateProfile) =
  ReferenceStatePerGridPoint.(press,
                              press_log,
                              temp,
                              Ref(press_min),
                              Ref(press_max),
                              Ref(temp_min),
                              Ref(temp_max),
                              Ref(press_log_delta),
                              Ref(temp_delta),
                              Ref(press_trop_log))

"""
    ReferenceState{FT}

Reference state variables for look-up tables / interpolation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ReferenceState{FT}
  vmr::AbstractArray{FT,3}   # vmr(lower or upper atmosphere, gas, temp)
  function ReferenceState(vmr_ref::Array{FT},
                          available_gases::Vector{AbstractGas},
                          gas_names::Vector{AbstractGas}) where {FT<:AbstractFloat}

    gas_is_present = map(x->x in available_gases, gas_names)
    ngas = count(gas_is_present)

    vmr_ref_red = OffsetArray{FT}(undef, 1:size(vmr_ref, 1),0:ngas, 1:size(vmr_ref, 3))

    # Gas 0 is used in single-key species method, set to 1.0 (col_dry)
    vmr_ref_red[:,0,:] = vmr_ref[:,1,:]
    for i = 1:ngas
      idx = loc_in_array(available_gases[i], gas_names)
      vmr_ref_red[:,i,:] = vmr_ref[:,idx+1,:]
    end
    return new{FT}(vmr_ref_red)
  end
end

"""
    get_press_min(ref::ReferenceStateProfile)

Minimum pressure on the interpolation grids
"""
get_press_min(ref::ReferenceStateProfile) = ref.press_min

end