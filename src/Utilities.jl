module Utilities

export loc_in_array, check_range, check_extent

"""
  loc_in_array

Location of element `e` in the array `array`. If the
element is not in the array, returns -1.
"""
loc_in_array(e, array) = e in array ? findfirst(i->i==e,array) : -1

"""
    check_extent(array, s, label)

Assert array has correct size
"""
check_extent(array, s, label) = @assert all(size(array).==s)

"""
    check_range(val, minV, maxV, label)

Assert range of values is valid
"""
check_range(val, minV, maxV, label) = any(val .< minV) || any(val .> maxV) ? error(strip(label) * " values out of range.") : nothing

end