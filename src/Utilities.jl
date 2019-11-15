module Utilities

export loc_in_array

"""
  loc_in_array

Location of element `e` in the array `array`. If the
element is not in the array, returns -1.
"""
loc_in_array(e, array) = e in array ? findfirst(i->i==e,array) : -1

end