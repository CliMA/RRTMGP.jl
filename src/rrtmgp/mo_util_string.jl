module mo_util_string

using ..fortran_intrinsics

export string_in_array, string_loc_in_array

"""
  string_in_array

Checks if a string `s` is in the array `array`.
"""
string_in_array(s, array) = s in array

"""
  string_in_array

Location of string `s` in the array `array`. If the
string is not in the array, returns -1.
"""
function string_loc_in_array(s, array)
  i = findfirst(s,array)
  i==nothing ? -1 : i
end

end