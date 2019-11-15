module mo_util_string

using ..fortran_intrinsics

export string_loc_in_array

"""
  string_loc_in_array

Location of string `s` in the array `array`. If the
string is not in the array, returns -1.
"""
string_loc_in_array(s, array) = s in array ? findfirst(i->i==s,array) : -1

end