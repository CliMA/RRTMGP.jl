module mo_util_string

using ..fortran_intrinsics

export string_loc_in_array

"""
  string_loc_in_array

Location of string `s` in the array `array`. If the
string is not in the array, returns -1.
"""
function string_loc_in_array(s, array)
  s_loc_in_array = -1
  lc_string = lowercase(trim(s))
  for i in eachindex(array)
    if lc_string == lowercase(trim(array[i]))
      s_loc_in_array = i
      break
    end
  end
  return s_loc_in_array
end

end