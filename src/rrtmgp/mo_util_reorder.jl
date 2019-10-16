module mo_util_reorder

export reorder123x312!, reorder123x321!

"""
    reorder123x312!

Remap indexes via: (x,y,z) -> (z,x,y)
"""
function reorder123x312!(array::Array{FT}, array_out::Array{FT}) where FT
  reorder_123x312_kernel!(array, array_out)
end
"""
    reorder123x321!

Remap indexes via: (x,y,z) -> (z,y,x)
"""
function reorder123x321!(array::Array{FT}, array_out::Array{FT}) where FT
  reorder_123x321_kernel!(array, array_out)
end

function reorder_123x312_kernel!(array_in, array_out)
  for i2 = 1:size(array_in, 2)
    for i1 = 1:size(array_in, 1)
      for i3 = 1:size(array_in, 3)
        array_out[i3,i1,i2] = array_in[i1,i2,i3]
      end
    end
  end
end

function reorder_123x321_kernel!(array_in, array_out)
  for i1 = 1:size(array_in, 1)
    for i2 = 1:size(array_in, 2)
      for i3 = 1:size(array_in, 3)
        array_out[i3,i2,i1] = array_in[i1,i2,i3]
      end
    end
  end
end

end