"""
    mo_util_reorder

Routines for permuting arrays:
  - (x,y,z) -> (z,x,y)
  - (x,y,z) -> (z,y,x)
"""
module mo_util_reorder

export reorder123x312!, reorder123x321!

"""
    reorder123x312!(src, dst)

Reorder indexes using mapping: (x,y,z) -> (z,x,y)
"""
reorder123x312!(src, dst) = permutedims!(dst, src, [3,1,2])

"""
    reorder123x321!(src, dst)

Reorder indexes using mapping: (x,y,z) -> (z,y,x)
"""
reorder123x321!(src, dst) = permutedims!(dst, src, [3,2,1])

end
