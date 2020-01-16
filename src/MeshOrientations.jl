"""
    MeshOrientations

Mesh orientation.
"""
module MeshOrientations

using DocStringExtensions

export MeshOrientation
export nhat, binary
export lev_range, lev_range_reversed
export ilev_top, ilev_bot, ilay_top, ilay_bot

"""
    MeshOrientation{FT}

Mesh orientation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct MeshOrientation{I}
  "Indicates whether arrays are ordered in the vertical with 1 at the top or the bottom of the domain."
  top_at_1::Bool
  "Number of layers"
  nlay::I
  MeshOrientation(top_at_1::Bool, nlay::I) where {I<:Int} = new{I}(top_at_1, nlay)
end

""" i-th level at the top of the domain"""
ilev_top(mo::MeshOrientation) = mo.top_at_1 ? 1 : mo.nlay+1

"""i-th level at the bottom of the domain"""
ilev_bot(mo::MeshOrientation) = !mo.top_at_1 ? 1 : mo.nlay+1

"""i-th layer at the top of the domain"""
ilay_top(mo::MeshOrientation) =  mo.top_at_1 ? 1 : mo.nlay

"""i-th layer at the bottom of the domain"""
ilay_bot(mo::MeshOrientation) = !mo.top_at_1 ? 1 : mo.nlay

"""
    lev_range(mo::MeshOrientation)

"""
lev_range(mo::MeshOrientation) = mo.top_at_1 ? (2:mo.nlay+1) : (mo.nlay:-1:1)

"""
    lev_range_reversed(mo::MeshOrientation)

"""
lev_range_reversed(mo::MeshOrientation) = mo.top_at_1 ? (mo.nlay:-1:1) : (2:mo.nlay+1)

binary(::Val{false}) = 0
binary(::Val{true}) = 1
binary(b::Bool) = binary(Val(b))
binary(mo::MeshOrientation) = binary(mo.top_at_1)

nhat(::Val{false}) = -1
nhat(::Val{true}) = 1
nhat(b::Bool) = nhat(Val(b))
nhat(mo::MeshOrientation) = nhat(mo.top_at_1)

end #module