"""
    MeshOrientations

Mesh orientation.
"""
module MeshOrientations

using DocStringExtensions

export MeshOrientation

"""
    MeshOrientation{FT}

Mesh orientation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct MeshOrientation{I}
  "Indicates whether arrays are ordered in the vertical with 1 at the top or the bottom of the domain."
  top_at_1::Bool
  "i-th level at the top of the domain"
  ilev_top::I
  "i-th level at the bottom of the domain"
  ilev_bot::I
  "i-th layer at the top of the domain"
  ilay_top::I
  "i-th layer at the bottom of the domain"
  ilay_bot::I
  function MeshOrientation(top_at_1::Bool, nlay::I) where {I<:Int}
    ilev_top =  top_at_1 ? 1 : nlay+1
    ilev_bot = !top_at_1 ? 1 : nlay+1

    ilay_top =  top_at_1 ? 1 : nlay
    ilay_bot = !top_at_1 ? 1 : nlay

    return new{I}(top_at_1, ilev_top, ilev_bot, ilay_top, ilay_bot)
  end
end

end #module