#####
##### Canonical MERRA aerosol species ordering
#####
#
# This is the single source of truth for the order of aerosol species in the
# `AerosolState` arrays (`aero_mass`, `aero_size`), i.e., the meaning of the
# first index of those arrays. It MUST agree with:
#
#  - the integer indices hard-coded by the aerosol optics kernel
#    `compute_lookup_aerosol` in `src/optics/aerosol_optics.jl`, and
#  - the `idx_aerosol` map built when the MERRA aerosol lookup table is loaded
#    (`ext/lookup_constructors.jl`, `LookUpAerosolMerra`).
#
# Agreement is locked by `test/aerosol_name_map.jl` (and, when the lookup table
# is loaded, by an assertion in `ext/RRTMGPNCDatasetsExt.jl`). Do not reorder
# these without updating the kernel indices.

const AEROSOL_IDX = Dict{String, Int}(
    "dust1" => 1,
    "sea_salt1" => 2,
    "sulfate" => 3,
    "black_carbon_rh" => 4,
    "black_carbon" => 5,
    "organic_carbon_rh" => 6,
    "organic_carbon" => 7,
    "dust2" => 8,
    "dust3" => 9,
    "dust4" => 10,
    "dust5" => 11,
    "sea_salt2" => 12,
    "sea_salt3" => 13,
    "sea_salt4" => 14,
    "sea_salt5" => 15,
)

# Aerosol species names, ordered by their index (so `_AEROSOL_NAMES[i]` is the
# species stored at index `i` of the aerosol state arrays).
const _AEROSOL_NAMES = [
    "dust1",
    "sea_salt1",
    "sulfate",
    "black_carbon_rh",
    "black_carbon",
    "organic_carbon_rh",
    "organic_carbon",
    "dust2",
    "dust3",
    "dust4",
    "dust5",
    "sea_salt2",
    "sea_salt3",
    "sea_salt4",
    "sea_salt5",
]

"""
    aerosol_index_map()

Return the canonical mapping from aerosol species name to its index in the
`AerosolState` arrays (`aero_mass`, `aero_size`). This ordering is the single
source of truth shared by the optics kernel, the lookup-table loader, and the
state accessors. Returns a fresh copy — mutating it does not affect RRTMGP.
"""
aerosol_index_map() = copy(AEROSOL_IDX)

"""
    aerosol_names()

Return the canonical aerosol species names, ordered by their index, so that
`aerosol_names()[i]` is the species stored at index `i` of the `AerosolState`
arrays. Returns a fresh copy — mutating it does not affect RRTMGP.
"""
aerosol_names() = copy(_AEROSOL_NAMES)

aerosol_names_docs() = map(x -> "`$x`", aerosol_names())

"""
    canonical_aerosol_name(name)

Resolve an aerosol `name` to its canonical name. Throws if `name` is unknown.
"""
function canonical_aerosol_name(name::AbstractString)
    haskey(AEROSOL_IDX, name) && return String(name)
    return error(
        "unknown aerosol name $(repr(name)); known names are " *
        "$(sort(collect(keys(AEROSOL_IDX))))",
    )
end

"""
    aerosol_index(name)

Return the index, into the `AerosolState` arrays, of aerosol `name`.
"""
aerosol_index(name::AbstractString) = AEROSOL_IDX[canonical_aerosol_name(name)]
