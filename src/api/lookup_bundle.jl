#####
##### Typed container for the lookup tables and index maps a solver carries
#####

import Adapt
import Serialization

"""
    LookupBundle

Typed, immutable bundle of everything `RRTMGPSolver` needs from the lookup
artifacts: the gas/cloud/aerosol lookup tables for each band (absent entries
are `nothing`), the gas and aerosol name→index maps, and the band/gas counts.
Built by [`lookup_tables`](@ref); pass a prebuilt bundle back to the
`RRTMGPSolver` constructor via `lookups = ...` to avoid a second NetCDF read,
and use [`save_lookup_tables`](@ref)/[`load_lookup_tables`](@ref) to cache it
on disk (e.g. for standalone use without NCDatasets).

# Fields
- `lookup_lw`, `lookup_sw`: gas-optics lookup tables (`nothing` for gray).
- `lookup_lw_cld`, `lookup_sw_cld`: cloud-optics tables (all-sky methods only).
- `lookup_lw_aero`, `lookup_sw_aero`: aerosol tables (`aerosol_radiation` only).
- `idx_gases_lw`, `idx_gases_sw`: gas name → index maps.
- `idx_aerosol_lw`, `idx_aerosol_sw`, `idx_aerosize_lw`, `idx_aerosize_sw`:
  aerosol name → index and size-bin maps.
- `nbnd_lw`, `nbnd_sw`: band counts (1 for gray).
- `ngas_lw`, `ngas_sw`: gas counts (0 for gray).
"""
struct LookupBundle{LLW, LSW, LLWC, LSWC, LLWA, LSWA, IGL, IGS, IAL, IAS, ZL, ZS}
    lookup_lw::LLW
    lookup_sw::LSW
    lookup_lw_cld::LLWC
    lookup_sw_cld::LSWC
    lookup_lw_aero::LLWA
    lookup_sw_aero::LSWA
    idx_gases_lw::IGL
    idx_gases_sw::IGS
    idx_aerosol_lw::IAL
    idx_aerosol_sw::IAS
    idx_aerosize_lw::ZL
    idx_aerosize_sw::ZS
    nbnd_lw::Int
    nbnd_sw::Int
    ngas_lw::Int
    ngas_sw::Int
end
Adapt.@adapt_structure LookupBundle

function LookupBundle(;
    lookup_lw = nothing,
    lookup_sw = nothing,
    lookup_lw_cld = nothing,
    lookup_sw_cld = nothing,
    lookup_lw_aero = nothing,
    lookup_sw_aero = nothing,
    idx_gases_lw = nothing,
    idx_gases_sw = nothing,
    idx_aerosol_lw = nothing,
    idx_aerosol_sw = nothing,
    idx_aerosize_lw = nothing,
    idx_aerosize_sw = nothing,
    nbnd_lw::Int,
    nbnd_sw::Int,
    ngas_lw::Int = 0,
    ngas_sw::Int = 0,
)
    return LookupBundle(
        lookup_lw,
        lookup_sw,
        lookup_lw_cld,
        lookup_sw_cld,
        lookup_lw_aero,
        lookup_sw_aero,
        idx_gases_lw,
        idx_gases_sw,
        idx_aerosol_lw,
        idx_aerosol_sw,
        idx_aerosize_lw,
        idx_aerosize_sw,
        nbnd_lw,
        nbnd_sw,
        ngas_lw,
        ngas_sw,
    )
end

"""
    save_lookup_tables(path, lookups::LookupBundle)

Serialize `lookups` to `path` (host-side copies of any device arrays), so a
later session can [`load_lookup_tables`](@ref) without NCDatasets — e.g. for
standalone/classroom use of the spectral methods, or to skip the NetCDF read.
Uses Julia's `Serialization` stdlib: the file is tied to the Julia version and
package layout that wrote it and is a *cache*, not an interchange format — the
NetCDF artifacts remain the source of truth. Returns `path`.
"""
function save_lookup_tables(path::AbstractString, lookups::LookupBundle)
    Serialization.serialize(path, Adapt.adapt(Array, lookups))
    return path
end

"""
    load_lookup_tables(path, grid_params::RRTMGPGridParams)

Load a [`LookupBundle`](@ref) previously written by
[`save_lookup_tables`](@ref) and move its tables to the device described by
`grid_params`. See `save_lookup_tables` for the cache-format caveats.
"""
function load_lookup_tables(path::AbstractString, grid_params::RRTMGPGridParams)
    bundle = Serialization.deserialize(path)
    bundle isa LookupBundle || error(
        "$(path) did not deserialize to a LookupBundle; it was not written by \
         save_lookup_tables (or was written by an incompatible version).",
    )
    DA = ClimaComms.array_type(grid_params)
    return Adapt.adapt(DA, bundle)
end
