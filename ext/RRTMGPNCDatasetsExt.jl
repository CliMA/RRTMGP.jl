module RRTMGPNCDatasetsExt

import NCDatasets as NC
using ClimaComms
using RRTMGP
using RRTMGP: RRTMGPGridParams
import RRTMGP: RRTMGPSolver
import RRTMGP: AbstractRRTMGPMethod
import RRTMGP:
    GrayRadiation,
    ClearSkyRadiation,
    AllSkyRadiation,
    AllSkyRadiationWithClearSkyDiagnostics
import RRTMGP: lookup_tables, LookupBundle

# (The gray method needs no tables and lives in src/api/solver.jl; only the
# spectral methods below require NCDatasets.)
lookup_tables(
    grid_params::RRTMGPGridParams,
    radiation_method::AbstractRRTMGPMethod,
) = lookup_tables(
    radiation_method,
    ClimaComms.device(grid_params),
    eltype(grid_params),
)

function lookup_tables(
    radiation_method::ClearSkyRadiation,
    device::ClimaComms.AbstractDevice,
    ::Type{FT},
) where {FT}
    DA = ClimaComms.array_type(device)

    # Call functions in lookup_constructors.jl
    artifact(t, b, n) =
        NC.Dataset(RRTMGP.ArtifactPaths.get_lookup_filename(t, b)) do ds
            getproperty(RRTMGP.LookUpTables, n)(ds, FT, DA)
        end
    lookup_lw, idx_gases_lw = artifact(:gas, :lw, :LookUpLW)

    nbnd_lw = RRTMGP.LookUpTables.get_n_bnd(lookup_lw)
    ngas_lw = RRTMGP.LookUpTables.get_n_gases(lookup_lw)

    lookup_lw_aero, idx_aerosol_lw, idx_aerosize_lw =
        if radiation_method.aerosol_radiation
            artifact(:aerosol, :lw, :LookUpAerosolMerra)
        else
            (nothing, nothing, nothing)
        end

    lookup_sw, idx_gases_sw = artifact(:gas, :sw, :LookUpSW)
    @assert sort(collect(keys(idx_gases_sw))) == sort(RRTMGP.gas_names_sw()) "shortwave lookup gas names do not match RRTMGP.gas_names_sw()"

    nbnd_sw = RRTMGP.LookUpTables.get_n_bnd(lookup_sw)
    ngas_sw = RRTMGP.LookUpTables.get_n_gases(lookup_sw)

    lookup_sw_aero, idx_aerosol_sw, idx_aerosize_sw =
        if radiation_method.aerosol_radiation
            artifact(:aerosol, :sw, :LookUpAerosolMerra)
        else
            (nothing, nothing, nothing)
        end
    if !isnothing(idx_aerosol_sw)
        # Lock the loaded aerosol ordering to the canonical map (the single
        # source of truth in src/api/aerosols.jl).
        @assert idx_aerosol_sw == RRTMGP.aerosol_index_map() "shortwave aerosol lookup ordering does not match RRTMGP.aerosol_index_map() (src/api/aerosols.jl)"
        @assert idx_aerosol_lw == RRTMGP.aerosol_index_map() "longwave aerosol lookup ordering does not match RRTMGP.aerosol_index_map() (src/api/aerosols.jl)"
    end

    @assert RRTMGP.LookUpTables.get_n_gases(lookup_lw) ==
            RRTMGP.LookUpTables.get_n_gases(lookup_sw)
    @assert lookup_lw.p_ref_min == lookup_sw.p_ref_min
    return LookupBundle(;
        lookup_lw,
        lookup_sw,
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

function lookup_tables(
    radiation_method::Union{
        AllSkyRadiation,
        AllSkyRadiationWithClearSkyDiagnostics,
    },
    device::ClimaComms.AbstractDevice,
    ::Type{FT},
) where {FT}
    DA = ClimaComms.array_type(device)

    # Call functions in lookup_constructors.jl
    artifact(t, b, n) =
        NC.Dataset(RRTMGP.ArtifactPaths.get_lookup_filename(t, b)) do ds
            getproperty(RRTMGP.LookUpTables, n)(ds, FT, DA)
        end
    lookup_lw_cld = artifact(:cloud, :lw, :LookUpCld)
    lookup_sw_cld = artifact(:cloud, :sw, :LookUpCld)
    b = lookup_tables(
        ClearSkyRadiation(radiation_method.aerosol_radiation),
        device,
        FT,
    )
    return LookupBundle(;
        b.lookup_lw,
        b.lookup_sw,
        lookup_lw_cld,
        lookup_sw_cld,
        b.lookup_lw_aero,
        b.lookup_sw_aero,
        b.idx_gases_lw,
        b.idx_gases_sw,
        b.idx_aerosol_lw,
        b.idx_aerosol_sw,
        b.idx_aerosize_lw,
        b.idx_aerosize_sw,
        b.nbnd_lw,
        b.nbnd_sw,
        b.ngas_lw,
        b.ngas_sw,
    )
end

include("lookup_constructors.jl")

end # module
