module RRTMGPNCDatasetsExt

import NCDatasets as NC
using ClimaComms
using RRTMGP
using RRTMGP: RRTMGPGridParams
import RRTMGP: RRTMGPSolver
import RRTMGP: AbstractRRTMGPMethod
import RRTMGP: GrayRadiation, ClearSkyRadiation, AllSkyRadiation, AllSkyRadiationWithClearSkyDiagnostics
import RRTMGP: lookup_tables

lookup_tables(grid_params::RRTMGPGridParams, radiation_method::AbstractRRTMGPMethod) =
    lookup_tables(radiation_method, ClimaComms.device(grid_params), eltype(grid_params))

lookup_tables(radiation_method::GrayRadiation, device::ClimaComms.AbstractDevice, ::Type{FT}) where {FT} =
    (; lookups = (;), lu_kwargs = (; nbnd_lw = 1, nbnd_sw = 1))

function lookup_tables(radiation_method::ClearSkyRadiation, device::ClimaComms.AbstractDevice, ::Type{FT}) where {FT}
    DA = ClimaComms.array_type(device)

    # Call functions in lookup_constructors.jl
    artifact(t, b, n) =
        NC.Dataset(RRTMGP.ArtifactPaths.get_lookup_filename(t, b)) do ds
            getproperty(RRTMGP.LookUpTables, n)(ds, FT, DA)
        end
    lookup_lw, idx_gases_lw = artifact(:gas, :lw, :LookUpLW)

    nbnd_lw = RRTMGP.LookUpTables.get_n_bnd(lookup_lw)
    ngas_lw = RRTMGP.LookUpTables.get_n_gases(lookup_lw)

    lookup_lw_aero, idx_aerosol_lw, idx_aerosize_lw = if radiation_method.aerosol_radiation
        artifact(:aerosol, :lw, :LookUpAerosolMerra)
    else
        (nothing, nothing, nothing)
    end

    lookup_sw, idx_gases_sw = artifact(:gas, :sw, :LookUpSW)
    @assert sort(collect(keys(idx_gases_sw))) == sort(RRTMGP.gas_names_sw())

    nbnd_sw = RRTMGP.LookUpTables.get_n_bnd(lookup_sw)
    ngas_sw = RRTMGP.LookUpTables.get_n_gases(lookup_sw)

    lookup_sw_aero, idx_aerosol_sw, idx_aerosize_sw = if radiation_method.aerosol_radiation
        artifact(:aerosol, :sw, :LookUpAerosolMerra)
    else
        (nothing, nothing, nothing)
    end
    @assert sort(collect(keys(idx_aerosol_sw))) == sort(RRTMGP.aerosol_names())

    lookups = (;
        idx_aerosize_lw,
        idx_aerosize_sw,
        idx_aerosol_lw,
        idx_aerosol_sw,
        idx_gases_lw,
        idx_gases_sw,
        lookup_lw,
        lookup_lw_aero,
        lookup_sw,
        lookup_sw_aero,
    )

    @assert RRTMGP.LookUpTables.get_n_gases(lookup_lw) == RRTMGP.LookUpTables.get_n_gases(lookup_sw)
    @assert lookup_lw.p_ref_min == lookup_sw.p_ref_min
    return (; lookups, lu_kwargs = (; nbnd_lw, ngas_lw, nbnd_sw, ngas_sw))
end

function lookup_tables(
    radiation_method::Union{AllSkyRadiation, AllSkyRadiationWithClearSkyDiagnostics},
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
    clear_sky_lookups = lookup_tables(ClearSkyRadiation(radiation_method.aerosol_radiation), device, FT)
    lookups = (; clear_sky_lookups.lookups..., lookup_lw_cld, lookup_sw_cld)
    (; lu_kwargs) = clear_sky_lookups

    return (; lookups, lu_kwargs)
end

include("lookup_constructors.jl")

end # module
