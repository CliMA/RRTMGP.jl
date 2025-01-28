module ArtifactPaths

using Artifacts
import ..get_artifact_path

export get_lookup_filename, get_input_filename

"""
    get_lookup_filename(optics_type::Symbol, λ::Symbol)

This function generates the file names for lookup table files, for a given optics type, for the
shortwave and longwave solvers.

- `:gas`, `:cloud` and `:aerosol` optics types are supported for the longwave and shortwave solvers. 
  - The `:gas` option provides the file that is used to load either the `LookUpLW` or `LookUpSW` struct in `LookUpTables.jl`.
  - The `:cloud` option provides the file that is used to load the `LookUpCld` struct in `LookUpTables.jl`.
  - The `:aerosol` option provides the file that is used to load the `LookUpAerosolMerra` struct in `LookUpTables.jl`.
- `:lw` (longwave) and `:sw` (shortwave) wavelength types are supported.

These artifacts are obtained from "Pincus, R., Mlawer, E. J., Delamere, J., Iacono, M. J., & Pernak, R. (2023). RRTMGP data (Version 1.7) [Data set]. https://github.com/earth-system-radiation/rrtmgp-data" 

The file "rrtmgp-sw-inputs-aerosol-optics.nc" overrides the lookup table available from the artifacts. This table corrects an error in the array ordering for the aerosol optics lookup table for the shortwave sea-salt data (‘aero_salt_tbl’). This file is provided by Michael Iacono at Atmospheric and Environmental Research, via personal communication. This file is expected to replace the currently existing lookup table in the `rrtmgp-data` repository in their next public release.
"""
function get_lookup_filename(optics_type::Symbol, λ::Symbol)
    @assert optics_type ∈ (:gas, :cloud, :aerosol)
    @assert λ ∈ (:lw, :sw)
    basedir = get_artifact_path()
    currdir = @__DIR__
    config = (optics_type, λ)

    config == (:gas, :lw) && return joinpath(basedir, "rrtmgp-gas-lw-g256.nc")
    config == (:gas, :sw) && return joinpath(basedir, "rrtmgp-gas-sw-g224.nc")

    config == (:cloud, :lw) && return joinpath(basedir, "rrtmgp-clouds-lw-bnd.nc")
    config == (:cloud, :sw) && return joinpath(basedir, "rrtmgp-clouds-sw-bnd.nc")

    config == (:aerosol, :lw) && return joinpath(basedir, "rrtmgp-aerosols-merra-lw.nc")
    config == (:aerosol, :sw) && return joinpath(basedir, "rrtmgp-aerosols-merra-sw.nc")
end

"""
    get_input_filename(problemtype::Symbol, λ::Symbol)

This function generates the file names for input files for tests for a given problem type and wavelength type.
`:gas`, and `:gas_clouds` and `:gas_clouds_aerosols` problem types are supported.
`:lw` (longwave) and `:sw` (shortwave) wavelength types are supported.

This file provides data for loading the `AtmosphericState` struct, the provides the atmospheric conditions for computing optical properties.

- The `:gas` option is used for the `clear sky` test.
- The `:gas_clouds` option is used for the `all sky` test.
- The `:gas_clouds_aerosols` option is used for the `all sky with aerosols` test.

While these are primarily intended for the tests, some of this input data is also used in `ClimaAtmos.jl` and is therefore provided here for user's convenience.

These artifacts are obtained from "Pincus, R., Mlawer, E. J., Delamere, J., Iacono, M. J., & Pernak, R. (2023). RRTMGP data (Version 1.7) [Data set]. https://github.com/earth-system-radiation/rrtmgp-data" 
"""
function get_input_filename(problemtype::Symbol, λ::Symbol)
    @assert problemtype ∈ (:gas, :gas_clouds, :gas_clouds_aerosols)
    @assert λ ∈ (:lw, :sw)

    basedir = get_artifact_path()

    if problemtype == :gas
        return joinpath(
            basedir,
            "examples",
            "rfmip-clear-sky",
            "inputs",
            "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc",
        )
    elseif problemtype == :gas_clouds
        dir = joinpath(basedir, "examples", "all-sky", "reference")
        return λ == :lw ? joinpath(dir, "rrtmgp-allsky-lw-no-aerosols.nc") :
               joinpath(dir, "rrtmgp-allsky-sw-no-aerosols.nc")
    else # :gas_clouds_aerosols
        dir = joinpath(basedir, "examples", "all-sky", "reference")
        return λ == :lw ? joinpath(dir, "rrtmgp-allsky-lw.nc") : joinpath(dir, "rrtmgp-allsky-sw.nc")
    end
end

end
