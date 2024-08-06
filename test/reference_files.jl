using Artifacts
import RRTMGP: get_artifact_path

"""
    get_reference_filename(problemtype::Symbol, λ::Symbol, flux_up_dn::Symbol)

This function generates the filenames for reference files for comparison for various tests. Specifically, these files provide the upwelling (`flux_up`) and downwelling (`flux_dn`) flux data generated from [RTE-RRTMGP](https://github.com/earth-system-radiation/rte-rrtmgp) runs for various production configurations.

- `:gas`, and `:gas_clouds` and `:gas_clouds_aerosols` problem types are supported.
- `:lw` (longwave) and `:sw` (shortwave) wavelength types are supported.
- `flux_up` and `flux_dn` options are supported for flux types.

These artifacts are obtained from https://github.com/earth-system-radiation/rrtmgp-data package.
"""
function get_reference_filename(problemtype::Symbol, λ::Symbol, flux_up_dn::Symbol)
    @assert problemtype ∈ (:gas, :gas_clouds, :gas_clouds_aerosols)
    @assert λ ∈ (:lw, :sw)
    @assert flux_up_dn ∈ (:flux_up, :flux_dn)

    basedir = get_artifact_path()

    if problemtype == :gas
        dir = joinpath(basedir, "examples", "rfmip-clear-sky", "reference")
        if flux_up_dn == :flux_up
            if λ == :lw
                return joinpath(dir, "rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc")
            else # λ == :sw
                return joinpath(dir, "rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc")
            end
        else
            if λ == :lw
                return joinpath(dir, "rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc")
            else # λ == :sw
                return joinpath(dir, "rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc")
            end
        end
    elseif problemtype == :gas_clouds
        dir = joinpath(basedir, "examples", "all-sky", "reference")
        return λ == :lw ? joinpath(dir, "rrtmgp-allsky-lw-no-aerosols.nc") :
               joinpath(dir, "rrtmgp-allsky-sw-no-aerosols.nc")
    else # :gas_clouds_aerosols
        dir = joinpath(basedir, "examples", "all-sky", "reference")
        return λ == :lw ? joinpath(dir, "rrtmgp-allsky-lw.nc") : joinpath(dir, "rrtmgp-allsky-sw.nc")
    end
end
