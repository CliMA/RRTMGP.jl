# Dataset files

# This function generates the file names for lookup table files.
function get_lookup_filename(optics_type, λ)
    @assert optics_type ∈ (:gas, :cloud, :aerosol)
    @assert λ ∈ (:lw, :sw)
    dir = joinpath(artifact"rrtmgp-data", "rrtmgp-data-1.8.1")

    return optics_type == :gas ? joinpath(dir, "rrtmgp-gas-" * (λ == :lw ? "lw-g256.nc" : "sw-g224.nc")) :
           (
        optics_type == :cloud ? joinpath(dir, "rrtmgp-clouds-" * (λ == :lw ? "lw.nc" : "sw.nc")) :
        joinpath(dir, "rrtmgp-aerosols-merra-" * (λ == :lw ? "lw.nc" : "sw.nc"))
    )
end

# This function generates the file names for input files for tests.
function get_input_filename(problemtype, λ)
    @assert problemtype ∈ (:gas, :gas_clouds, :gas_clouds_aerosols)
    @assert λ ∈ (:lw, :sw)

    dir = joinpath(artifact"rrtmgp-data", "rrtmgp-data-1.8.1")

    if problemtype == :gas
        return joinpath(
            dir,
            "examples",
            "rfmip-clear-sky",
            "inputs",
            "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc",
        )
    elseif problemtype == :gas_clouds
        dir = joinpath(dir, "examples", "all-sky", "reference")
        return λ == :lw ? joinpath(dir, "rrtmgp-allsky-lw-no-aerosols.nc") :
               joinpath(dir, "rrtmgp-allsky-sw-no-aerosols.nc")
    else # :gas_clouds_aerosols
        dir = joinpath(dir, "examples", "all-sky", "reference")
        return λ == :lw ? joinpath(dir, "rrtmgp-allsky-lw.nc") : joinpath(dir, "rrtmgp-allsky-sw.nc")
    end
end

# This function generates the reference files for comparison for various tests.
function get_reference_filename(problemtype, λ, flux_up_dn)
    @assert problemtype ∈ (:gas, :gas_clouds, :gas_clouds_aerosols)
    @assert λ ∈ (:lw, :sw)
    @assert flux_up_dn ∈ (:flux_up, :flux_dn)

    dir = joinpath(artifact"rrtmgp-data", "rrtmgp-data-1.8.1")

    if problemtype == :gas
        dir = joinpath(dir, "examples", "rfmip-clear-sky", "reference")
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
        dir = joinpath(dir, "examples", "all-sky", "reference")
        return λ == :lw ? joinpath(dir, "rrtmgp-allsky-lw-no-aerosols.nc") :
               joinpath(dir, "rrtmgp-allsky-sw-no-aerosols.nc")
    else # :gas_clouds_aerosols
        dir = joinpath(dir, "examples", "all-sky", "reference")
        return λ == :lw ? joinpath(dir, "rrtmgp-allsky-lw.nc") : joinpath(dir, "rrtmgp-allsky-sw.nc")
    end
end
