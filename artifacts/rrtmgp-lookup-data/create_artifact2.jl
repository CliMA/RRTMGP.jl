using Pkg.Artifacts, Downloads

rrtmgp_commit_sha1 = "df02975ab93165b34a59f0d04b4ae6148fe5127c"
debug = false
rrtmgp_data(filename, subpath = "") =
    "https://github.com/earth-system-radiation/rrtmgp-data/raw/$rrtmgp_commit_sha1/$subpath$filename"

artifact_tree_sha1 = create_artifact() do dir
    for f in [
        "rrtmgp-aerosols-merra-lw.nc",
        "rrtmgp-aerosols-merra-sw.nc",
        "rrtmgp-clouds-lw.nc",
        "rrtmgp-clouds-sw.nc",
        "rrtmgp-gas-lw-g128.nc",
        "rrtmgp-gas-lw-g256.nc",
        "rrtmgp-gas-sw-g112.nc",
        "rrtmgp-gas-sw-g224.nc",
    ]
        subpath = ""
        debug && @info "Downloading $(rrtmgp_data(f, subpath))..."
        Downloads.download(rrtmgp_data(f, subpath), joinpath(dir, f))
    end
    for f in [
        "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc",
        "rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
        "rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
        "rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
        "rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
    ]
        subpath = "examples/rfmip-clear-sky/inputs/"
        debug && @info "Downloading $(rrtmgp_data(f, subpath))..."
        Downloads.download(rrtmgp_data(f, subpath), joinpath(dir, f))
    end
    for f in [
        "rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
        "rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
        "rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
        "rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
    ]
        subpath = "examples/rfmip-clear-sky/reference/"
        debug && @info "Downloading $(rrtmgp_data(f, subpath))..."
        Downloads.download(rrtmgp_data(f, subpath), joinpath(dir, f))
    end

    for f in [
        "rrtmgp-allsky-lw-no-aerosols.nc",
        "rrtmgp-allsky-lw.nc",
        "rrtmgp-allsky-sw-no-aerosols.nc",
        "rrtmgp-allsky-sw.nc",
    ]
        subpath = "examples/all-sky/reference/"
        debug && @info "Downloading $(rrtmgp_data(f, subpath))..."
        Downloads.download(rrtmgp_data(f, subpath), joinpath(dir, f))
    end

    cp(joinpath(@__DIR__, "README.md"), joinpath(dir, "README.md"))
end

artifact_path = joinpath(@__DIR__, "rrtmgp-lookup-data.tar.gz")

artifact_archive_sha256 = archive_artifact(artifact_tree_sha1, artifact_path)

@info "Created artifact" artifact_path artifact_tree_sha1 artifact_archive_sha256
