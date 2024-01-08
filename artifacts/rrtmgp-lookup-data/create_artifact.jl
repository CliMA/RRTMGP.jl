using Pkg.Artifacts, Downloads

rrtmgp_repo = "https://github.com/earth-system-radiation/rte-rrtmgp"
rrtmgp_commit_sha1 = "ed5b0113109fcd23a010a90c61f21bad551146ef" # v1.0


artifact_tree_sha1 = create_artifact() do dir
    Downloads.download(
        "$rrtmgp_repo/raw/$rrtmgp_commit_sha1/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc",
        joinpath(dir, "clearsky_lw.nc"),
    )
    Downloads.download(
        "$rrtmgp_repo/raw/$rrtmgp_commit_sha1/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc",
        joinpath(dir, "clearsky_sw.nc"),
    )
    Downloads.download(
        "$rrtmgp_repo/raw/$rrtmgp_commit_sha1/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc",
        joinpath(dir, "cloudysky_lw.nc"),
    )
    Downloads.download(
        "$rrtmgp_repo/raw/$rrtmgp_commit_sha1/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc",
        joinpath(dir, "cloudysky_sw.nc"),
    )

    cp(joinpath(@__DIR__, "README.md"), joinpath(dir, "README.md"))
end

artifact_path = joinpath(@__DIR__, "rrtmgp-lookup-data.tar.gz")

artifact_archive_sha256 = archive_artifact(artifact_tree_sha1, artifact_path)

@info "Created artifact" artifact_path artifact_tree_sha1 artifact_archive_sha256
