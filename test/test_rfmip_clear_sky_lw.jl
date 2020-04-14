using Profile
include("rfmip_clear_sky_lw.jl")
include("rfmip_clear_sky_lw_pgp.jl")
include("data_set_files.jl")

@testset "rfmip clear sky longwave driver" begin
    datafolder = RRTMGP.data_folder_rrtmgp()

    ds = dataset_dict(data_files_dict(datafolder, "lw"))

    Δt_all = Dict()

    if !pgp_only
        rfmip_clear_sky_lw(ds, OneScalar; compile_first = true)
        rfmip_clear_sky_lw(ds, TwoStream; compile_first = true)
    end

    rfmip_clear_sky_lw_pgp(ds, OneScalar; compile_first = true)
    rfmip_clear_sky_lw_pgp(ds, TwoStream; compile_first = true)

    if !pgp_only
        Δt_all["clear_sky_lw", "1scl"] = @elapsed rfmip_clear_sky_lw(
            ds,
            OneScalar,
        )
        Δt_all["clear_sky_lw", "2str"] = @elapsed rfmip_clear_sky_lw(
            ds,
            TwoStream,
        )
    end

    Δt_all["clear_sky_lw_pgp", "1scl"] = @elapsed rfmip_clear_sky_lw_pgp(
        ds,
        OneScalar,
    )
    Δt_all["clear_sky_lw_pgp", "2str"] = @elapsed rfmip_clear_sky_lw_pgp(
        ds,
        TwoStream,
    )

    for (case, Δt) in Δt_all
        @show case, Δt
    end

    close.(values(ds))

end
