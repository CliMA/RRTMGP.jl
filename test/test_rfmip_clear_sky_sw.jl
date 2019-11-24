using Profile
include("rfmip_clear_sky_sw.jl")
include("DataSetFiles.jl")

@testset "rfmip clear sky shortwave driver" begin
  datafolder = RRTMGP.data_folder_rrtmgp()

  ds = dataset_dict(data_files_dict(datafolder, "sw"))

  Δt_all = Dict()

  rfmip_clear_sky_sw(ds, OneScalar; compile_first=true)
  rfmip_clear_sky_sw(ds, TwoStream; compile_first=true)

  Δt_all["clear_sky_sw", "1scl"] = @elapsed rfmip_clear_sky_sw(ds, OneScalar)
  Δt_all["clear_sky_sw", "2str"] = @elapsed rfmip_clear_sky_sw(ds, TwoStream)

  for (case,Δt) in Δt_all
    @show case, Δt
  end

  close.(values(ds))

end
