using Profile
include("rfmip_clear_sky_lw.jl")
include("DataSetFiles.jl")

@testset "rfmip clear sky longwave driver" begin
  datafolder = RRTMGP.data_folder_rrtmgp()

  ds = dataset_dict(data_files_dict(datafolder, "lw"))

  Δt_all = Dict()

  rfmip_clear_sky_lw(ds, ty_optical_props_1scl; compile_first=true)
  rfmip_clear_sky_lw(ds, ty_optical_props_2str; compile_first=true)

  Δt_all["clear_sky_lw", "1scl"] = @elapsed rfmip_clear_sky_lw(ds, ty_optical_props_1scl)
  Δt_all["clear_sky_lw", "2str"] = @elapsed rfmip_clear_sky_lw(ds, ty_optical_props_2str)

  for (case,Δt) in Δt_all
    @show case, Δt
  end

  close.(values(ds))

end
