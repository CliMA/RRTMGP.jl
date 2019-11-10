using Profile
include("rfmip_clear_sky_sw.jl")
include("DataSetFiles.jl")

@testset "rfmip clear sky shortwave driver" begin
  datafolder = JRRTMGP.data_folder_rrtmgp()

  ds = dataset_dict(data_files_dict(datafolder, "sw"))

  rfmip_clear_sky_sw(ds, ty_optical_props_1scl; compile_first=true)
  rfmip_clear_sky_sw(ds, ty_optical_props_2str; compile_first=true)

  rfmip_clear_sky_sw(ds, ty_optical_props_1scl)
  rfmip_clear_sky_sw(ds, ty_optical_props_2str)

  close.(values(ds))

end
