using Profile
include("allsky.jl")
include("DataSetFiles.jl")

@testset "All sky" begin
  datafolder = JRRTMGP.data_folder_rrtmgp()

  ds_sw = dataset_dict(data_files_dict(datafolder, "sw"; allsky=true))
  ds_lw = dataset_dict(data_files_dict(datafolder, "lw"; allsky=true))

  all_sky(ds_lw; use_luts=false, λ_string = "lw", compile_first=true)
  all_sky(ds_sw; use_luts=false, λ_string = "sw", compile_first=true)
  all_sky(ds_lw; use_luts=true, λ_string = "lw", compile_first=true)
  all_sky(ds_sw; use_luts=true, λ_string = "sw", compile_first=true)

  all_sky(ds_lw; use_luts=false, λ_string = "lw")
  all_sky(ds_sw; use_luts=false, λ_string = "sw")
  all_sky(ds_lw; use_luts=true, λ_string = "lw")
  all_sky(ds_sw; use_luts=true, λ_string = "sw")

  close.(values(ds_lw))
  close.(values(ds_sw))

end

