using Profile
include("allsky.jl")
include("DataSetFiles.jl")

@testset "All sky" begin
  datafolder = JRRTMGP.data_folder_rrtmgp()

  ds_sw = dataset_dict(data_files_dict(datafolder, "sw"; allsky=true))
  ds_lw = dataset_dict(data_files_dict(datafolder, "lw"; allsky=true))

  Δt_all = Dict()

  all_sky(ds_lw; use_luts=false, λ_string = "lw", compile_first=true)
  all_sky(ds_sw; use_luts=false, λ_string = "sw", compile_first=true)
  all_sky(ds_lw; use_luts=true, λ_string = "lw", compile_first=true)
  all_sky(ds_sw; use_luts=true, λ_string = "sw", compile_first=true)

  Δt_all["all_sky_lw","use_luts_false"] = @elapsed all_sky(ds_lw; use_luts=false, λ_string = "lw")
  Δt_all["all_sky_sw","use_luts_false"] = @elapsed all_sky(ds_sw; use_luts=false, λ_string = "sw")
  Δt_all["all_sky_lw","use_luts_true"] = @elapsed all_sky(ds_lw; use_luts=true, λ_string = "lw")
  Δt_all["all_sky_sw","use_luts_true"] = @elapsed all_sky(ds_sw; use_luts=true, λ_string = "sw")

  for (case,Δt) in Δt_all
    @show case, Δt
  end

  close.(values(ds_lw))
  close.(values(ds_sw))

end

