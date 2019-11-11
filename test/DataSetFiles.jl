#### Dataset files

function n_g_points(λ)
  "sw" == λ && return "224"
  "lw" == λ && return "256"
  error("Bad $(λ) in n_g_points")
end

function data_files_dict(datafolder::AbstractString,
                         λ::AbstractString;
                         allsky=false)
  @assert λ=="sw" || λ=="lw"
  file = Dict()
  file[:k_dist] = joinpath(datafolder, "rrtmgp", "data", "rrtmgp-data-$(λ)-g$(n_g_points(λ))-2018-12-04.nc")
  if allsky
    file[:cloud_optics] = joinpath(datafolder, "extensions", "cloud_optics", "rrtmgp-cloud-optics-coeffs-$(λ).nc")
    file[:input] = joinpath(datafolder, "examples", "all-sky", "rrtmgp-allsky.nc")
    file[:ref] = joinpath(datafolder, "examples", "all-sky", "ref", "rrtmgp-allsky.nc")
  else
    clear_sky_dir = joinpath(datafolder, "examples","rfmip-clear-sky")
    file[:rfmip] = joinpath(clear_sky_dir, "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc")
    file[:flx_up] = joinpath(clear_sky_dir, "r$(λ[1])u_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc")
    file[:flx_dn] = joinpath(clear_sky_dir, "r$(λ[1])d_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc")
  end
  return file
end

dataset_dict(files) = Dict([k => Dataset(v, "r") for (k,v) in files])
