# This file assumes that the directories are organized as follows:
#
# ./SomeDirectory/RRTMGP.jl/
# ./SomeDirectory/rte-rrtmgp/
#

include("ConfigRRTMGP.jl")

cd(joinpath(root_dir,"rte-rrtmgp","build")) do
  rm("Makefile.libs"; force=true)
  rm("Makefile.conf"; force=true)
  rm("Makefile.rules"; force=true)
end

cd(joinpath(root_dir,"rte-rrtmgp","build")) do
  run(`make clean`)
end

cd(joinpath(root_dir,"rte-rrtmgp")) do
  rm(".DS_store"; force=true)
end
cd(joinpath(root_dir,"rte-rrtmgp","examples")) do
  rm(".DS_store"; force=true)
end

cd(rte_rrtmgp_rfmip_clear_sky) do
  conds_file = joinpath(".", "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc")
  rm(conds_file; force=true)
  rm(".DS_store"; force=true)
  rm("generate-output-file-templates.py"; force=true)
  rm("rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"; force=true)
  rm("rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"; force=true)
  rm("rrtmgp_rfmip_lw"; force=true)
  rm("rrtmgp_rfmip_sw"; force=true)
  rm("rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"; force=true)
  rm("rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc"; force=true)
  run(`make clean`)
end

cd(rte_rrtmgp_rfmip_all_sky) do
  run(`make clean`)
end
