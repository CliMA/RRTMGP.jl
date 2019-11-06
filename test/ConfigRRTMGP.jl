# This file assumes that the directories are organized as follows:
#
# ./SomeDirectory/JRRTMGP
# ./SomeDirectory/rte-rrtmgp
#

if Sys.iswindows()
  # ENV["NCHOME"] = joinpath("C:","Program Files","netCDF 4.7.1","bin")
  # ENV["NCHOME"] = joinpath("C:","Program Files","netCDF 4.7.1","lib")
  ENV["NCHOME"] = joinpath("C:","Program Files","netCDF 4.7.1")
  # ENV["NFHOME"] = joinpath("C:","MinGW","bin")
  ENV["NFHOME"] = joinpath("C:","MinGW")
else
# ENV["NCHOME"] = "/usr/local/Cellar/netcdf/4.6.1_4"
  ENV["NCHOME"] = "/usr/local"
  ENV["NFHOME"] = "/usr/local"
  # ENV["NCHOME"] = joinpath("/usr","local","Cellar","netcdf","4.6.1_4")
  # ENV["NFHOME"] = joinpath("/usr","local","Cellar","netcdf","4.6.1_4")
end

root_dir = dirname(pwd())
rte_rrtmgp_dir = joinpath(root_dir,"rte-rrtmgp")
rte_rrtmgp_rfmip_clear_sky = joinpath(root_dir,"rte-rrtmgp","examples","rfmip-clear-sky")
rte_rrtmgp_rfmip_all_sky = joinpath(root_dir,"rte-rrtmgp","examples","all-sky")
rte_rrtmgp_test_val_dir = joinpath(root_dir,"rte-rrtmgp-test-val")

# rte-rrtmgp ENV variables
ENV["RTE_KERNELS"] = joinpath(rte_rrtmgp_dir,"rte","kernels")

# rte-rrtmgp and rte-rrtmgp-test-val ENV variables
ENV["RTE_DIR"] = joinpath(rte_rrtmgp_dir,"rte")
ENV["RTE_KERNEL_DIR"] = joinpath(rte_rrtmgp_dir,"rte","kernels")
ENV["RRTMGP_KERNEL_DIR"] = joinpath(rte_rrtmgp_dir,"rrtmgp","kernels")
ENV["RRTMGP_DIR"] = joinpath(rte_rrtmgp_dir,"build")
ENV["FCFLAGS"]="-fimplicit-none -ffree-line-length-none -Wuninitialized -cpp -std=gnu"
ENV["FC"] = "gfortran"

# rte-rrtmgp-test-val ENV variables
ENV["TEST_ROOT"] = rte_rrtmgp_test_val_dir
ENV["RRTMGP_ROOT"] = rte_rrtmgp_dir
ENV["RRTMGP_BUILD"] = joinpath(rte_rrtmgp_dir,"build")
