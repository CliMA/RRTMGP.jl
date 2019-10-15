

# PS = Sys.iswindows() ? "\\" : "/"
# PS = "/"
if Sys.iswindows()
  # ENV["NCHOME"] = joinpath("C:","Program Files","netCDF 4.7.1","bin")
  # ENV["NCHOME"] = joinpath("C:","Program Files","netCDF 4.7.1","lib")
  ENV["NCHOME"] = joinpath("C:","Program Files","netCDF 4.7.1")
  # ENV["NFHOME"] = joinpath("C:","MinGW","bin")
  ENV["NFHOME"] = joinpath("C:","MinGW")
else
# ENV["NCHOME"] = "/usr/local/Cellar/netcdf/4.6.1_4"
  ENV["NCHOME"] = "/usr/local/"
  ENV["NFHOME"] = "/usr/local/"
end

cd("../rte-rrtmgp/build") do
  run(`make clean`)
  run(`make`)
end

cd("../rte-rrtmgp/examples/rfmip-clear-sky") do
  run(`make clean`)
  run(`make`)
  run(`python stage_files.py`)
  run(`python run-rfmip-examples.py`)
end

