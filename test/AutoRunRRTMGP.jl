

if Sys.iswindows()
  ENV["NCHOME"] = "/usr/local/"
  ENV["NFHOME"] = "/usr/local/"
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

