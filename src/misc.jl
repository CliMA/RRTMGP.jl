using DataDeps

"""
    extract_targz(file)

Platform-independent file extraction
"""
function extract_targz(file)
  if Sys.iswindows()
    run(pipeline(`7z x -tgzip -so $file`, `7z x -si -ttar`))
  else
    run(`tar -xzf $file`)
  end
end

function data_folder_rrtmgp()
  register(DataDep("rte-rrtmgp",
                   "data for rte-rrtmgp repo",
                   "https://caltech.box.com/shared/static/8x6dsqt9puv8cxsg875dkp45f3g5dvug.gz",
                   "d3d25ea7c9382fe2016257fe5cb96565da8dc2faa3cd6e3aab6a5503176886dc",
                   post_fetch_method=extract_targz))
  datafolder = datadep"rte-rrtmgp"
  return datafolder
end
function data_folder_rrtmgp_test_val()
  register(DataDep("rte-rrtmgp-test-val",
                   "data for rte-rrtmgp-test-val repo",
                   "https://caltech.box.com/shared/static/j1mksrrq0bfhsdaeupi9vri2jtjfepfv.gz",
                   "5f43b57fd021703a557e98161944b012dcaded22a48bc3437848a43219fb1b72",
                   post_fetch_method=extract_targz))
  datafolder = datadep"rte-rrtmgp-test-val"
  return datafolder
end