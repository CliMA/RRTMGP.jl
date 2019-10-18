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
                   "https://caltech.box.com/shared/static/7uzxzvdvl0hi91eyd1unws8qlojueqq7.gz",
                   "300a9f398c3239ff8d16adba50d6f8ef87ebf5e32c9308ab6aa3887ba603c55e",
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