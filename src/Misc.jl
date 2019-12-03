using DataDeps

using Pkg
export haspkg
"""
    haspkg(pkgname::String)::Bool

Determines if the package `pkgname` available in the current environment.
"""
haspkg(pkgname::String) = haskey(Pkg.installed(), pkgname)

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
                   # "https://caltech.box.com/shared/static/8x6dsqt9puv8cxsg875dkp45f3g5dvug.gz", # commit ed5b0113109fcd23a010a90c61f21bad551146ef (original)
                   # "https://caltech.box.com/shared/static/j64qnw96izvheukk3avvz1j3js2g9kpw.gz", # commit ce295e55713550723797c128fadad7964576617d (cloud_optics added)
                   "https://caltech.box.com/shared/static/urnr3dlv0ff5uwt94t9g8zoon9t5jroa.gz", # commit ce295e55713550723797c128fadad7964576617d (cloud_optics added, with all-sky reference data)
                   "a30ef1e22aeb1245d838b870449cebc4ea52c3c883f8114f124b0d06457cacc0",
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