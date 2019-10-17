ENV["DATADEPS_ALWAYS_ACCEPT"] = true

using DataDeps
using JRRTMGP
register(DataDep("rte-rrtmgp-test-data",
                 "Reference files for rte-rrtmgp tests",
                 "https://caltech.box.com/shared/static/j1mksrrq0bfhsdaeupi9vri2jtjfepfv.gz",
                 "5f43b57fd021703a557e98161944b012dcaded22a48bc3437848a43219fb1b72",
                 post_fetch_method=JRRTMGP.extract_targz))

testdata = datadep"rte-rrtmgp-test-data"

