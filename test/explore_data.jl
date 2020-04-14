using RRTMGP
using NCDatasets
using NCDatasets

include("data_set_files.jl")

datafolder = RRTMGP.data_folder_rrtmgp()

data = Dict()
data["clearsky_sw"] = dataset_dict(data_files_dict(datafolder, "sw"))
data["clearsky_lw"] = dataset_dict(data_files_dict(datafolder, "lw"))
data["allsky_sw"] =
    dataset_dict(data_files_dict(datafolder, "sw"; allsky = true))
data["allsky_lw"] =
    dataset_dict(data_files_dict(datafolder, "lw"; allsky = true))
