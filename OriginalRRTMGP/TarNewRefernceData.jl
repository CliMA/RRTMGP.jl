using JRRTMGP

# include("AutoRunRRTMGP.jl")

root_dir = dirname(pwd())
rte_rrtmgp_folder = "rte-rrtmgp"
rte_rrtmgp_dir = joinpath(root_dir,rte_rrtmgp_folder)

data_dir = rte_rrtmgp_dir*"-data"
mkpath(data_dir)
data_folder = rte_rrtmgp_folder*"-data"
data_dir_tar = data_folder*".tar.gz"

for (root, dirs, files) in walkdir(rte_rrtmgp_dir)
  for file in files
    if endswith(file, ".nc")
      new_root = data_dir*replace(root,rte_rrtmgp_dir => "")
      mkpath(new_root)
      cp(joinpath(root, file), joinpath(new_root, file); force=true)
    end
  end
end

cd(data_dir) do
  rm(data_dir_tar; force=true)
  JRRTMGP.compress_targz(data_dir_tar)
end
