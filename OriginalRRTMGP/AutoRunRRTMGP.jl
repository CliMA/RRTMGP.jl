# This file assumes that the directories are organized as follows:
#
# ./SomeDirectory/RRTMGP.jl/
# ./SomeDirectory/rte-rrtmgp/
#

include("ConfigRRTMGP.jl")
include("CleanRRTMGP.jl")
traverse_input_params = false
run_clear_sky = false
run_all_sky = true

cd(joinpath(root_dir,"rte-rrtmgp","build")) do
  run(`make`)
end

if run_clear_sky
  cd(rte_rrtmgp_rfmip_clear_sky) do
    run(`make`)
    if !traverse_input_params
      # stage_files.py calls urllib.request.urlretrieve, so call sparingly (can get blocked out by too many calls)
      run(`python stage_files.py`)
      run(`python run-rfmip-examples.py`)
    else
      println("Running RFMIP drivers")
      conds_file = joinpath(".", "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc")
      lw_coeffs = joinpath(rte_rrtmgp_dir, "rrtmgp", "data", "rrtmgp-data-lw-g256-2018-12-04.nc")
      sw_coeffs = joinpath(rte_rrtmgp_dir, "rrtmgp", "data", "rrtmgp-data-sw-g224-2018-12-04.nc")
      for wavelength in ["l","s"]
        for direction in ["u","d"]
          for forcing_index in ["1",]
          # for forcing_index in ["1","2","3"]
            for physics_index in ["1",]
            # for physics_index in ["1","2"]
              @show wavelength,direction,forcing_index,physics_index
              # Arguments: block size, input conditions, coefficient files, forcing index, physics index
              run(`$(joinpath(".", "rrtmgp_rfmip_lw")) 8 $(conds_file) $(lw_coeffs) $(forcing_index) $(physics_index)`)
              run(`$(joinpath(".", "rrtmgp_rfmip_sw")) 8 $(conds_file) $(sw_coeffs) $(forcing_index) $(physics_index)`)
            end
          end
        end
      end

    end
  end
end

if run_all_sky
  cd(rte_rrtmgp_rfmip_all_sky) do
    # subprocess.run([all_sky_exe_name, atmos_file, lw_gas_coeffs_file, lw_clouds_coeff_file, ncol_str, nloops_str])
    run(`make`)
    run(`python $(joinpath(".", "run-allsky-example.py"))`)
    # run(`python $(joinpath(".", "compare-to-reference.py"))`)
  end
end
