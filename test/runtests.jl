using Test
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# include(joinpath("OriginalRRTMGP","AutoRunRRTMGP.jl"))
for submodule in [
                  joinpath("..","PartialImplementations", "PartialImplementations"), # compile only
                  "DataDeps",
                  "test_allsky",
                  "test_rfmip_clear_sky_lw",
                  "test_rfmip_clear_sky_sw",
                  "solar_zenith_angle",
                  # "Benchmarks",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
