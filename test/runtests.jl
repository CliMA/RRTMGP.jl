using Test
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# Only run the per-grid-point tests
pgp_only = false

# include(joinpath("OriginalRRTMGP","AutoRunRRTMGP.jl"))
for submodule in [
                  "test_allsky",
                  "test_rfmip_clear_sky_lw",
                  "test_rfmip_clear_sky_sw",
                  "solar_zenith_angle",
                  # "Benchmarks",
                  joinpath("..","PartialImplementations", "PartialImplementations"), # compile only
                  "DataDeps",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
