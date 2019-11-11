using Test
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# include("AutoRunRRTMGP.jl")
for submodule in [
                  "PartialImplementations", # compile only
                  "DataDeps",
                  "test_rfmip_clear_sky_sw",
                  "test_rfmip_clear_sky_lw",
                  "test_allsky",
                  "Benchmarks",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
