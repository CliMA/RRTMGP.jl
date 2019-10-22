using Test

# include("AutoRunRRTMGP.jl")
for submodule in [
                  "DataDeps",
                  "DriverSW",
                  "DriverLW",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
