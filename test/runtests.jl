using Test

for submodule in [
                  # "adding",
                  "DriverSW",
                  "DriverLW",
                  "RTE",
                  "RRTMG",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
