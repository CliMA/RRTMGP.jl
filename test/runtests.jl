using Test

for submodule in [
                  "adding",
                  "DriverSW",
                  "RTE",
                  "RRTMG",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
