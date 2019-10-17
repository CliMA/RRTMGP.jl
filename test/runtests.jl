using Test

for submodule in [
                  # "adding",
                  "DriverSW",
                  "test_two_stream",
                  "RTE",
                  "RRTMG",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
