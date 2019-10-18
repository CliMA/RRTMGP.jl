using Test

for submodule in [
                  "DataDeps",
                  "DriverSW",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
