using Test

for submodule in [
                  "RTE",
                  "RRTMG",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
