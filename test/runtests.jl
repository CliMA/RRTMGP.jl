using Test

for submodule in [
                  "adding",
                  "RTE",
                  "RRTMG",
                  ]

  println("Testing $submodule")
  include(joinpath(submodule*".jl"))
end
