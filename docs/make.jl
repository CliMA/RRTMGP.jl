Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[]) # JuliaLang/julia/pull/28625

using RRTMGP, Documenter

makedocs(
  sitename = "RRTMGP",
  doctest = false,
  strict = false,
  format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  clean = false,
  modules = [Documenter, RRTMGP],
  pages = Any[
    "Home" => "index.md",
    "RRTMGP" => Any[
    "Gas Concentrations" => "RRTMGP/GasConcs.md"
    "Gas Optics" => "RRTMGP/GasOptics.md"
    ],
    "RTE" => Any[
    "Optical Properties" => "RTE/OpticalProps.md"
    "Fluxes" => "RTE/Fluxes.md"
    "RTE Solver" => "RTE/RTESolver.md"
    ],
    "References" => "References.md",
  ],
)

deploydocs(
           repo = "github.com/climate-machine/RRTMGP.jl.git",
           target = "build",
          )
