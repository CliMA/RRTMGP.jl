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
    "Mathematical Formulation" => "MathFormulation.md",
    "RRTMGP" => Any[
    "Gas Concentrations" => "RRTMGP/GasConcs.md"
    "Gas Optics" => "RRTMGP/GasOptics.md"
    "References States" => "RRTMGP/ReferenceStates.md"
    "Atmospheric States" => "RRTMGP/AtmosphericStates.md"
    "Physical Constants" => "RRTMGP/PhysicalConstants.md"
    "Optical Properties" => "OpticalProps.md"
    ],
    "RTE" => Any[
    "Optical Properties" => "OpticalProps.md"
    "Fluxes" => "RTE/Fluxes.md"
    "RTE Solver" => "RTE/RTESolver.md"
    "Solar Zenith Angle" => "RTE/SolarZenithAngle.md"
    ],
    "References" => "References.md",
  ],
)

deploydocs(
           repo = "github.com/climate-machine/RRTMGP.jl.git",
           target = "build",
          )
