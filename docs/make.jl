using RRTMGP, Documenter, LaTeXStrings

makedocs(
    sitename = "RRTMGP.jl",
    doctest = false,
    strict = true,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict(),
            ),
        )),
    ),
    clean = true,
    modules = [RRTMGP],
    pages = Any[
        "Home"=>"index.md",
        "Mathematical Formulation"=>"MathFormulation.md",
        "RRTMGP"=>Any[
            "Gas Concentrations" => "RRTMGP/GasConcs.md"
            "Gases" => "RRTMGP/Gases.md"
            "Gas Optics" => "RRTMGP/GasOptics.md"
            "References States" => "RRTMGP/ReferenceStates.md"
            "Source Functions" => "RRTMGP/SourceFunctions.md"
            "K-Distribution" => "RRTMGP/KDistribution.md"
            "Atmospheric States" => "RRTMGP/AtmosphericStates.md"
            "Optical Properties" => "OpticalProps.md"
        ],
        "RTE"=>Any[
            "Optical Properties" => "OpticalProps.md"
            "Fluxes" => "RTE/Fluxes.md"
            "Boundary Conditions" => "RTE/BoundaryConditions.md"
            "Angular Discretizations" => "RTE/AngularDiscretizations.md"
            "RTE Solver" => "RTE/RTESolver.md"
        ],
        "References"=>"References.md",
    ],
)

deploydocs(repo = "github.com/CliMA/RRTMGP.jl.git", target = "build")
