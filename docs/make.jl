using Documenter
using LaTeXStrings
using RRTMGP
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

makedocs(;
    plugins = [bib],
    sitename = "RRTMGP.jl",
    doctest = false,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(
            Dict(
                :TeX => Dict(
                    :equationNumbers => Dict(:autoNumber => "AMS"),
                    :Macros => Dict(),
                ),
            ),
        ),
    ),
    clean = true,
    modules = [RRTMGP],
    pages = Any[
        "Home" => "index.md",
        "Functional core" => "functional_core.md",
        "RTE" => "RTE.md",
        "Optics" => "Optics.md",
        "Example" => "Example.md",
        "Fortran and paper concordance" => "concordance.md",
        "API" => Any[
            "API" => "api.md"
            "Getter contract" => "getters.md"
            "Angular Discretization" => "optics/AngularDiscretizations.md"
            "Artifacts" => "ArtifactPaths.md"
            "Atmospheric State" => "optics/AtmosphericStates.md"
            "Boundary Conditions" => "optics/BCs.md"
            "Fluxes" => "optics/Fluxes.md"
            "Utilities for gray radiation simulation" => "optics/GrayAtmosphere.md"
            "Lookup Tables" => "optics/LookUpTables.md"
            "Optics" => "optics/Optics.md"
            "Solver struct" => "optics/RTE.md"
            "Source Functions" => "optics/Sources.md"
            "Volume Mixing Ratios" => "optics/VolumeMixingRatios.md"
            "RTE Solver" => "rte/RTESolver.md"
        ],
        "References" => "References.md",
    ],
)

deploydocs(
    repo = "github.com/CliMA/RRTMGP.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
