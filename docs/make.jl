using Documenter
using LaTeXStrings
using RRTMGP
using DocumenterCitations
import Literate

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

# Generate the tutorial pages from the Literate.jl scripts in docs/tutorials.
# The generated markdown contains `@example` blocks, so Documenter executes
# the tutorials during the docs build (they cannot silently rot).
const tutorials_src = joinpath(@__DIR__, "tutorials")
const tutorials_out = joinpath(@__DIR__, "src", "tutorials")
for script in readdir(tutorials_src)
    endswith(script, ".jl") || continue
    Literate.markdown(
        joinpath(tutorials_src, script),
        tutorials_out;
        documenter = true,
        execute = false,
    )
end

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
        "Tutorials" => Any[
            "A first radiation calculation" => "tutorials/getting_started.md",
            "Radiative-convective equilibrium" => "tutorials/manabe_rce.md",
        ],
        "How-to guides" => Any[
            "Drive RRTMGP from a host model" => "howto/host_model.md",
            "Run on GPUs" => "howto/gpu.md",
            "Cache the lookup tables" => "howto/lookup_cache.md",
            "Get per-band (spectral) fluxes" => "howto/spectral_fluxes.md",
            "Run the validated test problems" => "Example.md",
        ],
        "Explanation" => Any[
            "Functional core" => "functional_core.md",
            "RTE" => "RTE.md",
            "Optics" => "Optics.md",
            "Float32 and Float64" => "precision.md",
            "Fortran and paper concordance" => "concordance.md",
        ],
        "Reference" => Any[
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
