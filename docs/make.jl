Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[]) # JuliaLang/julia/pull/28625

using JRRTMGP, Documenter

makedocs(
  sitename = "JRRTMGP",
  doctest = false,
  strict = false,
  format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  clean = false,
  modules = [Documenter, JRRTMGP],
  pages = Any[
    "Home" => "index.md",
    "RTE" => "RTE.md",
    "RRTMGP" => Any[ "gas optics rrtmgp" => "RRTMGP/mo_gas_optics_rrtmgp.md" ],
  ],
)

deploydocs(
           repo = "github.com/climate-machine/JRRTMGP.git",
           target = "build",
          )
