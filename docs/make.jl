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
    "RTE" => "RTE.md",
    "RRTMGP" => Any[ "gas optics" => "RRTMGP/mo_gas_optics_rrtmgp.md" ],
  ],
)

deploydocs(
           repo = "github.com/climate-machine/RRTMGP.git",
           target = "build",
          )
