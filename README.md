# RRTMGP.jl

Julia implementation of Rapid and accurate Radiative Transfer Model for General Circulation Models.

|||
|---------------------:|:----------------------------------------------|
| **Documentation**    | [![latest][docs-latest-img]][docs-latest-url] |
| **Azure Build**      | [![azure][azure-img]][azure-url]              |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **Bors**             | [![Bors enabled][bors-img]][bors-url]         |
| **Travis Build**     | [![travis][travis-img]][travis-url]           |
| **AppVeyor Build**   | [![appveyor][appveyor-img]][appveyor-url]     |

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://CliMA.github.io/RRTMGP.jl/latest/

[azure-img]: https://dev.azure.com/climate-machine/RRTMGP.jl/_apis/build/status/climate-machine.RRTMGP.jl?branchName=master
[azure-url]: https://dev.azure.com/climate-machine/RRTMGP.jl/_build/latest?definitionId=1&branchName=master

[codecov-img]: https://codecov.io/gh/CliMA/RRTMGP.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/RRTMGP.jl

[bors-img]: https://bors.tech/images/badge_small.svg
[bors-url]: https://app.bors.tech/repositories/24778

[travis-img]: https://travis-ci.org/CliMA/RRTMGP.jl.svg?branch=master
[travis-url]: https://travis-ci.org/CliMA/RRTMGP.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/c6eykd0w94pmyjt8/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/climate-machine/rrtmgp-jl/branch/master

This is based off of the [rte-rrtmgp](https://github.com/RobertPincus/rte-rrtmgp) repository.

# Install

RRTMGP.jl is registered in the general Julia registry. To install, enter the package manager by typing `]` in the Julia REPL, and then type:

```julia
pkg> add RRTMGP
```

Then, to use

```julia
julia> using RRTMGP
```

# Acknowledgments

 - [Robert Pincus](https://github.com/RobertPincus) for his invaluable help (and developing RRTMGP in the first place).
 - NASA for images of the sun (for our [logo](https://climate-machine.github.io/RRTMGP.jl/latest/assets/logo.png))
