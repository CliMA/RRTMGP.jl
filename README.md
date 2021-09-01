# RRTMGP.jl

Julia implementation of Rapid and accurate Radiative Transfer Model for General Circulation Models.

|||
|---------------------:|:----------------------------------------------|
| **Documentation**    | [![latest][docs-latest-img]][docs-latest-url] |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **BuildKite**        | [![buildkite][buildkite-img]][buildkite-url]  |
| **Bors**             | [![Bors enabled][bors-img]][bors-url]         |


[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://CliMA.github.io/RRTMGP.jl/latest/

[codecov-img]: https://codecov.io/gh/CliMA/RRTMGP.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/RRTMGP.jl

[buildkite-img]: https://badge.buildkite.com/ee3a0c43cf4925ee14a966f794ac85d0b9439244d23e43b308.svg
[buildkite-url]: https://buildkite.com/clima/rrtmgp-ci

[bors-img]: https://bors.tech/images/badge_small.svg
[bors-url]: https://app.bors.tech/repositories/24778


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
 - NASA for images of the sun (for our [logo](https://clima.github.io/RRTMGP.jl/latest/assets/logo.png))
