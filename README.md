# RRTMGP.jl

Julia implementation of Radiative Transfer for Energetics and the RRTMGP optics.

|||
|---------------------:|:----------------------------------------------|
| **Documentation**    | [![latest][docs-latest-img]][docs-latest-url] |
| **Code Coverage**    | [![codecov][codecov-img]][codecov-url]        |
| **BuildKite**        | [![buildkite][buildkite-img]][buildkite-url]  |
| **Downloads**        | [![downloads][downloads-img]][downloads-url]  |


[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://CliMA.github.io/RRTMGP.jl/latest/

[codecov-img]: https://codecov.io/gh/CliMA/RRTMGP.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/RRTMGP.jl

[buildkite-img]: https://badge.buildkite.com/ee3a0c43cf4925ee14a966f794ac85d0b9439244d23e43b308.svg
[buildkite-url]: https://buildkite.com/clima/rrtmgp-ci

[downloads-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FRRTMGP&query=total_requests&suffix=%2Ftotal&label=Downloads
[downloads-url]: http://juliapkgstats.com/pkg/RRTMGP

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

 - [Robert Pincus](https://github.com/RobertPincus) for his invaluable help. 
 - The authors of the [Fortran implementation](https://github.com/earth-system-radiation/rte-rrtmgp) of RTE-RRTMGP on which this code is based
 - NASA for images of the sun (for our [logo](https://clima.github.io/RRTMGP.jl/latest/assets/logo.png))
