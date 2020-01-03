```math
\newcommand{\DEL}{\nabla}
\newcommand{\DOT}{\cdot}
\newcommand{\DIV}{\DEL \DOT}
\newcommand{\PD}{\partial}
\newcommand{\F}{\mathcal F}
\newcommand{\sum}{\Sigma}
\newcommand{\FR}{F_{\mathrm{rad}}}
\newcommand{\SV}{\mathcal S}
\newcommand{\SMAP}{\mathcal M^{\lambda}}
\newcommand{\OD}[1]{\tau_{#1}}
\newcommand{\SSA}{\bar{\omega}}
\newcommand{\FS}{\mathcal F_{\mathcal S}}
\newcommand{\KD}{\mathcal K_{\mathrm{dist}}}
\newcommand{\AS}{\mathcal S_{\mathrm{atmos}}}
\newcommand{\SOD}{\OD{\lambda}}
\newcommand{\Sg}{g_{\lambda}}
\newcommand{\SSSA}{\SSA_{\lambda}}
```

# Radiation in the context of the thermal energy equation

Radiative fluxes enter the thermal energy equation as follows:

```math
\begin{align}
\PD_t (\rho e) = - \DIV \F + \SV
\end{align}
```

Here, $\F = \sum_i F_i$, with $F_j = \FR$ for some $j$.

!!! note
    There are two implementations in `RRTMGP.jl`: one that operates on 1) all layers and multiple columns and 2) a single grid point.

# Map atmospheric state to optical properties

To solve for $\FR$, one must choose an approximation for optical properties. In `RRTMGP.jl`, choices include

 - [`OneScalar`](@ref): optical depth ($\OD{}$) only
 - [`TwoStream`](@ref): optical depth ($\OD{}$), single scattering albedo ($\SSA$), asymmetry factor (``g``)

To compute either set of optical properties, one must call [`gas_optics!`](@ref) with the following arguments:

 - A [K-Distribution](@ref) $\KD$, `KDistributionLongwave` or `KDistributionShortwave`: typically read from [Data files](@ref)
 - An [`AtmosphericState`](@ref) $\AS$, defined by: temperature $T$, pressure $p$, set of gas concentrations ([`AbstractGas`](@ref))
 - [optionally] An [`AbstractSourceFunc`](@ref) $\FS$, `SourceFuncLongWave` or `SourceFuncShortWave`: source functions

Calling `gas_optics!` computes the optical depths via a 3D interpolation routine in the following simplified steps:

 - Initialize and compute [`InterpolationCoefficients`](@ref)
 - Compute absorption optical depth $\OD{absorption}$ for major and minor gas species (which has binary variation as a function of height `itropo=1,2`)
 - Compute Rayleigh scattering optical depth $\OD{Rayleigh}$
 - Compute total optical depth: $\OD{total} = \OD{absorption}+\OD{Rayleigh}$
 - [optionally] Compute $\SSA = \OD{Rayleigh}/\OD{total}$
 - [optionally] Compute source (if long-wave radiation)

Mathematically, let's write this entire computation as a map (a spectral map, since it depends on wavelength):

```math
\begin{align}
\OD{\lambda} = \SMAP_1(\KD,\AS[,\FS]) \\
(\SOD,\SSSA,\Sg) = \SMAP_2(\KD,\AS[,\FS]) \\
\end{align}
```

# Computing fluxes given optical depth

## One-Scalar case


## Two-Stream case


# Data files

 - [allsky](https://owncloud.gwdg.de/index.php/s/OjbNzRTlXUk0G5w/download) (has both up/down fluxes)
 - [clear-sky long-wave downward flux](https://owncloud.gwdg.de/index.php/s/kbhl3JOSccGtR0m/download)
 - [clear-sky long-wave upward flux](https://owncloud.gwdg.de/index.php/s/5DbhryVSfztioPG/download)
 - [clear-sky short-wave downward flux](https://owncloud.gwdg.de/index.php/s/uCemCHlGxbGK0gJ/download)
 - [clear-sky short-wave upward flux](https://owncloud.gwdg.de/index.php/s/l8ZG28j9ttZWD9r/download)
 - [RFMIP profiles](`http://aims3.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/atmos/fx/multiple/none/v20190401/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc)

See `data_files_dict` to see how files are saved.
