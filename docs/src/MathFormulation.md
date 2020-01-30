```math
\newcommand{\π}{\pi}                                  % Math
\newcommand{\DEL}{\nabla}                             % Math
\newcommand{\DOT}{\cdot}                              % Math
\newcommand{\DIV}{\DEL \DOT}                          % Math
\newcommand{\PD}{\partial}                            % Math
\newcommand{\sum}{\Sigma}                             % Math
```

```math
\newcommand{\F}{\mathcal F}                           % Generic fluxes in thermal energy equation
\newcommand{\FR}{F_{\mathrm{rad}}}                    % Radiative flux
\newcommand{\τ}[1]{\tau_{#1}}                         % Optical depth
\newcommand{\SSA}{\tilde{\omega}}                     % Single scattering albedo
\newcommand{\λ}{\lambda}                              % Wavelength
\newcommand{\gpt}{\xi}                                % G-point variable
\newcommand{\FS}{\mathcal F_{\mathcal S}}             % Symbol for source function
\newcommand{\SMAP}{\mathcal M^{\lambda}}              % Symbol for spectral map
\newcommand{\KD}{\mathcal K_{\mathrm{dist}}}          % Symbol for k-distribution
\newcommand{\SV}{\mathcal S}                          % Symbol for thermal energy equation source term
\newcommand{\AS}{\mathcal S_{\mathrm{atmos}}}         % Symbol for atmospheric state
\newcommand{\ASY}{g}                                  % Asymmetry factor
\newcommand{\θ}{\theta}                               % Solar zenith angle
\newcommand{\μ}{\mu}                                  % Cosine of solar zenith angle
\newcommand{\ϕ}{\phi}                                 % Azimuthal angle
\newcommand{\ρ}{\rho}                                 % Density ?
\newcommand{\ν}{\nu}                                  % Frequency
\newcommand{\RAD}{I}                                  % Radiative intensity
\newcommand{\ks}{k_{\mathrm{scatter}}}                % Absorption coefficient for scattering
\newcommand{\ka}{k_{\mathrm{abs}}}                    % Absorption coefficient
\newcommand{\DA}{D}                                   % Secant of propagation angle
\newcommand{\transmittance}{T}                        % Transmittance
\newcommand{\temperature}{\Theta}                     % Temperature
\newcommand{\pressure}{p}                             % Pressure
\newcommand{\PhaseFun}[2]{P(#1,#2)}                   % Phase function
\newcommand{\SpectralPhaseFun}[2]{P_{\ν}(#1,#2)}      % Spectral Phase function
\newcommand{\ScatPhaseFun}[1]{P(#1)}                  % Scattering phase function
\newcommand{\PlanckF}[1]{B(#1)}                       % Planck function
\newcommand{\SpecPlanckF}[1]{B_{\ν}(#1)}                 % Spectral Planck function
```

```@meta
CurrentModule = RRTMGP
```

# Radiation in the context of the thermal energy equation

Radiative fluxes enter the thermal energy equation as follows:

```math
\begin{align}
\PD_t (\rho e) = - \DIV \F + \SV
\end{align}
```

Here, ``\F = \sum_i F_i``, with ``F_j = \FR`` for some ``j``.

!!! note
    There are two implementations in `RRTMGP.jl`: one that operates on 1) all layers and multiple columns and 2) a single grid point.

# Map atmospheric state to optical properties

To solve for ``\FR``, one must choose an approximation for optical properties. In `RRTMGP.jl`, choices include

 - [`OneScalar`](@ref OpticalProps.OneScalar): optical depth (``\τ{}``) only
 - [`TwoStream`](@ref OpticalProps.TwoStream): optical depth (``\τ{}``), single scattering albedo (``\SSA``), asymmetry factor (``\ASY``)

To compute either set of optical properties, one must call [`gas_optics!`](@ref GasOptics.gas_optics!) with the following arguments:

 - A [K-Distribution](@ref) ``\KD``, `KDistributionLongwave` or `KDistributionShortwave`: typically read from [Data files](@ref)
 - An [`AtmosphericState`](@ref AtmosphericStates.AtmosphericState) ``\AS``, defined by: temperature ``\temperature``, pressure ``\pressure``, set of gas concentrations ([`AbstractGas`](@ref Gases.AbstractGas))
 - [optionally] An [`AbstractSourceFunc`](@ref SourceFunctions.AbstractSourceFunc) ``\FS``, `SourceFuncLongWave` or `SourceFuncShortWave`: source functions

Calling `GasOptics.gas_optics!` computes the optical depths via a 3D interpolation routine in the following simplified steps:

 - Initialize and compute [`InterpolationCoefficients`](@ref GasOptics.InterpolationCoefficients)
 - Compute absorption optical depth ``\τ{absorption}`` for major and minor gas species (which has binary variation as a function of height `itropo=1,2`)
 - Compute Rayleigh scattering optical depth ``\τ{Rayleigh}``
 - Compute total optical depth: ``\τ{total} = \τ{absorption}+\τ{Rayleigh}``
 - [optionally] Compute ``\SSA = \τ{Rayleigh}/\τ{total}``
 - [optionally] Compute source (if long-wave radiation)

Mathematically, let's write this entire computation as a map (a spectral map, since it maps wavelength/frequency to g-points):

```math
\begin{align}
\τ{\gpt} = \SMAP_1(\KD,\AS[,\FS]) \\
(\τ{\gpt},\SSA_{\gpt},\ASY_{\gpt}) = \SMAP_2(\KD,\AS[,\FS]) \\
\end{align}
```

# Computing fluxes given optical depth

The radiative transfer equation under Thermodynamic Equilibrium (TE) in differential form is

```math
\begin{align}
d \frac{d\RAD_{\ν}}{\ρ k(\ν) ds} = - \RAD_{\ν} + J_{\ν} \\
\end{align}
```
here ``\RAD_{\ν}`` is monochromatic intensity, ``J`` is the source function and ``k(\ν) = \ka(\ν)+\ks(\ν)``. Define the single scattering albedo and the normalized phase function as below:
```math
\begin{align}
\SSA(\nu) = \frac{\ks(\nu)}{k(\nu)}\\
\frac{1}{4\pi}\int \SpectralPhaseFun{\Omega}{\Omega^{'}} d\Omega^{'}=1
\end{align}
```

The phase function tells the fraction of radiation scattering by an individual particle from ``\Omega^{'}`` to ``\Omega``, the direction of ``I_{\nu}``. Then ``J`` can be written as
```math
\begin{align}
    J_{\nu} = (1-\SSA_{\nu})\SpecPlanckF{T} + \frac{\SSA_{\nu}}{4\pi}\int_{4\pi} I_{\nu}(\Omega^{'})\SpectralPhaseFun{\Omega}{\Omega^{'}}d\Omega^{'}
\end{align}
```
If we define the optical depth ``d\τ{}`` as opposite to ``ds``, i.e., ``d\τ{} = -\rho k(\nu)ds``. Therefore, the radiative transfer equation is equivalent to
```math
\begin{align}
    \frac{d\RAD_{\nu}(\Omega)}{d\τ{}} = \RAD_{\nu}(\Omega) -
    (1-\SSA_{\nu})\SpecPlanckF{T} -
    \frac{\SSA_{\nu}}{4\pi}\int_{4\pi}
    \RAD_{\nu}(\Omega^{'}) \SpectralPhaseFun{\Omega}{\Omega^{'}} d\Omega^{'}
\end{align}
```
The terms of right hand side is the attenuation by absorption and scattering, emission, and scattering, respectively.

## Compactly

The general form of the radiative transfer equation is

```math
\begin{align}
\μ \frac{d \RAD(\τ{},\μ)}{d\τ{}} = \RAD(\τ{},\μ) - (1-\SSA) \PlanckF{T} - \frac{\SSA}{2} \int_{-1}^{1} \RAD(\τ{},\μ) P d\μ' \\
\μ = \cos(\θ) \\
\PhaseFun{\μ}{\μ'} = \frac{1}{2\π} \int_{2\π}^{0} \ScatPhaseFun{\cos(\Theta)} d\ϕ \\
\cos(\Theta) = \μ\μ' + (1-{\μ'}^2)^{1/2}(1-\μ^2)^{1/2}\cos(\ϕ) \\
\end{align}
```

## Short-wave, no scattering (one-scalar)

```math
\begin{align}
\FR = \sum \RAD_{\gpt} \\
\RAD_{\gpt} = \RAD_{\gpt}^{TOA} e^{-\τ{} / \μ} \\
\end{align}
```

## Long-wave, no scattering (one-scalar)

```math
\begin{align}
\FR = \sum \RAD_{\gpt} \\
\RAD_{\gpt}^{\pm} = \pm 2 \π \sum \RAD_{\gpt,\μ} w_{\μ} \\
\transmittance = e^{-\τ{} D}
\RAD_{\gpt,\μ} = (1 - \transmittance)
\left\{
  B_U + 2 (\bar{B} - B_U)
\right\}
\end{align}
```


## Two-Stream case


# Data files

 - [allsky](https://owncloud.gwdg.de/index.php/s/OjbNzRTlXUk0G5w/download) (has both up/down fluxes)
 - [clear-sky long-wave downward flux](https://owncloud.gwdg.de/index.php/s/kbhl3JOSccGtR0m/download)
 - [clear-sky long-wave upward flux](https://owncloud.gwdg.de/index.php/s/5DbhryVSfztioPG/download)
 - [clear-sky short-wave downward flux](https://owncloud.gwdg.de/index.php/s/uCemCHlGxbGK0gJ/download)
 - [clear-sky short-wave upward flux](https://owncloud.gwdg.de/index.php/s/l8ZG28j9ttZWD9r/download)
 - [RFMIP profiles](http://aims3.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/atmos/fx/multiple/none/v20190401/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc)

See `data_files_dict` to see how files are saved.
