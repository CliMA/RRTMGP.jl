```math
\newcommand{\π}{\pi}                                  % Math
\newcommand{\DEL}{\nabla}                             % Math
\newcommand{\DOT}{\cdot}                              % Math
\newcommand{\DIV}{\DEL \DOT}                          % Math
\newcommand{\PD}{\partial}                            % Math
\newcommand{\dag}{\dagger}                            % Math
```

```math
\newcommand{\dir}{\dag}                               % Direct
\newcommand{\dif}{\star}                              % Diffuse
\newcommand{\SFC}{\mathrm{sfc}}                       % Surface
\newcommand{\TOA}{\mathrm{TOA}}                       % Top of atmosphere
\newcommand{\F}{\mathcal F}                           % Generic fluxes in thermal energy equation
\newcommand{\FR}[1]{F^{#1}_{\mathrm{rad}}}            % Radiative flux
\newcommand{\τ}[1]{\tau_{#1}}                         % Optical depth
\newcommand{\threshtau}[1]{\tau_{#1\mathrm{thresh}}}  % Optical depth
\newcommand{\SSA}{\tilde{\omega}}                     % Single scattering albedo
\newcommand{\SA}{\varepsilon_{\SFC}}                  % Surface albedo
\newcommand{\EMIS}{\epsilon}                          % Emissivity
\newcommand{\λ}{\lambda}                              % Wavelength
\newcommand{\gpt}{\xi}                                % G-point variable
\newcommand{\FS}{\mathcal F_{\mathcal S}}             % Symbol for source function
\newcommand{\OPMAP}[1]{\mathcal M_{#1}^{\lambda}}     % Optical properties map
\newcommand{\SMAP}{\Upsilon^{\lambda}}                % Radiative sources map
\newcommand{\KD}{\mathcal K_{\mathrm{dist}}}          % Symbol for k-distribution
\newcommand{\SV}{\mathcal S}                          % Symbol for thermal energy equation source term
\newcommand{\AS}{\mathcal S_{\mathrm{atmos}}}         % Symbol for atmospheric state
\newcommand{\ASY}{g}                                  % Asymmetry factor
\newcommand{\θ}{\theta}                               % Solar zenith angle
\newcommand{\μ}{\mu}                                  % Cosine of solar zenith angle
\newcommand{\ϕ}{\phi}                                 % Azimuthal angle
\newcommand{\ρ}{\rho}                                 % Density ?
\newcommand{\ν}{\nu}                                  % Frequency
\newcommand{\RAD}[1]{I^{#1}}                          % Radiative intensity
\newcommand{\ks}{k_{\mathrm{scatter}}}                % Absorption coefficient for scattering
\newcommand{\ka}{k_{\mathrm{abs}}}                    % Absorption coefficient
\newcommand{\DA}{D}                                   % Secant of propagation angle
\newcommand{\transmittance}{T}                        % Transmittance
\newcommand{\temperature}{\Theta}                     % Temperature
\newcommand{\pressure}{p}                             % Pressure
\newcommand{\PhaseFun}[2]{P(#1,#2)}                   % Phase function
\newcommand{\SpectralPhaseFun}[2]{P_{\ν}(#1,#2)}      % Spectral Phase function
\newcommand{\ScatPhaseFun}[1]{P(#1)}                  % Scattering phase function
\newcommand{\PlanckF}[1]{\hat{B}(#1)}                 % Planck function
\newcommand{\PlanckS}[1]{B^{#1}}                      % Layer Planck source
\newcommand{\SpecPlanckF}[1]{B_{\ν}(#1)}              % Spectral Planck function
\newcommand{\RSRC}[1]{S^{#1}}                         % Radiation source
\newcommand{\SRSRC}[1]{\underline{S}^{#1}}            % Surface radiation source
\newcommand{\BS}[1]{S^{#1}}                           % Level radiation source
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

Here, ``\F = \sum_i F_i``, with ``F_j = \FR{net} = \FR{+} - \FR{-}`` for some ``j``. Note that ``\FR{+}`` and ``\FR{-}`` denote the sum of diffuse (``\dif``) and direct (``\dir``) upward and downward fluxes respectively.

!!! note
    There are two implementations in `RRTMGP.jl`: one that operates on 1) all layers and multiple columns and 2) a single grid point.

# Map atmospheric state to optical properties

To solve for ``\FR{net}``, one must choose an approximation for optical properties. In `RRTMGP.jl`, choices include

 - [`OneScalar`](@ref Optics.OneScalar): optical depth (``\τ{}``) only
 - [`TwoStream`](@ref Optics.TwoStream): optical depth (``\τ{}``), single scattering albedo (``\SSA``), asymmetry factor (``\ASY``)

To compute either set of optical properties, one must call [`compute_optical_props!`](@ref Optics.compute_optical_props!) with the following arguments:

 - A K-distribution ``\KD``, [`LookUpLW`](@ref LookUpTables.LookUpLW) or [`LookUpSW`](@ref LookUpTables.LookUpSW): typically read from [Data files](@ref)
 - An [`AtmosphericState`](@ref AtmosphericStates.AtmosphericState) ``\AS``, defined by: temperature ``\temperature``, pressure ``\pressure``, set of gas volume mixing ratios ([`AbstractVmr`](@ref Vmrs.AbstractVmr))
 - [optionally] A source function ``\FS``, [`SourceLWNoScat`](@ref Sources.SourceLWNoScat), [`SourceLW2Str`](@ref Sources.SourceLW2Str) or [`SourceSW2Str`](@ref Sources.SourceSW2Str): source functions

Calling `Optics.compute_optical_props!` computes the optical depths via a 3D interpolation routine in the following simplified steps:

 - Initialize and compute [`InterpolationCoefficients`](@ref Optics.compute_interp_fractions)
 - Compute absorption optical depth ``\τ{absorption}`` for major and minor gas species (which has binary variation as a function of height `itropo=1,2`)
 - Compute Rayleigh scattering optical depth ``\τ{Rayleigh}``
 - Compute total optical depth: ``\τ{total} = \τ{absorption}+\τ{Rayleigh}``
 - [optionally] Compute ``\SSA = \τ{Rayleigh}/\τ{total}``
 - [optionally] Compute source (if long-wave radiation)

Mathematically, let's write this entire computation as a map (a spectral map, since it maps wavelength/frequency to g-points):

```math
\begin{align}
\τ{\gpt} = \OPMAP{1}(\KD,\AS[,\FS]) \\
(\τ{\gpt},\SSA_{\gpt},\ASY_{\gpt}) = \OPMAP{2}(\KD,\AS[,\FS]) \\
\end{align}
```

In addition, `Optics.compute_optical_props!` computes radiative sources (currently Planck sources), which are computed the same regardless of the type of optical properties:

```math
\begin{align}
\PlanckS{*}_{\gpt},\PlanckS{+}_{\gpt},\PlanckS{-}_{\gpt} = \SMAP(\KD,\AS,\FS) \\
\end{align}
```

Using the denoting notation:
```math
\begin{align}
* &= \quad \text{at average temperature} \\
+ &= \quad \text{upward} \\
- &= \quad \text{downward}
\end{align}
```


# Computing fluxes given optical depth

!!! note
    This section is under development

The radiative transfer equation under Thermodynamic Equilibrium (TE) in differential form is

```math
\begin{align}
d \frac{d\RAD{}_{\ν}}{\ρ k(\ν) ds} = - \RAD{}_{\ν} + J_{\ν} \\
\end{align}
```
here ``\RAD{}_{\ν}`` is monochromatic intensity, ``J`` is the source function and ``k(\ν) = \ka(\ν)+\ks(\ν)``. Define the single scattering albedo and the normalized phase function as below:
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
    \frac{d\RAD{}_{\nu}(\Omega)}{d\τ{}} = \RAD{}_{\nu}(\Omega) -
    (1-\SSA_{\nu})\SpecPlanckF{T} -
    \frac{\SSA_{\nu}}{4\pi}\int_{4\pi}
    \RAD{}_{\nu}(\Omega^{'}) \SpectralPhaseFun{\Omega}{\Omega^{'}} d\Omega^{'}
\end{align}
```
The terms of right hand side is the attenuation by absorption and scattering, emission, and scattering, respectively.

# Equations derived from code

## Short-wave, without scattering (one-scalar)

```math
\begin{align}
\FR{+} = 0 \\
\FR{-} = \FR{\dir} = \sum_{\gpt} \RAD{\dir}_{\gpt} \\
\RAD{\dir}_{\gpt} = \RAD{\dir,\TOA}_{\gpt} e^{-\τ{} / \μ} \\
\end{align}
```

## Long-wave, without scattering (one-scalar)

```math
\begin{align}
\FR{\pm} = \sum_{\gpt} \RAD{\pm}_{\gpt} \\
\RAD{\pm}_{\gpt} = 2 \π \sum_{\μ} \RAD{\pm}_{\gpt,\μ} w_{\μ} \\
D = \text{quadrature secants} \\
\transmittance = e^{-\τ{} D} \\
f = \begin{cases}
  \frac{1 - \transmittance}{\τ{}} - \transmittance & \τ{} > \threshtau{} \\
  \τ{} \left(\frac{1}{2} - \frac{\τ{}}{3} \right)  & \text{otherwise} \\
\end{cases} \\
\RSRC{\pm}_{\gpt} = (1 - \transmittance) \PlanckS{\pm}_{\gpt} + 2 f (\PlanckS{*}_{\gpt} - \PlanckS{\pm}_{\gpt}) \\
\RAD{+}_{\gpt,\μ} = \transmittance \RAD{\SFC}_{\gpt,\μ} + \RSRC{+}_{\gpt} \\
\RAD{-}_{\gpt,\μ} = \transmittance \RAD{\TOA}_{\gpt,\μ} + \RSRC{-}_{\gpt} \\
\end{align}
```

## Short-wave, with scattering (two-stream)

## Long-wave, with scattering (two-stream)


# Data files

 - [allsky](https://owncloud.gwdg.de/index.php/s/OjbNzRTlXUk0G5w/download) (has both up/down fluxes)
 - [clear-sky long-wave downward flux](https://owncloud.gwdg.de/index.php/s/kbhl3JOSccGtR0m/download)
 - [clear-sky long-wave upward flux](https://owncloud.gwdg.de/index.php/s/5DbhryVSfztioPG/download)
 - [clear-sky short-wave downward flux](https://owncloud.gwdg.de/index.php/s/uCemCHlGxbGK0gJ/download)
 - [clear-sky short-wave upward flux](https://owncloud.gwdg.de/index.php/s/l8ZG28j9ttZWD9r/download)
 - [RFMIP profiles](http://aims3.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/RFMIP/UColorado/UColorado-RFMIP-1-2/atmos/fx/multiple/none/v20190401/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc)

See `data_files_dict` to see how files are saved.
