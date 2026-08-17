# Float32 and Float64

Every RRTMGP.jl kernel is generic in the working precision: the float type
chosen at construction (through [`RRTMGPGridParams`](@ref
RRTMGP.RRTMGPGridParams), [`solve`](@ref RRTMGP.solve), or [`solve_gray`](@ref
RRTMGP.solve_gray)) flows through the lookup tables, the optics, and the RTE
solve. `Float64` is the reference; `Float32` halves the memory traffic of the
bandwidth-bound kernels and is the precision at which the CliMA atmosphere runs
radiation on GPUs.

## Accuracy achieved at Float32

A consistency test (`test/float32_consistency.jl`) solves the same gray,
clear-sky, and cloudy problems at both precisions, with the `Float32` states
rounded from the `Float64` ones, and measures the end-to-end difference: input
rounding, gas and cloud optics, and the RTE solve together. The measured maximum
broadband flux differences are:

- every longwave path is at its interpolation-noise floor of about 2–5 × 10⁻⁴
  W/m², 25–90 times smaller than before the single-precision measures described
  below;
- shortwave differences are about 1–2 × 10⁻² W/m² in clear skies and about 5 ×
  10⁻² W/m² with clouds. These were verified to be per-g-point interpolation and
  coefficient noise rather than accumulation error: repeating the runs with full
  `Float64` broadband accumulation left them unchanged, so compensated summation
  was deliberately not added.

For scale, these differences are orders of magnitude below both typical flux
magnitudes (hundreds of W/m²) and the accuracy of the correlated-``k``
approximation itself; the Fortran reference (rte-rrtmgp) accepts an absolute
flux error of 0.35 W/m² in its own validation. CI enforces ratcheting thresholds
(gray 10⁻³ W/m² and 10⁻⁸ K/s for the heating rate; longwave 10⁻³ W/m²; shortwave
3 × 10⁻² clear and 1.2 × 10⁻¹ W/m² cloudy) that are tightened as numerics fixes
land and never loosened. The comparisons against the rte-rrtmgp reference
results also run at both precisions.

## How: the single-precision numerics

The `Float32` accuracy rests on evaluating exact reformulations instead of
cancellation-prone ones, and on switching to series expansions where
cancellation is unavoidable. The measures, roughly in the order radiation flows:

- **Gas optics.** The interpolation fraction in the binary species parameter
  ``η`` is computed relative to its clamped table cell, so it stays continuous
  at ``η = 1``, a value `Float32` rounds to more often than `Float64`.
- **Longwave, no scattering.** The linear-in-``τ`` source factor
  ``(1 - t)/\tau - t`` ([clough1992](@citet), Eq. 13) loses
  ``\sim\!\mathrm{eps}/\tau^2`` to cancellation for thin layers; it switches to
  its third-order series at ``\tau = \mathrm{eps}^{1/4}``, where the two error
  curves cross.
- **Two-stream coefficients.** The diffusion eigenvalue
  ``k = \sqrt{(\gamma_1-\gamma_2)(\gamma_1+\gamma_2)}`` uses exact
  ``\gamma_1 - \gamma_2`` identities rather than subtracting near-equal numbers,
  with ``k^2`` floored at ``\sqrt{\mathrm{eps}}``; thin-layer factors
  ``1 - e^{-k\tau}`` are evaluated with `expm1`.
- **Longwave two-stream source.** The [toon1989](@citet) source terms are
  factored so that no term divides by ``\tau`` and no small difference of
  near-equal quantities appears (exact
  ``1 \mp R_{\mathrm{dif}} - T_{\mathrm{dif}}`` factorizations).
- **Shortwave direct beam.** The direct reflectance and transmittance of
  [meador1980](@citet) (Eqs. 14–15) have a removable singularity at
  ``k\mu_0 = 1``; ``k\mu_0`` is nudged off resonance symmetrically within a
  ``\sqrt{\mathrm{eps}}`` window, keeping numerator and denominator mutually
  consistent. An energy-conservation clamp then bounds the direct reflectance
  and transmittance by the energy the unscattered beam leaves behind, following
  [ukkonen2024](@citet); and the direct-beam profile itself is built by
  accumulating optical depth downward rather than by repeated multiplication.
- **Delta scaling.** The ``f = g^2`` similarity scaling of cloud and aerosol
  optics ([joseph1976](@citet)) is evaluated in exact, non-cancelling forms (for
  example ``g' = g/(1+g)`` instead of ``(g - g^2)/(1 - g^2)``).

All the guard constants involved (`k_min`, `τ_thresh`, `resonance_window`,
`μ₀_min`) live in the `RRTMGP.Numerics` module, each with the derivation of its
value, and all scale with `eps(FT)`, so the kernels behave consistently at
either precision, and the same policy would extend to other float types. Where a
threshold matches the Fortran reference (rte-rrtmgp), the docstring says so.
