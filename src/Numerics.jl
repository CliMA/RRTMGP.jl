"""
    Numerics

The numerical policy of RRTMGP.jl in one place: every floating-point guard
constant used by the optics/RTE kernels, with its derivation, plus small
numerical utilities. All constants scale with the working precision `FT`
(compare `epsilon(1._wp)`-based parameters in the Fortran reference,
rte-rrtmgp), so kernels behave consistently at Float32 and Float64.
"""
module Numerics

"""
    k_min(::Type{FT})

Floor for the two-stream diffusion eigenvalue `k = sqrt((γ1−γ2)(γ1+γ2))`,
applied under the square root. `k → 0` for isotropic, conservative scattering
(ssa → 1, the (γ1−γ2) factor vanishes); flooring `k²` at `sqrt(eps(FT))` keeps
`k` ≥ eps^(1/4) so the `e^{−2kτ}` cancellations in the diffuse
reflectance/transmittance stay resolvable. The relative error with respect to
the conservative-scattering limit is < 0.1% down to τ ~ 1e-9 (see the analogous
`min_k = 1e4·epsilon` bound and note in rte-rrtmgp, credited to Chiel van
Heerwaarden; swirl-lm floors `k` itself at 1e-2 for the same reason).
"""
@inline k_min(::Type{FT}) where {FT} = sqrt(eps(FT))

"""
    τ_thresh(::Type{FT})

Optical-depth threshold at which the longwave no-scattering source factor
`fact = (1 − T)/τ − T` (Clough et al. 1992, Eq 13) switches to its 3rd-order
series `τ(1/2 − τ/3 + τ²/8)`. The direct form loses ~eps/τ² to cancellation;
the series truncates at ~τ³. The two error curves cross at τ ~ eps^(1/4)
(≈ 1.9e-2 at Float32, ≈ 1.2e-4 at Float64) — matching rte-rrtmgp's
`tau_thresh = sqrt(sqrt(epsilon(tau)))`, credited to Peter Blossey and
Dmitry Alexeev.
"""
@inline τ_thresh(::Type{FT}) where {FT} = sqrt(sqrt(eps(FT)))

"""
    resonance_window(::Type{FT})

Half-width of the window around the removable singularity of the shortwave
direct reflectance/transmittance (Meador & Weaver 1980, Eqs 14–15) at
k·μ₀ = 1, inside which k·μ₀ is nudged off resonance. At the window edge
|1 − k²μ₀²| = sqrt(eps), the directly computed denominator still carries
≤ sqrt(eps) relative rounding, matching the O(sqrt(eps)) perturbation from
the nudge — the balanced choice.
"""
@inline resonance_window(::Type{FT}) where {FT} = sqrt(eps(FT))

"""
    μ₀_min(::Type{FT})

Floor for the cosine of the solar zenith angle wherever the shortwave solvers
divide by μ₀ (`exp(−τ/μ₀)` arguments). Columns with μ₀ ≤ 0 are excluded from
the solve and zeroed, so the floor only guards against division by a
zero/denormal μ₀ at the day–night terminator; the resulting exp argument
overflows negative and the flux underflows to zero, as it should. (The Fortran
reference uses sqrt(eps) for its round-earth path, where deep layers with
μ₀ ≤ 0 are computed nominally and masked afterwards; RRTMGP.jl skips those
columns instead, so a smaller floor suffices.)
"""
@inline μ₀_min(::Type{FT}) where {FT} = eps(FT)

"""
    pow_fast(x, y)

`x^y` via `exp(y·log(x))` for `x > 0`. Julia's generic `x^y` hits a slow path
for bases very close to 1, which the gray optical-depth profile evaluates
often.
"""
@inline pow_fast(x, y) = exp(y * log(x))

end # module
