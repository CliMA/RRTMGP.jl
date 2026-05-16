"""
    planck_interp(T::FT, t_planck, totplnk_band) where {FT}

Interpolate band-integrated Planck function B(T) from RRTMGP lookup tables.

`t_planck` is the temperature grid, `totplnk_band` is the total Planck emission
per band at each temperature grid point.
"""
function planck_interp(T::FT, t_planck, totplnk_band) where {FT}
    T < t_planck[1] && return FT(totplnk_band[1])
    T > t_planck[end] && return FT(totplnk_band[end])
    dT = t_planck[2] - t_planck[1]
    loc = min(floor(Int, (T - t_planck[1]) / dT) + 1, length(t_planck) - 1)
    frac = (T - t_planck[loc]) / dT
    return totplnk_band[loc] * (FT(1) - frac) + totplnk_band[loc + 1] * frac
end

"""
    compute_beta_dir(μ₀::FT, R_b::FT, T_b::FT, leaf_angle_fn; nθ=90, nφ=360) where {FT}

Direct-beam upscatter fraction β₀(μ₀) for a canopy layer.

Integrates over the leaf angle distribution to find the fraction of single-scattered
radiation from a collimated beam at zenith angle acos(μ₀) that goes upward.
For each leaf orientation (θ_L, φ_L), computes the fraction of Lambertian-reflected
and transmitted light entering the upper hemisphere.

`leaf_angle_fn(θ_L)` should return the leaf angle probability density at inclination `θ_L`.

Returns β₀ in [0,1]: fraction of scattered radiation going upward.
  - β₀ → R/(R+T) when μ₀ → 1 (overhead sun, horizontal leaves dominate)
  - β₀ → 0.5 when μ₀ → 0 (grazing sun, up/down symmetry)
"""
function compute_beta_dir(μ₀::FT, R_b::FT, T_b::FT, leaf_angle_fn; nθ=90, nφ=360) where {FT}
    ω = R_b + T_b
    ω < eps(FT) && return FT(0.5)

    sin_θ₀ = sqrt(max(FT(0), FT(1) - μ₀^2))
    dθ = FT(π / 2) / nθ
    dφ = FT(2π) / nφ

    num = FT(0)  # ∫ upscatter × |cos γ| × f_L dΩ
    den = FT(0)  # ∫ |cos γ| × f_L dΩ

    for iθ in 1:nθ
        θ_L = (iθ - FT(0.5)) * dθ
        sin_θ_L = sin(θ_L)
        cos_θ_L = cos(θ_L)
        f_L = FT(leaf_angle_fn(θ_L))

        for iφ in 1:nφ
            φ_L = (iφ - FT(0.5)) * dφ
            # Projection of beam onto leaf normal: cos γ = n̂·(-Ω̂₀)
            cos_γ = cos_θ_L * μ₀ - sin_θ_L * cos(φ_L) * sin_θ₀
            abs_cos_γ = abs(cos_γ)

            # Lambertian scattering from tilted leaf:
            # For reflected light: emitted from the hit surface
            #   cos γ > 0 (upper surface hit): normal up → (1+cos θ_L)/2 goes up
            #   cos γ < 0 (lower surface hit): normal down → (1-cos θ_L)/2 goes up
            # For transmitted light: emitted from the opposite surface
            s = cos_γ > 0 ? FT(1) : FT(-1)
            up_R = (FT(1) + s * cos_θ_L) / FT(2)
            up_T = (FT(1) - s * cos_θ_L) / FT(2)
            upscatter = R_b * up_R + T_b * up_T

            w = abs_cos_γ * f_L * dθ * dφ
            num += upscatter * w
            den += w
        end
    end

    # β₀ = total upscatter / (ω × total intercepted)
    return den > eps(FT) ? num / (ω * den) : FT(0.5)
end
