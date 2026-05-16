"""
    canopy_2stream_props(R_b::FT, T_b::FT, G_val::FT, LAI_layer::FT) where {FT}

Map canopy leaf properties to TwoStream parameters for one layer and band.

Returns `(τ, ssa, g_eff)` where:
  - `τ = G(μ₀) × LAI_layer` (optical depth)
  - `ssa = R + T` (single scattering albedo)
  - `g_eff` maps the PIFM backward-scatter fraction to canopy `β = R/(R+T)`

The asymmetry parameter `g_eff` is derived from the PIFM relation between
`γ₁`, `γ₂` and the hemispheric backscatter fraction `β = R/(R+T)`:
  `g = 1 - β(8 - 2ω) / (3ω)`
"""
function canopy_2stream_props(R_b::FT, T_b::FT, G_val::FT, LAI_layer::FT) where {FT}
    τ = G_val * LAI_layer
    ω = R_b + T_b
    if ω > eps(FT)
        β = R_b / ω
        g_eff = clamp(FT(1) - β * (FT(8) - FT(2) * ω) / (FT(3) * ω), FT(-1), FT(1))
    else
        g_eff = FT(0)
    end
    return (τ, ω, g_eff)
end

"""
    canopy_sw_2stream_coeffs(τ::FT, ssa::FT, g::FT, μ₀::FT, β_dir::FT) where {FT}

Two-stream coefficients for canopy layers with physically-computed γ₃.

Like `sw_2stream_coeffs` but uses `β_dir` (directional backscatter from leaf angle
geometry) as γ₃ instead of the PIFM approximation `γ₃ = (2 - 3μ₀g)/4`.
γ₁ and γ₂ still use `g_eff` (correct for diffuse hemispherically-averaged scattering).

Returns `(Rdir, Tdir, Tnoscat, Rdif, Tdif)` — same tuple as `sw_2stream_coeffs`.
"""
function canopy_sw_2stream_coeffs(τ::FT, ssa::FT, g::FT, μ₀::FT, β_dir::FT) where {FT}
    k_min = sqrt(eps(FT))
    # PIFM γ₁, γ₂ from g_eff (hemispheric, correct for diffuse)
    γ1 = (FT(8) - ssa * (FT(5) + FT(3) * g)) * FT(0.25)
    γ2 = FT(3) * (ssa * (FT(1) - g)) * FT(0.25)
    # Direct beam: use β_dir from leaf angle geometry instead of PIFM
    γ3 = β_dir
    γ4 = FT(1) - γ3
    α1 = γ1 * γ4 + γ2 * γ3
    α2 = γ1 * γ3 + γ2 * γ4
    k = sqrt(max((γ1 - γ2) * (γ1 + γ2), k_min))

    exp_minusktau = exp(-τ * k)
    exp_minus2ktau = exp_minusktau * exp_minusktau

    RT_term = FT(1) / (k * (FT(1) + exp_minus2ktau) + γ1 * (FT(1) - exp_minus2ktau))

    Rdif = RT_term * γ2 * (FT(1) - exp_minus2ktau)
    Tdif = RT_term * FT(2) * k * exp_minusktau

    T₀ = Tnoscat = exp(-τ / max(μ₀, eps(FT)))

    k_μ = k * μ₀
    k_γ3 = k * γ3
    k_γ4 = k * γ4

    if abs(FT(1) - k_μ * k_μ) ≥ eps(FT)
        RT_term = ssa * RT_term / (FT(1) - k_μ * k_μ)
    else
        RT_term = ssa * RT_term / eps(FT)
    end

    Rdir_unconstrained =
        RT_term * (
            (FT(1) - k_μ) * (α2 + k_γ3) - (FT(1) + k_μ) * (α2 - k_γ3) * exp_minus2ktau -
            FT(2) * (k_γ3 - α2 * k_μ) * exp_minusktau * T₀
        )
    Tdir_unconstrained =
        -RT_term * (
            (FT(1) + k_μ) * (α1 + k_γ4) * T₀ - (FT(1) - k_μ) * (α1 - k_γ4) * exp_minus2ktau * T₀ -
            FT(2) * (k_γ4 + α1 * k_μ) * exp_minusktau
        )
    Rdir = max(FT(0), min(Rdir_unconstrained, (FT(1) - T₀)))
    Tdir = max(FT(0), min(Tdir_unconstrained, (FT(1) - T₀ - Rdir)))
    return (Rdir, Tdir, Tnoscat, Rdif, Tdif)
end

"""
    fill_canopy_sw_optics!(op_exp, op_atm, canopy, ibnd, gcol, nlay_atm)

Fill expanded column optical properties for one SW g-point and column.

Copies atmospheric optics from `op_atm` into the upper layers of `op_exp`
and fills the lower canopy layers from `canopy` leaf properties for band `ibnd`.
"""
function fill_canopy_sw_optics!(
    op_exp::TwoStream,
    op_atm::TwoStream,
    canopy::CanopyState,
    ibnd::Int,
    gcol::Int,
    nlay_atm::Int,
)
    ncanlay = canopy.ncanlay
    FT = eltype(op_exp.τ)

    # Copy atmospheric optics to upper layers of expanded column
    @inbounds for ilay in 1:nlay_atm
        op_exp.τ[ncanlay + ilay, gcol] = op_atm.τ[ilay, gcol]
        op_exp.ssa[ncanlay + ilay, gcol] = op_atm.ssa[ilay, gcol]
        op_exp.g[ncanlay + ilay, gcol] = op_atm.g[ilay, gcol]
    end

    # Fill canopy layers from leaf properties
    R_b = canopy.leaf_r[ibnd, gcol]
    T_b = canopy.leaf_t[ibnd, gcol]
    G_dir = canopy.G_dir[gcol]

    @inbounds for ilay in 1:ncanlay
        lai = canopy.lai_layer[ilay, gcol]
        τ_c, ω_c, g_c = canopy_2stream_props(FT(R_b), FT(T_b), FT(G_dir), FT(lai))
        op_exp.τ[ilay, gcol] = τ_c
        op_exp.ssa[ilay, gcol] = ω_c
        op_exp.g[ilay, gcol] = g_c
    end
    return nothing
end

"""
    fill_canopy_lw_optics!(op_exp, op_atm, canopy, ibnd_lw, gcol, nlay_atm)

Fill expanded column optical properties for one LW g-point and column.

Uses hemispheric-mean G for optical depth (LW is diffuse-only) and
`leaf_omega_lw` for single scattering albedo. Asymmetry parameter is zero
(isotropic LW scattering: leaf R_lw = T_lw).
"""
function fill_canopy_lw_optics!(
    op_exp::TwoStream,
    op_atm::TwoStream,
    canopy::CanopyState,
    ibnd_lw::Int,
    gcol::Int,
    nlay_atm::Int,
)
    ncanlay = canopy.ncanlay
    FT = eltype(op_exp.τ)

    # Copy atmospheric optics to upper layers
    @inbounds for ilay in 1:nlay_atm
        op_exp.τ[ncanlay + ilay, gcol] = op_atm.τ[ilay, gcol]
        op_exp.ssa[ncanlay + ilay, gcol] = op_atm.ssa[ilay, gcol]
        op_exp.g[ncanlay + ilay, gcol] = op_atm.g[ilay, gcol]
    end

    # Fill canopy LW optics
    G_h = canopy.G_hemi[gcol]
    ω_lw = canopy.leaf_omega_lw[ibnd_lw, gcol]

    @inbounds for ilay in 1:ncanlay
        lai = canopy.lai_layer[ilay, gcol]
        op_exp.τ[ilay, gcol] = FT(G_h) * FT(lai)
        op_exp.ssa[ilay, gcol] = FT(ω_lw)
        op_exp.g[ilay, gcol] = FT(0)  # isotropic scattering
    end
    return nothing
end

"""
    fill_canopy_lw_sources!(sf_exp, sf_atm, canopy, ibnd_lw, gcol, nlay_atm, t_planck, totplnk_band)

Fill expanded column Planck sources for LW canopy layers.

Copies atmospheric level sources to upper levels of the expanded column and
scales canopy level sources by `B(T_canopy) / B(T_ref)` to preserve the
planckfrac weighting from the atmospheric solver.
"""
function fill_canopy_lw_sources!(
    sf_exp::SourceLW2Str,
    sf_atm::SourceLW2Str,
    canopy::CanopyState,
    ibnd_lw::Int,
    gcol::Int,
    nlay_atm::Int,
    t_planck,
    totplnk_band,
)
    ncanlay = canopy.ncanlay
    FT = eltype(sf_exp.lev_source)

    # Copy atmospheric sources to upper levels of expanded column
    @inbounds for ilev in 1:(nlay_atm + 1)
        sf_exp.lev_source[ncanlay + ilev, gcol] = sf_atm.lev_source[ilev, gcol]
    end

    # Scale canopy levels by B(T_can)/B(T_ref)
    # T_ref = atmospheric bottom level temperature = canopy top
    T_ref = canopy.t_canopy_lev[ncanlay + 1, gcol]
    B_ref = planck_interp(FT(T_ref), t_planck, totplnk_band)

    @inbounds for ilev in 1:ncanlay
        T_can = canopy.t_canopy_lev[ilev, gcol]
        B_can = planck_interp(FT(T_can), t_planck, totplnk_band)
        sf_exp.lev_source[ilev, gcol] = sf_atm.lev_source[1, gcol] * (B_can / B_ref)
    end

    # Surface Planck source
    sf_exp.sfc_source[gcol] = sf_atm.sfc_source[gcol]
    return nothing
end
