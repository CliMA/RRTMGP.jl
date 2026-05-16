module RRTMGPCanopyOpticsExt

import CanopyOptics
import CanopyOptics: LeafDistribution, LeafProspectProProperties, LeafOpticalProperties
import CanopyOptics: prospect, createLeafOpticalStruct, G
using Unitful: ustrip, @u_str

using RRTMGP
using RRTMGP.Canopy: CanopyState, compute_beta_dir
import RRTMGP.Canopy: fill_leaf_optics!, fill_G!, fill_beta_dir!

# Access pdf from Distributions via CanopyOptics (which depends on it)
const _pdf = CanopyOptics.Distributions.pdf

"""
    fill_leaf_optics!(canopy::CanopyState, leaf, opti, bnd_lims_wn)

Band-average PROSPECT leaf R/T onto RRTMGP SW bands and fill `canopy.leaf_r`, `canopy.leaf_t`.
Also computes per-band chlorophyll and carotenoid absorption fractions (`f_cab_band`, `f_car_band`).

`leaf` is a `LeafProspectProProperties`, `opti` is from `createLeafOpticalStruct`,
and `bnd_lims_wn` is `(2, nbnd)` array of band wavenumber limits in cm⁻¹.

Returns a `BitVector` indicating which bands overlap with the PROSPECT wavelength range.
"""
function fill_leaf_optics!(
    canopy::CanopyState,
    leaf::LeafProspectProProperties,
    opti::LeafOpticalProperties,
    bnd_lims_wn::AbstractMatrix,
)
    T_leaf, R_leaf = prospect(leaf, opti)
    λ_nm = ustrip.(u"nm", opti.λ)

    # Compute per-wavelength component absorption fractions
    (; N, Ccab, Ccar, Canth, Cbrown, Cw, Cm, Cprot, Ccbc) = leaf
    (; Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc) = opti
    K_unnorm = Ccab .* Kcab .+ Ccar .* Kcar .+ Canth .* Kant .+ Cbrown .* Kb .+
               Cw .* Kw .+ Cm .* Km .+ Cprot .* Kp .+ Ccbc .* Kcbc
    FT = eltype(canopy.leaf_r)
    f_cab_λ = (Ccab .* Kcab) ./ max.(K_unnorm, eps(FT))
    f_car_λ = (Ccar .* Kcar) ./ max.(K_unnorm, eps(FT))
    α_λ = FT(1) .- R_leaf .- T_leaf  # total leaf absorptance

    nbnd = size(bnd_lims_wn, 2)
    ncol = size(canopy.leaf_r, 2)
    band_in_prospect = falses(nbnd)

    for ibnd in 1:nbnd
        wn_lo, wn_hi = bnd_lims_wn[1, ibnd], bnd_lims_wn[2, ibnd]
        λ_hi_nm = 1e7 / wn_lo
        λ_lo_nm = 1e7 / wn_hi

        idx = findall(λ_lo_nm .≤ λ_nm .≤ λ_hi_nm)

        if !isempty(idx)
            R_avg = sum(R_leaf[idx]) / length(idx)
            T_avg = sum(T_leaf[idx]) / length(idx)
            band_in_prospect[ibnd] = true

            # Absorptance-weighted component fractions
            α_sum = sum(α_λ[idx])
            f_cab_avg = α_sum > eps(FT) ? sum(f_cab_λ[idx] .* α_λ[idx]) / α_sum : FT(0)
            f_car_avg = α_sum > eps(FT) ? sum(f_car_λ[idx] .* α_λ[idx]) / α_sum : FT(0)

            for gcol in 1:ncol
                canopy.leaf_r[ibnd, gcol] = R_avg
                canopy.leaf_t[ibnd, gcol] = T_avg
                canopy.f_cab_band[ibnd, gcol] = f_cab_avg
                canopy.f_car_band[ibnd, gcol] = f_car_avg
            end
        end
    end
    return band_in_prospect
end

"""
    fill_G!(canopy::CanopyState, ld::LeafDistribution, cos_zenith)

Fill `G_dir` from leaf angle distribution and zenith angle, set `G_hemi` to 0.5.
"""
function fill_G!(canopy::CanopyState, ld::LeafDistribution, cos_zenith::AbstractVector)
    ncol = length(cos_zenith)
    G_vals = G(cos_zenith, ld)
    for gcol in 1:ncol
        canopy.G_dir[gcol] = G_vals[gcol]
        canopy.G_hemi[gcol] = eltype(canopy.G_hemi)(0.5)
    end
    return nothing
end

"""
    fill_beta_dir!(canopy::CanopyState, ld::LeafDistribution, cos_zenith; nθ=90, nφ=360)

Precompute directional backscatter fraction `β_dir` for each SW band and column.
"""
function fill_beta_dir!(
    canopy::CanopyState,
    ld::LeafDistribution,
    cos_zenith::AbstractVector;
    nθ = 90,
    nφ = 360,
)
    FT = eltype(canopy.beta_dir)
    nbnd = size(canopy.beta_dir, 1)
    ncol = size(canopy.beta_dir, 2)

    # Build leaf angle density callable from LeafDistribution
    leaf_angle_fn = θ_L -> FT(_pdf(ld.LD, 2θ_L / π) * ld.scaling)

    for gcol in 1:ncol
        μ₀ = FT(cos_zenith[gcol])
        for ibnd in 1:nbnd
            R_b = FT(canopy.leaf_r[ibnd, gcol])
            T_b = FT(canopy.leaf_t[ibnd, gcol])
            ω = R_b + T_b
            if ω > eps(FT)
                canopy.beta_dir[ibnd, gcol] = compute_beta_dir(μ₀, R_b, T_b, leaf_angle_fn; nθ, nφ)
            else
                canopy.beta_dir[ibnd, gcol] = FT(0.5)
            end
        end
    end
    return nothing
end

end # module
