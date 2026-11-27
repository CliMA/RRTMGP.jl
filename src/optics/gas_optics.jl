"""
    compute_col_gas_kernel!(
        col_gas,
        p_lev,
        mol_m_dry,
        mol_m_h2o,
        avogadro,
        helmert1,
        vmr_h2o,
        lat,
        glay, gcol,
    )

Compute the column amounts of dry or moist air.
"""
function compute_col_gas_kernel!(
    col_gas::AbstractArray{FT, 2},
    p_lev::AbstractArray{FT, 2},
    mol_m_dry::FT,
    mol_m_h2o::FT,
    avogadro::FT,
    helmert1::FT,
    vmr_h2o::Union{AbstractArray{FT, 2}, Nothing},
    lat::Union{AbstractArray{FT, 1}, Nothing},
    glay::Integer, # global lay id
    gcol::Integer, # global col id
) where {FT <: AbstractFloat}
    helmert2 = FT(0.02586)  # second constant of Helmert formula
    m2_to_cm2 = FT(100 * 100) # m^2 to cm^2
    g0 =
        lat isa AbstractArray ?
        helmert1 - helmert2 * cos(FT(2) * FT(π) * lat[gcol] / FT(180)) : # acceleration due to gravity [m/s^2]
        helmert1
    Δp = p_lev[glay, gcol] - p_lev[glay + 1, gcol]
    vmr_h2o_glaygcol = vmr_h2o isa AbstractArray ? vmr_h2o[glay, gcol] : FT(0)
    # Get average mass of moist air per mole of moist air
    m_air = (mol_m_dry + mol_m_h2o * vmr_h2o_glaygcol)
    # Hydrostatic equation
    col_gas[glay, gcol] = (Δp * avogadro / (m2_to_cm2 * m_air * g0)) # molecules/cm^2
    return nothing
end

"""
compute_relative_humidity_kernel!(
    rh::AbstractArray{FT, 2},
    p_lay::AbstractArray{FT, 2},
    t_lay::AbstractArray{FT, 2},
    vmr_h2o::AbstractArray{FT, 2},
    mwd::FT,
    t_ref::FT,
    q_lay_min::FT,
    glay::Integer,
    gcol::Integer,
) where {FT}

Compute the relative humidity.
"""
function compute_relative_humidity_kernel!(
    rh::AbstractArray{FT, 2},
    p_lay::AbstractArray{FT, 2},
    t_lay::AbstractArray{FT, 2},
    vmr_h2o::AbstractArray{FT, 2},
    mwd::FT,
    t_ref::FT,
    q_lay_min::FT,
    glay::Integer,
    gcol::Integer,
) where {FT}
    # Convert h2o vmr to mmr
    mmr_h2o = vmr_h2o[glay, gcol] * mwd
    q_lay = mmr_h2o / (FT(1) + mmr_h2o)
    q_tmp = max(q_lay_min, q_lay)
    es_tmp = exp(
        (FT(17.67) * (t_lay[glay, gcol] - t_ref)) /
        (t_lay[glay, gcol] - FT(29.65)),
    )
    rh[glay, gcol] =
        max(FT(0.01) * (FT(0.263) * p_lay[glay, gcol] * q_tmp) / es_tmp, FT(0))
    return nothing
end

"""
    compute_interp_frac_temp(t_ref, t_lay)

Compute interpolation fraction for temperature.
"""
@inline function compute_interp_frac_temp(t_ref, t_lay)
    @inbounds Δ_t_ref = t_ref[2] - t_ref[1]
    n_t_ref = length(t_ref)
    @inbounds jtemp = loc_lower(t_lay, Δ_t_ref, n_t_ref, t_ref)
    @inbounds ftemp = (t_lay - t_ref[jtemp]) / Δ_t_ref
    return (jtemp, ftemp)
end

"""
    compute_interp_frac_press(ln_p_ref, p_lay, tropo)

Compute interpolation fraction for pressure.
"""
@inline function compute_interp_frac_press(ln_p_ref, p_lay, tropo)
    @inbounds Δ_ln_p_ref = ln_p_ref[1] - ln_p_ref[2]
    n_p_ref = length(ln_p_ref)
    log_p_lay = log(p_lay)
    @inbounds jpress = Int(
        min(
            max(
                unsafe_trunc(Int, (ln_p_ref[1] - log_p_lay) / Δ_ln_p_ref) + 1,
                1,
            ),
            n_p_ref - 1,
        ) + 1,
    )
    @inbounds fpress = (ln_p_ref[jpress - 1] - log_p_lay) / Δ_ln_p_ref
    return (jpress + tropo - 1, fpress)
end

"""
    compute_interp_frac_η(
        n_η,
        ig,
        vmr_ref,
        (vmr1, vmr2),
        tropo,
        jtemp,
    )

Compute interpolation fraction for the binary species parameter ``η``.
"""
@inline function compute_interp_frac_η(
    n_η,
    ig,
    vmr_ref,
    (vmr1, vmr2),
    tropo,
    jtemp,
)
    FT = eltype(vmr1)
    itemp = 1
    @inbounds η_half =
        vmr_ref[tropo, ig[1] + 1, jtemp + itemp - 1] /
        vmr_ref[tropo, ig[2] + 1, jtemp + itemp - 1]
    col_mix1 = vmr1 + η_half * vmr2
    η = vmr1 * inv(col_mix1) # rewritten to ease register pressure
    if col_mix1 ≤ FT(0)
        η = FT(0.5)
    end
    #η = col_mix1 > 0 ? vmr1 / col_mix1 : FT(0.5) # rte-rrtmgp uses col_mix1 > tiny(col_mix1)
    loc_η = FT(η * (n_η - 1))
    jη1 = min(unsafe_trunc(Int, loc_η) + 1, n_η - 1)
    # Fraction relative to the (clamped) cell, so η = 1 lands on the last table
    # node (fη = 1) instead of jumping back a full cell (fη = 0). At Float32, η
    # rounds to exactly 1 far more often than at Float64.
    fη1 = loc_η - (jη1 - 1)

    itemp = 2
    @inbounds η_half =
        vmr_ref[tropo, ig[1] + 1, jtemp + itemp - 1] /
        vmr_ref[tropo, ig[2] + 1, jtemp + itemp - 1]
    col_mix2 = vmr1 + η_half * vmr2
    η = vmr1 * inv(col_mix2) # rewritten to ease register pressure
    if col_mix2 ≤ FT(0)
        η = FT(0.5)
    end
    #η = col_mix2 > 0 ? vmr1 / col_mix2 : FT(0.5) # rte-rrtmgp uses col_mix2 > tiny(col_mix2)
    loc_η = FT(η * (n_η - 1))
    jη2 = min(unsafe_trunc(Int, loc_η) + 1, n_η - 1)
    fη2 = loc_η - (jη2 - 1) # cell-relative; continuous at η = 1 (see fη1)

    return ((jη1, jη2, fη1, fη2), (col_mix1, col_mix2))
end

# The interpolation of layer `glay`'s gas optics that all g-points of band `ibnd`
# share: the troposphere flag (1 below the tropopause), the water vapor mixing
# ratio, the indices and fractions in temperature, pressure and `η`, and the
# column amounts of the band's key species. The fused GPU solves compute it once
# per band and layer.
@inline function gas_layer_interp(lkp, vmr, ibnd, p_lay, t_lay, glay, gcol)
    tropo = p_lay > lkp.p_ref_tropo ? 1 : 2
    vmr_h2o = get_vmr(vmr, lkp.idx_h2o, glay, gcol)
    jftemp = compute_interp_frac_temp(lkp.ref_points.t_ref, t_lay)
    jfpress = compute_interp_frac_press(lkp.ref_points.ln_p_ref, p_lay, tropo)
    ig = view(lkp.key_species, 1:2, tropo, ibnd)
    @inbounds vmr1 = get_vmr(vmr, ig[1], glay, gcol)
    @inbounds vmr2 = get_vmr(vmr, ig[2], glay, gcol)
    jfη, col_mix = compute_interp_frac_η(
        size(lkp.kmajor, 1),
        ig,
        lkp.ref_points.vmr_ref,
        (vmr1, vmr2),
        tropo,
        jftemp[1],
    )
    return (tropo, vmr_h2o, jftemp, jfpress, jfη, col_mix)
end

# The major and minor optical depths of g-point `igpt`, the `loc`-th of its band,
# in a layer interpolated as `layer` (see `gas_layer_interp`), from the minor
# contributors `bounds = (st, n)` of the layer's part of the atmosphere, scaled by
# `scaling(imnr)` (see `compute_τ_minor`).
@inline function gas_τ(
    lkp,
    (tropo, _, jftemp, jfpress, jfη, col_mix),
    col_dry,
    igpt,
    loc,
    (st, n),
    scaling::F,
    skip_absent::Val,
) where {F}
    @inbounds kmajor = view(lkp.kmajor, :, :, :, igpt)
    τ_major =
        interp3d(jfη..., jftemp..., jfpress..., kmajor, col_mix...) * col_dry
    lkp_minor = tropo == 1 ? lkp.minor_lower : lkp.minor_upper
    τ_minor = compute_τ_minor(
        lkp_minor,
        jftemp...,
        jfη...,
        loc,
        st,
        n,
        scaling,
        skip_absent,
    )
    return τ_major, τ_minor
end

# The Planck fraction of g-point `igpt` in a layer interpolated as `layer`.
@inline planck_fraction(lkp::LookUpLW, (_, _, jftemp, jfpress, jfη, _), igpt) =
    @inbounds interp3d(
        jfη...,
        jftemp...,
        jfpress...,
        view(lkp.planck.planck_fraction, :, :, :, igpt),
    )

# The Rayleigh optical depth of g-point `igpt` in a layer interpolated as `layer`.
@inline function rayleigh_τ(
    lkp::LookUpSW,
    (tropo, vmr_h2o, jftemp, _, jfη, _),
    col_dry,
    igpt,
)
    @inbounds rayleigh_coeff =
        tropo == 1 ? view(lkp.rayl_lower, :, :, igpt) :
        view(lkp.rayl_upper, :, :, igpt)
    return compute_τ_rayleigh(
        rayleigh_coeff,
        col_dry,
        vmr_h2o,
        jftemp...,
        jfη...,
        igpt,
    )
end

# The shortwave gas optical depth and single-scattering albedo.
@inline function gas_τ_ssa(τ_major, τ_minor, τ_ray)
    FT = typeof(τ_major)
    τ = max(τ_major + τ_minor + τ_ray, FT(0))
    ssa = τ_ray * inv(τ) # rewritten to ease register pressure
    if τ ≤ FT(0)
        ssa = FT(0)
    end
    return τ, ssa
end

# The major and minor optical depths of g-point `igpt` of band `ibnd` in a layer
# interpolated as `layer`, with the minor scalings computed on the fly.
@inline function serial_gas_τ(
    lkp,
    layer,
    vmr,
    col_dry,
    igpt,
    ibnd,
    p_lay,
    t_lay,
    glay,
    gcol,
)
    tropo, vmr_h2o = layer
    lkp_minor = tropo == 1 ? lkp.minor_lower : lkp.minor_upper
    # @inline, so that CPU solves hoist the loop-invariant divisions of
    # `minor_scaling` out of the loop over minor contributors
    @inline scaling(imnr) = minor_scaling(
        lkp_minor,
        vmr,
        vmr_h2o,
        col_dry,
        p_lay,
        t_lay,
        imnr,
        glay,
        gcol,
    )
    return gas_τ(
        lkp,
        layer,
        col_dry,
        igpt,
        igpt - lkp.band_data.bnd_lims_gpt[1, ibnd] + 1,
        LookUpTables.minor_bounds(lkp_minor, ibnd),
        scaling,
        Val(true),
    )
end

"""
    compute_gas_optics(
        lkp::Union{LookUpLW, LookUpSW},
        vmr,
        col_dry,
        igpt,
        ibnd,
        p_lay,
        t_lay,
        glay, gcol,
    )

Compute optical thickness, single scattering albedo, and asymmetry parameter.
"""
@inline function compute_gas_optics(
    lkp::LookUpLW,
    vmr,
    col_dry,
    igpt,
    ibnd,
    p_lay,
    t_lay,
    glay,
    gcol,
)
    layer = gas_layer_interp(lkp, vmr, ibnd, p_lay, t_lay, glay, gcol)
    τ_major, τ_minor = serial_gas_τ(
        lkp,
        layer,
        vmr,
        col_dry,
        igpt,
        ibnd,
        p_lay,
        t_lay,
        glay,
        gcol,
    )
    # τ_Rayleigh is zero for longwave
    τ = max(τ_major + τ_minor, zero(τ_major))
    return (τ, zero(τ), zero(τ), planck_fraction(lkp, layer, igpt))
end

@inline function compute_gas_optics(
    lkp::LookUpSW,
    vmr,
    col_dry,
    igpt,
    ibnd,
    p_lay,
    t_lay,
    glay,
    gcol,
)
    layer = gas_layer_interp(lkp, vmr, ibnd, p_lay, t_lay, glay, gcol)
    τ_major, τ_minor = serial_gas_τ(
        lkp,
        layer,
        vmr,
        col_dry,
        igpt,
        ibnd,
        p_lay,
        t_lay,
        glay,
        gcol,
    )
    τ, ssa = gas_τ_ssa(τ_major, τ_minor, rayleigh_τ(lkp, layer, col_dry, igpt))
    return (τ, ssa, zero(τ))
end

"""
    minor_scaling(lkp_minor, vmr, vmr_h2o, col_dry, p_lay, t_lay, imnr, glay, gcol)

The factor by which the absorption coefficient of minor contributor `imnr` of
`lkp_minor` is scaled in layer `glay` of column `gcol`: the column amount of its
gas, scaled with density and by a scaling gas or its complement as the table
directs, or zero where the gas is absent. It depends on the layer but not on the
g-point, so the fused solves compute it once per band and layer.
"""
@inline function minor_scaling(
    lkp_minor::LookUpMinor,
    vmr,
    vmr_h2o,
    col_dry,
    p_lay,
    t_lay,
    imnr,
    glay,
    gcol,
)
    FT = typeof(col_dry)
    idx_gas, idx_scaling_gas, scales_with_density, scale_by_complement =
        LookUpTables.get_minor_gas_data(lkp_minor, imnr)
    vmr_imnr = get_vmr(vmr, idx_gas, glay, gcol)
    vmr_imnr > 0 || return FT(0)
    pa2hpa = FT(0.01) # pascals to hectopascals
    dry_fact = inv(FT(1) + vmr_h2o)
    density_fact = pa2hpa * p_lay / t_lay
    scaling = vmr_imnr * col_dry
    if scales_with_density == 1
        scaling *= density_fact
        if idx_scaling_gas > 0
            if scale_by_complement == 1
                scaling *= (
                    FT(1) -
                    get_vmr(vmr, idx_scaling_gas, glay, gcol) * dry_fact
                )
            else
                scaling *= get_vmr(vmr, idx_scaling_gas, glay, gcol) * dry_fact
            end
        end
    end
    return scaling
end

"""
    compute_τ_minor(lkp_minor, jtemp, ftemp, jη1, jη2, fη1, fη2, loc, st, n, scaling, skip_absent)

The optical depth of the minor contributors `st` to `st + n - 1` of `lkp_minor`
at the `loc`-th g-point of their band, each scaled by `scaling(imnr)` (see
[`minor_scaling`](@ref)) and summed in order. With `skip_absent = Val(true)`,
the contributors that `scaling` scales by zero are skipped; otherwise they are
summed two at a time, so that their table loads are in flight together.
"""
@inline function compute_τ_minor(
    lkp_minor::LookUpMinor,
    jtemp,
    ftemp::FT,
    jη1,
    jη2,
    fη1::FT,
    fη2::FT,
    loc,
    st,
    n,
    scaling::F,
    ::Val{skip_absent},
) where {FT, F, skip_absent}
    # @inline, so that CPU solves do not call it for each contributor
    @inline term(imnr, s) = @inbounds interp2d(
        fη1,
        fη2,
        ftemp,
        view(
            lkp_minor.kminor,
            :,
            :,
            LookUpTables.kminor_index(lkp_minor, imnr, loc),
        ),
        jη1,
        jη2,
        jtemp,
    ) * s
    τ_minor = FT(0)
    if skip_absent
        for i in zero(n):(n - one(n))
            s = scaling(st + i)
            iszero(s) && continue
            τ_minor += term(st + i, s)
        end
    else
        i = zero(n)
        while i + one(n) < n
            τ1 = term(st + i, scaling(st + i))
            τ2 = term(st + i + one(n), scaling(st + i + one(n)))
            τ_minor += τ1
            τ_minor += τ2
            i += 2 * one(n)
        end
        i < n && (τ_minor += term(st + i, scaling(st + i)))
    end
    return τ_minor
end

"""
    compute_τ_rayleigh(
        rayleigh_coeff::AbstractArray{FT, 2},
        col_dry::FT,
        vmr_h2o::FT,
        jtemp::Int,
        ftemp::FT,
        jη1::Int,
        jη2::Int,
        fη1::FT,
        fη2::FT,
        igpt::Integer,
    ) where {FT<:AbstractFloat}

Compute Rayleigh scattering optical depths for shortwave problem
"""
@inline compute_τ_rayleigh(
    rayleigh_coeff::AbstractArray{FT, 2},
    col_dry::FT,
    vmr_h2o::FT,
    jtemp::Int,
    ftemp::FT,
    jη1::Int,
    jη2::Int,
    fη1::FT,
    fη2::FT,
    igpt::Integer,
) where {FT <: AbstractFloat} =
    interp2d(fη1, fη2, ftemp, rayleigh_coeff, jη1, jη2, jtemp) *
    (vmr_h2o + FT(1)) *
    col_dry
