"""
    compute_col_gas_kernel!(
        col_gas,
        p_lev,
        mol_m_dry,
        mol_m_h2o,
        avogadro
        helmert1,
        vmr_h2o
        lat,
        glay, gcol,
    )

This function computes the column amounts of dry or moist air.
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
    glay::Int, # global lay id
    gcol::Int, # global col id
) where {FT <: AbstractFloat}
    helmert2 = FT(0.02586)  # second constant of Helmert formula
    m2_to_cm2 = FT(100 * 100) # m^2 to cm^2
    g0 =
        lat isa AbstractArray ? helmert1 - helmert2 * cos(FT(2) * FT(π) * lat[gcol] / FT(180)) : # acceleration due to gravity [m/s^2]
        helmert1
    Δp = p_lev[glay, gcol] - p_lev[glay + 1, gcol]
    vmr_h2o_glaygcol = vmr_h2o isa AbstractArray ? vmr_h2o[glay, gcol] : FT(0)
    # Get average mass of moist air per mole of moist air
    m_air = (mol_m_dry + mol_m_h2o * vmr_h2o_glaygcol)
    # Hydrostatic equation
    col_gas[glay, gcol] = (Δp * avogadro / (m2_to_cm2 * m_air * g0)) # molecules/cm^2
end
"""
    compute_interp_fractions(
        lkp::AbstractLookUp{FT},
        vmr,
        p_lay,
        t_lay,
        tropo,
        ibnd,
        glay, gcol,
    ) where {FT<:AbstractFloat}

compute interpolation fractions for binary species parameter, pressure and temperature.
"""
@inline function compute_interp_fractions(
    lkp::AbstractLookUp{FT},
    vmr,
    p_lay,
    t_lay,
    tropo,
    ibnd,
    glay,
    gcol,
) where {FT <: AbstractFloat}
    jftemp = compute_interp_frac_temp(lkp, t_lay, glay, gcol)
    jfpress = compute_interp_frac_press(lkp, p_lay, tropo, glay, gcol)
    jfη, col_mix = compute_interp_frac_η(lkp, vmr, tropo, jftemp[1], ibnd, glay, gcol)
    return (jftemp, jfpress, jfη, col_mix)
end

"""
    compute_interp_frac_temp(
        lkp::AbstractLookUp{FT},
        t_lay,
        glay,
        gcol,
    ) where {FT<:AbstractFloat}

compute interpolation fraction for temperature.
"""
@inline function compute_interp_frac_temp(lkp::AbstractLookUp{FT}, t_lay, glay, gcol) where {FT <: AbstractFloat}
    (; Δ_t_ref, n_t_ref, t_ref) = lkp

    @inbounds jtemp = loc_lower(t_lay, Δ_t_ref, n_t_ref, t_ref)
    @inbounds ftemp = (t_lay - t_ref[jtemp]) / Δ_t_ref

    return (jtemp, ftemp)
end

"""
    compute_interp_frac_press(
        lkp::AbstractLookUp,
        p_lay,
        tropo,
        glay,
        gcol,
    )

Compute interpolation fraction for pressure.
"""
@inline function compute_interp_frac_press(lkp::AbstractLookUp, p_lay, tropo, glay, gcol)
    (; Δ_ln_p_ref, p_ref, n_p_ref) = lkp

    @inbounds jpress = Int(min(max(fld(log(p_ref[1]) - log(p_lay), Δ_ln_p_ref) + 1, 1), n_p_ref - 1) + 1)

    @inbounds fpress = (log(p_ref[jpress - 1]) - log(p_lay)) / Δ_ln_p_ref
    jpress = jpress + tropo - 1

    return (jpress, fpress)
end

"""
    compute_interp_frac_η(
        lkp::AbstractLookUp{FT},
        vmr,
        tropo,
        jtemp,
        ibnd,
        glay,
        gcol,
    ) where {FT<:AbstractFloat}

Compute interpolation fraction for binary species parameter.
"""
@inline function compute_interp_frac_η(
    lkp::AbstractLookUp{FT},
    vmr,
    tropo,
    jtemp,
    ibnd,
    glay,
    gcol,
) where {FT <: AbstractFloat}
    (; n_η, key_species, vmr_ref) = lkp
    ig = view(key_species, :, tropo, ibnd)

    vmr1 = get_vmr(vmr, ig[1], glay, gcol)
    vmr2 = get_vmr(vmr, ig[2], glay, gcol)

    itemp = 1
    @inbounds η_half = vmr_ref[tropo, ig[1] + 1, jtemp + itemp - 1] / vmr_ref[tropo, ig[2] + 1, jtemp + itemp - 1]
    col_mix1 = vmr1 + η_half * vmr2
    η = col_mix1 ≥ eps(FT) * 2 ? vmr1 / col_mix1 : FT(0.5)
    loc_η = FT(η * (n_η - 1))
    jη1 = min(unsafe_trunc(Int, loc_η) + 1, n_η - 1)
    #fη1 = loc_η % FT(1) # TODO: "%: operator seems unstable on GPU
    #fη1 = FT(loc_η % 1) # to be revisited
    fη1 = loc_η - fld(loc_η, FT(1))

    itemp = 2
    @inbounds η_half = vmr_ref[tropo, ig[1] + 1, jtemp + itemp - 1] / vmr_ref[tropo, ig[2] + 1, jtemp + itemp - 1]
    col_mix2 = vmr1 + η_half * vmr2
    η = col_mix2 ≥ eps(FT) * 2 ? vmr1 / col_mix2 : FT(0.5)
    loc_η = FT(η * (n_η - 1))
    jη2 = min(unsafe_trunc(Int, loc_η) + 1, n_η - 1)
    #fη2 = loc_η % FT(1) # TODO: "%" operator seems unstable on GPU
    #fη2 = FT(loc_η % 1) # to be revisited
    fη2 = loc_η - fld(loc_η, FT(1))

    return ((jη1, jη2, fη1, fη2), (col_mix1, col_mix2))#nothing
end

"""
    compute_τ_ssa_lw_src!(
        lkp::AbstractLookUp{FT},
        vmr,
        col_dry,
        igpt,
        ibnd,
        p_lay::FT,
        t_lay,
        glay, gcol,
        src_args...,
    ) where {FT<:AbstractFloat}

Compute optical thickness, single scattering albedo, asymmetry parameter 
and longwave sources whenever applicable.
"""
@inline function compute_τ_ssa_lw_src!(
    lkp::AbstractLookUp{FT},
    vmr,
    col_dry,
    igpt,
    ibnd,
    p_lay::FT,
    t_lay,
    glay,
    gcol,
    src_args...,
) where {FT <: AbstractFloat}
    # upper/lower troposphere
    tropo = p_lay > lkp.p_ref_tropo ? 1 : 2
    # volume mixing ratio of h2o
    vmr_h2o = get_vmr(vmr, lkp.idx_h2o, glay, gcol)

    jftemp, jfpress, jfη, col_mix = compute_interp_fractions(lkp, vmr, p_lay, t_lay, tropo, ibnd, glay, gcol)

    # computing τ_major
    τ_major = interp3d(jfη..., jftemp..., jfpress..., lkp.kmajor, igpt, col_mix...) * col_dry
    # computing τ_minor
    τ_minor =
        compute_τ_minor(lkp, tropo, vmr, vmr_h2o, col_dry, p_lay, t_lay, jftemp..., jfη..., igpt, ibnd, glay, gcol)
    # compute τ_Rayleigh
    τ_ray = compute_τ_rayleigh(lkp, tropo, col_dry, vmr_h2o, jftemp..., jfη..., igpt)
    τ = τ_major + τ_minor + τ_ray
    ssa = FT(0)
    if τ > 2 * eps(FT) # single scattering albedo
        ssa = τ_ray / τ
    end
    # computing Planck sources for longwave problem
    compute_lw_planck_src!(lkp, jfη..., jfpress..., jftemp..., t_lay, igpt, ibnd, glay, gcol, src_args...)
    return (τ, ssa)
end

"""
    compute_τ_minor(
        lkp::AbstractLookUp,
        tropo::Int,
        vmr,
        vmr_h2o::FT,
        col_dry,
        p_lay::FT,
        t_lay::FT,
        jtemp::Int,
        ftemp::FT,
        jη1::Int,
        jη2::Int,
        fη1::FT,
        fη2::FT,
        igpt,
        ibnd,
        glay,
        gcol,
    ) where {FT<:AbstractFloat}

Compute optical thickness contributions from minor gases.
"""
@inline function compute_τ_minor(
    lkp::AbstractLookUp,
    tropo::Int,
    vmr,
    vmr_h2o::FT,
    col_dry,
    p_lay::FT,
    t_lay::FT,
    jtemp::Int,
    ftemp::FT,
    jη1::Int,
    jη2::Int,
    fη1::FT,
    fη2::FT,
    igpt,
    ibnd,
    glay,
    gcol,
) where {FT <: AbstractFloat}

    if tropo == 1 # in lower atmosphere
        minor_bnd_st = lkp.minor_lower_bnd_st
        idx_gases_minor = lkp.idx_gases_minor_lower
        minor_scales_with_density = lkp.minor_lower_scales_with_density
        idx_scaling_gas = lkp.idx_scaling_gas_lower
        scale_by_complement = lkp.lower_scale_by_complement
        minor_gpt_sh = lkp.minor_lower_gpt_sh
        kminor = lkp.kminor_lower
    else # in upper atmosphere
        minor_bnd_st = lkp.minor_upper_bnd_st
        idx_gases_minor = lkp.idx_gases_minor_upper
        minor_scales_with_density = lkp.minor_upper_scales_with_density
        idx_scaling_gas = lkp.idx_scaling_gas_upper
        scale_by_complement = lkp.upper_scale_by_complement
        minor_gpt_sh = lkp.minor_upper_gpt_sh
        kminor = lkp.kminor_upper
    end

    τ_minor = FT(0)
    pa2hpa = FT(0.01) # pascals to hectopascals
    dry_fact = FT(1) / (FT(1) + vmr_h2o)

    @inbounds loc_in_bnd = igpt - (lkp.bnd_lims_gpt[1, ibnd] - 1)

    @inbounds for i in minor_bnd_st[ibnd]:(minor_bnd_st[ibnd + 1] - 1)
        vmr_imnr = get_vmr(vmr, idx_gases_minor[i], glay, gcol)
        if vmr_imnr > eps(FT) * 2
            scaling = vmr_imnr * col_dry

            if minor_scales_with_density[i] == 1
                scaling *= (pa2hpa * p_lay / t_lay)
                sgas = idx_scaling_gas[i]
                if sgas > 0
                    if scale_by_complement[i] == 1
                        scaling *= (FT(1) - get_vmr(vmr, sgas, glay, gcol) * dry_fact)
                    else
                        scaling *= get_vmr(vmr, sgas, glay, gcol) * dry_fact
                    end
                end
            end
            k_loc = minor_gpt_sh[i] + loc_in_bnd
            τ_minor += interp2d(fη1, fη2, ftemp, kminor, k_loc, jη1, jη2, jtemp) * scaling
        end
    end

    return τ_minor
end

"""
    compute_τ_rayleigh(
        lkp::LookUpSW,
        tropo::Int,
        col_dry::FT,
        vmr_h2o::FT,
        jtemp::Int,
        ftemp::FT,
        jη1::Int,
        jη2::Int,
        fη1::FT,
        fη2::FT,
        igpt::Int,
    ) where {FT<:AbstractFloat}

Compute Rayleigh scattering optical depths for shortwave problem
"""
@inline function compute_τ_rayleigh(
    lkp::LookUpSW,
    tropo::Int,
    col_dry::FT,
    vmr_h2o::FT,
    jtemp::Int,
    ftemp::FT,
    jη1::Int,
    jη2::Int,
    fη1::FT,
    fη2::FT,
    igpt::Int,
) where {FT <: AbstractFloat}
    if tropo == 1
        τ_ray = interp2d(fη1, fη2, ftemp, lkp.rayl_lower, igpt, jη1, jη2, jtemp) * (vmr_h2o + FT(1)) * col_dry
    else
        τ_ray = interp2d(fη1, fη2, ftemp, lkp.rayl_upper, igpt, jη1, jη2, jtemp) * (vmr_h2o + FT(1)) * col_dry
    end

    return τ_ray
end

@inline function compute_τ_rayleigh(lkp::LookUpLW{FT}, args...) where {FT <: AbstractFloat}
    return FT(0)
end

"""
    compute_lw_planck_src!(
        lkp::LookUpLW,
        jη1,
        jη2,
        fη1,
        fη2,
        jpresst,
        fpress,
        jtemp,
        ftemp,
        t_lay,
        igpt,
        ibnd,
        glay,
        gcol,
        sf,
        t_lev,
        t_sfc,
    )

Computes Planck sources for the longwave problem.

"""
@inline function compute_lw_planck_src!(
    lkp::LookUpLW,
    jη1,
    jη2,
    fη1,
    fη2,
    jpresst,
    fpress,
    jtemp,
    ftemp,
    t_lay,
    igpt,
    ibnd,
    glay,
    gcol,
    sf,
    t_lev,
    t_sfc,
)
    (; planck_fraction, t_planck, n_t_plnk, totplnk) = lkp
    (; lay_source, lev_source_inc, lev_source_dec, sfc_source) = sf
    # compute Planck fraction
    p_frac = interp3d(jη1, jη2, fη1, fη2, jtemp, ftemp, jpresst, fpress, planck_fraction, igpt)

    planck_args = (t_planck, totplnk, ibnd)
    # computing lay_source
    @inbounds lay_source[glay, gcol] = interp1d(t_lay, planck_args...) * p_frac
    # computing lev_source_inc
    @inbounds lev_source_inc[glay, gcol] = interp1d(t_lev[glay + 1, gcol], planck_args...) * p_frac
    # computing lev_source_dec
    @inbounds lev_source_dec[glay, gcol] = interp1d(t_lev[glay, gcol], planck_args...) * p_frac
    if glay == 1 # computing sfc_source
        @inbounds sfc_source[gcol] = interp1d(t_sfc, planck_args...) * p_frac
    end
    return nothing
end

@inline function compute_lw_planck_src!(lkp::LookUpSW, args...)
    return nothing
end
