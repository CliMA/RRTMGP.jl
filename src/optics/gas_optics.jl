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
    Δp = p_lev[gcol, glay] - p_lev[gcol, glay + 1]
    vmr_h2o_local = vmr_h2o isa AbstractArray ? vmr_h2o[gcol, glay] : FT(0)
    # Get average mass of moist air per mole of moist air
    m_air = (mol_m_dry + mol_m_h2o * vmr_h2o_local)
    # Hydrostatic equation
    col_gas[gcol, glay] = (Δp * avogadro / (m2_to_cm2 * m_air * g0)) # molecules/cm^2
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
    glay::Int,
    gcol::Int,
) where {FT}

This function computes the relative humidity.
"""
function compute_relative_humidity_kernel!(
    rh::AbstractArray{FT, 2},
    p_lay::AbstractArray{FT, 2},
    t_lay::AbstractArray{FT, 2},
    vmr_h2o::AbstractArray{FT, 2},
    mwd::FT,
    t_ref::FT,
    q_lay_min::FT,
    glay::Int,
    gcol::Int,
) where {FT}
    mmr_h2o =
    # Convert h2o vmr to mmr
        mmr_h2o = vmr_h2o[gcol, glay] * mwd
    q_lay = mmr_h2o / (FT(1) + mmr_h2o)
    q_tmp = max(q_lay_min, q_lay)
    es_tmp = exp((FT(17.67) * (t_lay[gcol, glay] - t_ref)) / (t_lay[gcol, glay] - FT(29.65)))
    rh[gcol, glay] = max(FT(0.01) * (FT(0.263) * p_lay[gcol, glay] * q_tmp) / es_tmp, FT(0))
    return nothing
end

"""
    compute_interp_frac_temp(Δ_t_ref, n_t_ref, t_ref, t_lay)

compute interpolation fraction for temperature.
"""
@inline function compute_interp_frac_temp(t_ref, t_lay)
    @inbounds Δ_t_ref = t_ref[2] - t_ref[1]
    n_t_ref = length(t_ref)
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
@inline function compute_interp_frac_press(ln_p_ref, p_lay, tropo)
    @inbounds Δ_ln_p_ref = ln_p_ref[1] - ln_p_ref[2]
    n_p_ref = length(ln_p_ref)
    log_p_lay = log(p_lay)
    @inbounds jpress = Int(min(max(unsafe_trunc(Int, (ln_p_ref[1] - log_p_lay) / Δ_ln_p_ref) + 1, 1), n_p_ref - 1) + 1)
    @inbounds fpress = (ln_p_ref[jpress - 1] - log_p_lay) / Δ_ln_p_ref
    return (jpress + tropo - 1, fpress)
end

"""
    compute_interp_frac_η(
        lkp::AbstractLookUp,
        vmr,
        tropo,
        jtemp,
        ibnd,
        glay,
        gcol,
    )

Compute interpolation fraction for binary species parameter.
"""
@inline function compute_interp_frac_η(n_η, ig, vmr_ref, (vmr1, vmr2), tropo, jtemp)
    FT = eltype(vmr1)
    itemp = 1
    @inbounds η_half = vmr_ref[tropo, ig[1] + 1, jtemp + itemp - 1] / vmr_ref[tropo, ig[2] + 1, jtemp + itemp - 1]
    col_mix1 = vmr1 + η_half * vmr2
    η = vmr1 * (FT(1) / col_mix1) # rewritten to ease register pressure
    if col_mix1 ≤ FT(0)
        η = FT(0.5)
    end
    #η = col_mix1 > 0 ? vmr1 / col_mix1 : FT(0.5) # rte-rrtmgp uses col_mix1 > tiny(col_mix1)
    loc_η = FT(η * (n_η - 1))
    jη1 = min(unsafe_trunc(Int, loc_η) + 1, n_η - 1)
    fη1 = loc_η - unsafe_trunc(Int, loc_η)

    itemp = 2
    @inbounds η_half = vmr_ref[tropo, ig[1] + 1, jtemp + itemp - 1] / vmr_ref[tropo, ig[2] + 1, jtemp + itemp - 1]
    col_mix2 = vmr1 + η_half * vmr2
    η = vmr1 * (FT(1) / col_mix2) # rewritten to ease register pressure
    if col_mix2 ≤ FT(0)
        η = FT(0.5)
    end
    #η = col_mix2 > 0 ? vmr1 / col_mix2 : FT(0.5) # rte-rrtmgp uses col_mix2 > tiny(col_mix2)
    loc_η = FT(η * (n_η - 1))
    jη2 = min(unsafe_trunc(Int, loc_η) + 1, n_η - 1)
    fη2 = loc_η - unsafe_trunc(Int, loc_η)

    return ((jη1, jη2, fη1, fη2), (col_mix1, col_mix2))
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
    ) where {FT<:AbstractFloat}

Compute optical thickness, single scattering albedo, and asymmetry parameter.
"""
@inline function compute_gas_optics(lkp::LookUpLW, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
    # upper/lower troposphere
    tropo = p_lay > lkp.p_ref_tropo ? 1 : 2
    # volume mixing ratio of h2o
    vmr_h2o = get_vmr(vmr, lkp.idx_h2o, glay, gcol)

    jftemp = compute_interp_frac_temp(lkp.ref_points.t_ref, t_lay)
    jfpress = compute_interp_frac_press(lkp.ref_points.ln_p_ref, p_lay, tropo)

    @inbounds kmajor = view(lkp.kmajor, :, :, :, igpt)
    n_η = size(kmajor, 1)

    ig = view(lkp.key_species, 1:2, tropo, ibnd)
    @inbounds vmr1 = get_vmr(vmr, ig[1], glay, gcol)
    @inbounds vmr2 = get_vmr(vmr, ig[2], glay, gcol)

    jfη, col_mix = compute_interp_frac_η(n_η, ig, lkp.ref_points.vmr_ref, (vmr1, vmr2), tropo, jftemp[1])
    # compute Planck fraction
    @inbounds planck_fraction = view(lkp.planck.planck_fraction, :, :, :, igpt)
    pfrac = interp3d(jfη..., jftemp..., jfpress..., planck_fraction)

    # computing τ_major
    τ_major = interp3d(jfη..., jftemp..., jfpress..., kmajor, col_mix...) * col_dry
    # computing τ_minor
    lkp_minor = tropo == 1 ? lkp.minor_lower : lkp.minor_upper
    τ_minor = compute_τ_minor(lkp_minor, vmr, vmr_h2o, col_dry, p_lay, t_lay, jftemp..., jfη..., igpt, ibnd, glay, gcol)
    # τ_Rayleigh is zero for longwave
    τ = max(τ_major + τ_minor, zero(τ_major))
    return (τ, zero(τ_major), zero(τ_major), pfrac) # initializing asymmetry parameter
end

@inline function compute_gas_optics(lkp::LookUpSW, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
    # upper/lower troposphere
    tropo = p_lay > lkp.p_ref_tropo ? 1 : 2
    # volume mixing ratio of h2o
    vmr_h2o = get_vmr(vmr, lkp.idx_h2o, glay, gcol)

    jftemp = compute_interp_frac_temp(lkp.ref_points.t_ref, t_lay)
    jfpress = compute_interp_frac_press(lkp.ref_points.ln_p_ref, p_lay, tropo)

    @inbounds kmajor = view(lkp.kmajor, :, :, :, igpt)
    n_η = size(kmajor, 1)

    ig = view(lkp.key_species, 1:2, tropo, ibnd)
    @inbounds vmr1 = get_vmr(vmr, ig[1], glay, gcol)
    @inbounds vmr2 = get_vmr(vmr, ig[2], glay, gcol)

    jfη, col_mix = compute_interp_frac_η(n_η, ig, lkp.ref_points.vmr_ref, (vmr1, vmr2), tropo, jftemp[1])

    # computing τ_major
    τ_major = interp3d(jfη..., jftemp..., jfpress..., kmajor, col_mix...) * col_dry
    # computing τ_minor
    lkp_minor = tropo == 1 ? lkp.minor_lower : lkp.minor_upper
    τ_minor = compute_τ_minor(lkp_minor, vmr, vmr_h2o, col_dry, p_lay, t_lay, jftemp..., jfη..., igpt, ibnd, glay, gcol)
    # compute τ_Rayleigh
    @inbounds rayleigh_coeff = tropo == 1 ? view(lkp.rayl_lower, :, :, igpt) : view(lkp.rayl_upper, :, :, igpt)
    τ_ray = compute_τ_rayleigh(rayleigh_coeff, col_dry, vmr_h2o, jftemp..., jfη..., igpt)
    FT = eltype(τ_major)
    τ = max(τ_major + τ_minor + τ_ray, FT(0))
    ssa = τ_ray * (FT(1) / τ)# rewritten to ease register pressure
    if τ ≤ FT(0)
        ssa = FT(0)
    end
    #ssa = τ > 0 ? τ_ray / τ : zero(τ) # single scattering albedo
    return (τ, ssa, zero(τ)) # initializing asymmetry parameter
end

"""
    compute_τ_minor(
        lkp::AbstractLookUp,
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
    lkp_minor::LookUpMinor,
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
) where {FT}
    (; kminor) = lkp_minor
    τ_minor = FT(0)
    (st_bnd, st_gpt, n) = LookUpTables.getbounds(lkp_minor, ibnd, igpt)

    if n > 0
        @inbounds begin
            pa2hpa = FT(0.01) # pascals to hectopascals
            dry_fact = FT(1) / (FT(1) + vmr_h2o)
            density_fact = pa2hpa * p_lay / t_lay

            for i in 0:(n - 1)
                idx_gas, idx_scaling_gas, scales_with_density, scale_by_complement =
                    LookUpTables.get_minor_gas_data(lkp_minor, st_bnd + i)

                vmr_imnr = get_vmr(vmr, idx_gas, glay, gcol)
                if vmr_imnr > 0
                    scaling = vmr_imnr * col_dry
                    if scales_with_density == 1
                        scaling *= density_fact
                        if idx_scaling_gas > 0
                            if scale_by_complement == 1
                                scaling *= (FT(1) - get_vmr(vmr, idx_scaling_gas, glay, gcol) * dry_fact)
                            else
                                scaling *= get_vmr(vmr, idx_scaling_gas, glay, gcol) * dry_fact
                            end
                        end
                    end
                    τ_minor += interp2d(fη1, fη2, ftemp, view(kminor, :, :, st_gpt + i), jη1, jη2, jtemp) * scaling
                end
            end
        end
    end
    return τ_minor
end

"""
    compute_τ_rayleigh(
        lkp::LookUpSW,
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
    igpt::Int,
) where {FT <: AbstractFloat} = interp2d(fη1, fη2, ftemp, rayleigh_coeff, jη1, jη2, jtemp) * (vmr_h2o + FT(1)) * col_dry
