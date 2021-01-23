
@kernel function compute_col_dry_kernel!(
    ::Val{nlay},
    ::Val{ncol},
    p_lev::FTA2D,
    t_lay::FTA2D,
    col_dry::FTA2D,
    mol_m_dry::FT,
    mol_m_h2o::FT,
    avogadro::FT,
    vmr_h2o::Union{FTA2D,Nothing} = nothing,
    lat::Union{FTA1D,Nothing} = nothing,
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    nlay,
    ncol,
}

    glay, gcol = @index(Global, NTuple)  # global col & lay ids

    # first and second term of Helmert formula
    helmert1 = FT(9.80665)
    helmert2 = FT(0.02586)

    if lat ≠ nothing
        g0 = helmert1 - helmert2 * cos(FT(2) * π * lat[gcol] / FT(180)) # acceleration due to gravity [m/s^2]
    else
        g0 = helmert1
    end

    Δp = p_lev[glay, gcol] - p_lev[glay+1, gcol]


    if vmr_h2o == nothing
        col_dry[glay, gcol] =
            (FT(10) * Δp * avogadro / (FT(1000) * mol_m_dry * FT(100) .* g0))
    else
        # Get average mass of moist air per mole of moist air
        m_air =
            (mol_m_dry + mol_m_h2o * vmr_h2o[glay, gcol]) /
            (1 + vmr_h2o[glay, gcol])
        # Hydrostatic equation
        col_dry[glay, gcol] =
            (FT(10) * Δp * avogadro / (FT(1000) * m_air * FT(100) .* g0)) /
            (FT(1) + vmr_h2o[glay, gcol])
    end

end

@kernel function compute_optical_props_kernel!(
    ::Val{nlay},
    ::Val{ncol},
    as::AtmosphericState{FT,FTA1D,FTA2D,VMR,I},
    lkp::AbstractLookUp{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D},
    op::AbstractOpticalProps{FT,FTA2D},
    igpt::Int,
    islw::Bool,
    sf::Union{Nothing,AbstractSource},
) where {
    I<:Int,
    FT<:AbstractFloat,
    UI8<:UInt8,
    UI8A1D<:AbstractArray{UI8,1},
    IA1D<:AbstractArray{I,1},
    IA2D<:AbstractArray{I,2},
    IA3D<:AbstractArray{I,3},
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    FTA4D<:AbstractArray{FT,4},
    VMR<:Vmr{FT,FTA1D,FTA2D},
    nlay,
    ncol,
}

    glay, gcol = @index(Global, NTuple)  # global col & lay ids

    fmajor = @private FT (2, 2, 2) # (η,p,T)
    fminor = @private FT (2, 2)    # (η,T)
    col_mix = @private FT (2)
    jη = @private I (2)
    ig = @private I (2)
    τ_maj = @private FT (1)
    τ_min = @private FT (1)
    τ_ray = @private FT (1)
    loc_in_bnd = @private I (1)

    # setting references
    Δlnp = lkp.Δ_ln_p_ref
    Δ_t_ref = lkp.Δ_t_ref
    Δ_ln_p_ref = lkp.Δ_ln_p_ref
    n_t_ref = lkp.n_t_ref
    t_ref = lkp.t_ref
    p_ref = lkp.p_ref
    key_species = lkp.key_species
    p_ref_tropo = lkp.p_ref_tropo
    kmajor = lkp.kmajor
    major_gpt2bnd = lkp.major_gpt2bnd
    minor_lower_bnd_st = lkp.minor_lower_bnd_st
    minor_upper_bnd_st = lkp.minor_upper_bnd_st
    idx_gases_minor_lower = lkp.idx_gases_minor_lower
    idx_gases_minor_upper = lkp.idx_gases_minor_upper
    bnd_lims_gpt = lkp.bnd_lims_gpt
    minor_lower_scales_with_density = lkp.minor_lower_scales_with_density
    minor_upper_scales_with_density = lkp.minor_upper_scales_with_density
    idx_scaling_gas_lower = lkp.idx_scaling_gas_lower
    idx_scaling_gas_upper = lkp.idx_scaling_gas_upper
    lower_scale_by_complement = lkp.lower_scale_by_complement
    upper_scale_by_complement = lkp.upper_scale_by_complement
    kminor_lower = lkp.kminor_lower
    kminor_upper = lkp.kminor_upper
    minor_lower_gpt_sh = lkp.minor_lower_gpt_sh
    minor_upper_gpt_sh = lkp.minor_upper_gpt_sh
    vmr_ref = lkp.vmr_ref
    idx_h2o = lkp.idx_h2o
    t_lay = as.t_lay
    p_lay = as.p_lay
    t_lev = as.t_lev
    t_sfc = as.t_sfc
    col_dry = as.col_dry
    vmr = as.vmr
    #-----------------------
    pa2hpa = FT(1) / FT(100) # pascals to hectopascals
    npres = length(p_ref)
    n_η = size(kmajor, 2)
    bnd = major_gpt2bnd[igpt]

    tropo = p_lay[glay, gcol] > p_ref_tropo ? 1 : 2

    jtemp = max(
        min(
            unsafe_trunc(I, (t_lay[glay, gcol] - t_ref[1]) / Δ_t_ref) + 1,
            n_t_ref - 1,
        ),
        1,
    )
    ftemp = (t_lay[glay, gcol] - t_ref[jtemp]) / Δ_t_ref

    jpress =
        min(
            max(
                unsafe_trunc(
                    I,
                    (log(p_ref[1]) - log(p_lay[glay, gcol])) / Δ_ln_p_ref,
                ) + 1,
                1,
            ),
            npres - 1,
        ) + 1
    fpress = (log(p_ref[jpress-1]) - log(p_lay[glay, gcol])) / Δ_ln_p_ref

    ig[1] = key_species[1, tropo, bnd]
    ig[2] = key_species[2, tropo, bnd]

    τ_maj[1] = FT(0)
    τ_min[1] = FT(0)
    τ_ray[1] = FT(0)

    if ig[1] ≠ 0 || ig[2] ≠ 0
        # τ_major
        for itemp = 1:2
            η_half =
                vmr_ref[tropo, ig[1]+1, jtemp+itemp-1] /
                vmr_ref[tropo, ig[2]+1, jtemp+itemp-1]
            col_mix[itemp] =
                get_vmr(vmr, ig[1], glay, gcol) +
                η_half * get_vmr(vmr, ig[2], glay, gcol)
            if col_mix[itemp] ≥ eps(FT) * 2
                η = get_vmr(vmr, ig[1], glay, gcol) / col_mix[itemp]
            else
                η = FT(0.5)
            end
            loc_η = FT(η * (n_η - 1))
            jη[itemp] = min(unsafe_trunc(I, loc_η) + 1, n_η - 1)
            fη = loc_η % FT(1)

            # compute interpolation fractions needed for minor species
            # ftemp_term = (1._wp-ftemp(icol,ilay)) for itemp = 1, ftemp(icol,ilay) for itemp=2
            if itemp == 1
                ftemp_term = FT(1) - ftemp
            else
                ftemp_term = ftemp
            end
            fminor[1, itemp] = (FT(1) - fη) * ftemp_term
            fminor[2, itemp] = fη * ftemp_term
            # compute interpolation fractions needed for major species
            fmajor[1, 1, itemp] = (FT(1) - fpress) * fminor[1, itemp]
            fmajor[2, 1, itemp] = (FT(1) - fpress) * fminor[2, itemp]
            fmajor[1, 2, itemp] = fpress * fminor[1, itemp]
            fmajor[2, 2, itemp] = fpress * fminor[2, itemp]
        end
        # computing τ_major
        jpresst = jpress + tropo - 1
        τ_maj[1] =
            (
                col_mix[1] * (
                    fmajor[1, 1, 1] * kmajor[igpt, jη[1], jpresst-1, jtemp] +
                    fmajor[2, 1, 1] * kmajor[igpt, jη[1]+1, jpresst-1, jtemp] +
                    fmajor[1, 2, 1] * kmajor[igpt, jη[1], jpresst, jtemp] +
                    fmajor[2, 2, 1] * kmajor[igpt, jη[1]+1, jpresst, jtemp]
                ) +
                col_mix[2] * (
                    fmajor[1, 1, 2] * kmajor[igpt, jη[2], jpresst-1, jtemp+1] +
                    fmajor[2, 1, 2] *
                    kmajor[igpt, jη[2]+1, jpresst-1, jtemp+1] +
                    fmajor[1, 2, 2] * kmajor[igpt, jη[2], jpresst, jtemp+1] +
                    fmajor[2, 2, 2] * kmajor[igpt, jη[2]+1, jpresst, jtemp+1]
                )
            ) * col_dry[glay, gcol]
        #-----------------------------------------------------------------------------------
        # τ_minor
        loc_in_bnd[1] = igpt - (bnd_lims_gpt[1, bnd] - 1)

        τ_min[1] = FT(0)
        if tropo == 1 # lower
            for i = minor_lower_bnd_st[bnd]:minor_lower_bnd_st[bnd+1]-1
                imnr = idx_gases_minor_lower[i]
                vmr_imnr = get_vmr(vmr, imnr, glay, gcol)
                if vmr_imnr > eps(FT) * 2
                    scaling =
                        get_vmr(vmr, imnr, glay, gcol) * col_dry[glay, gcol]

                    if minor_lower_scales_with_density[i] == 1
                        scaling *=
                            (pa2hpa * p_lay[glay, gcol] / t_lay[glay, gcol])
                        sgas = idx_scaling_gas_lower[i]
                        if sgas > 0
                            dry_fact =
                                FT(1) /
                                (FT(1) + get_vmr(vmr, idx_h2o, glay, gcol))
                            if lower_scale_by_complement[i] == 1
                                scaling *= (
                                    FT(1) -
                                    get_vmr(vmr, sgas, glay, gcol) * dry_fact
                                )
                            else
                                scaling *=
                                    get_vmr(vmr, sgas, glay, gcol) * dry_fact
                            end
                        end
                    end
                    k_loc = minor_lower_gpt_sh[i] + loc_in_bnd[1]
                    τ_min[1] +=
                        (
                            fminor[1, 1] * kminor_lower[k_loc, jη[1], jtemp] +
                            fminor[2, 1] * kminor_lower[k_loc, jη[1]+1, jtemp] +
                            fminor[1, 2] * kminor_lower[k_loc, jη[2], jtemp+1] +
                            fminor[2, 2] *
                            kminor_lower[k_loc, jη[2]+1, jtemp+1]
                        ) * scaling
                end
            end
        else # upper
            for i = minor_upper_bnd_st[bnd]:minor_upper_bnd_st[bnd+1]-1
                imnr = idx_gases_minor_upper[i]
                vmr_imnr = get_vmr(vmr, imnr, glay, gcol)
                if vmr_imnr > eps(FT) * 2
                    scaling = vmr_imnr * col_dry[glay, gcol]

                    if minor_upper_scales_with_density[i] == 1
                        scaling *=
                            (pa2hpa * p_lay[glay, gcol] / t_lay[glay, gcol])
                        sgas = idx_scaling_gas_upper[i]
                        if sgas > 0
                            dry_fact =
                                FT(1) /
                                (FT(1) + get_vmr(vmr, idx_h2o, glay, gcol))
                            if upper_scale_by_complement[i] == 1
                                scaling *= (
                                    FT(1) -
                                    get_vmr(vmr, sgas, glay, gcol) * dry_fact
                                )
                            else
                                scaling *=
                                    get_vmr(vmr, sgas, glay, gcol) * dry_fact
                            end
                        end
                    end
                    k_loc = minor_upper_gpt_sh[i] + loc_in_bnd[1]
                    τ_min[1] +=
                        (
                            fminor[1, 1] * kminor_upper[k_loc, jη[1], jtemp] +
                            fminor[2, 1] * kminor_upper[k_loc, jη[1]+1, jtemp] +
                            fminor[1, 2] * kminor_upper[k_loc, jη[2], jtemp+1] +
                            fminor[2, 2] *
                            kminor_upper[k_loc, jη[2]+1, jtemp+1]
                        ) * scaling
                end
            end
        end

        # τ_Rayleigh
        if hasfield(typeof(lkp), :rayl_lower)
            if tropo == 1 # lower
                krayl = lkp.rayl_lower
            else # upper
                krayl = lkp.rayl_upper
            end
            scaling =
                (get_vmr(vmr, idx_h2o, glay, gcol) + FT(1)) *
                col_dry[glay, gcol]
            τ_ray[1] +=
                (
                    fminor[1, 1] * kminor_upper[igpt, jη[1], jtemp] +
                    fminor[2, 1] * kminor_upper[igpt, jη[1]+1, jtemp] +
                    fminor[1, 2] * kminor_upper[igpt, jη[2], jtemp+1] +
                    fminor[2, 2] * kminor_upper[igpt, jη[2]+1, jtemp+1]
                ) * scaling
        end
        #-----------------------------------------------------------------------------
    end

    op.τ[glay, gcol] = τ_maj[1] + τ_min[1] + τ_ray[1]

    # computing Planck sources for longwave problem
    if hasfield(typeof(lkp), :totplnk) && islw # compute source functions
        p_frac_ref = lkp.planck_fraction
        t_plnk_ref = lkp.t_planck
        n_t_plnk = lkp.n_t_plnk
        totplnk = lkp.totplnk
        lay_source = sf.lay_source
        lev_source_inc = sf.lev_source_inc
        lev_source_dec = sf.lev_source_dec
        sfc_source = sf.sfc_source
        t_plnk_min = t_plnk_ref[1]
        t_plnk_Δt = t_plnk_ref[2] - t_plnk_ref[1]
        if ig[1] ≠ 0 || ig[2] ≠ 0
            # compute Planck fraction
            p_frac =
                fmajor[1, 1, 1] * p_frac_ref[igpt, jη[1], jpresst-1, jtemp] +
                fmajor[2, 1, 1] * p_frac_ref[igpt, jη[1]+1, jpresst-1, jtemp] +
                fmajor[1, 2, 1] * p_frac_ref[igpt, jη[1], jpresst, jtemp] +
                fmajor[2, 2, 1] * p_frac_ref[igpt, jη[1]+1, jpresst, jtemp] +
                fmajor[1, 1, 2] * p_frac_ref[igpt, jη[2], jpresst-1, jtemp+1] +
                fmajor[2, 1, 2] *
                p_frac_ref[igpt, jη[2]+1, jpresst-1, jtemp+1] +
                fmajor[1, 2, 2] * p_frac_ref[igpt, jη[2], jpresst, jtemp+1] +
                fmajor[2, 2, 2] * p_frac_ref[igpt, jη[2]+1, jpresst, jtemp+1]

            # computing lay_source
            temp = t_lay[glay, gcol]
            jtemp_pl = max(
                min(
                    unsafe_trunc(I, (temp - t_plnk_min) / t_plnk_Δt) + 1,
                    n_t_plnk - 1,
                ),
                1,
            )
            ftemp_pl = (temp - t_plnk_ref[jtemp_pl]) / t_plnk_Δt
            lay_source[glay, gcol] =
                (
                    totplnk[jtemp_pl, bnd] * (FT(1) - ftemp_pl) +
                    totplnk[jtemp_pl+1, bnd] * ftemp_pl
                ) * p_frac

            # computing lev_source_inc
            temp = t_lev[glay+1, gcol]
            jtemp_pl = max(
                min(
                    unsafe_trunc(I, (temp - t_plnk_min) / t_plnk_Δt) + 1,
                    n_t_plnk - 1,
                ),
                1,
            )
            ftemp_pl = (temp - t_plnk_ref[jtemp_pl]) / t_plnk_Δt
            lev_source_inc[glay, gcol] =
                (
                    totplnk[jtemp_pl, bnd] * (FT(1) - ftemp_pl) +
                    totplnk[jtemp_pl+1, bnd] * ftemp_pl
                ) * p_frac

            # computing lev_source_dec
            temp = t_lev[glay, gcol]
            jtemp_pl = max(
                min(
                    unsafe_trunc(I, (temp - t_plnk_min) / t_plnk_Δt) + 1,
                    n_t_plnk - 1,
                ),
                1,
            )
            ftemp_pl = (temp - t_plnk_ref[jtemp_pl]) / t_plnk_Δt
            lev_source_dec[glay, gcol] =
                (
                    totplnk[jtemp_pl, bnd] * (FT(1) - ftemp_pl) +
                    totplnk[jtemp_pl+1, bnd] * ftemp_pl
                ) * p_frac

            # computing sfc_source 
            if glay == 1
                temp = t_sfc[gcol]
                jtemp_pl = max(
                    min(
                        unsafe_trunc(I, (temp - t_plnk_min) / t_plnk_Δt) + 1,
                        n_t_plnk - 1,
                    ),
                    1,
                )
                ftemp_pl = (temp - t_plnk_ref[jtemp_pl]) / t_plnk_Δt
                sfc_source[gcol] =
                    (
                        totplnk[jtemp_pl, bnd] * (FT(1) - ftemp_pl) +
                        totplnk[jtemp_pl+1, bnd] * ftemp_pl
                    ) * p_frac
            end
        else
            lay_source[glay, gcol] = FT(0)
            lev_source_inc[glay, gcol] = FT(0)
            lev_source_dec[glay, gcol] = FT(0)
            if glay == 1
                sfc_source[gcol] = FT(0)
            end
        end
    end

    #@cuprintln("τ_major = $(τ_maj[1]); τ_minor = $(τ_min[1]); glay = $glay; gcol = $gcol;")
    @synchronize

end

# This functions calculates the optical thickness based on pressure 
# and lapse rate for a gray atmosphere. 
# See Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
# https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2

@kernel function compute_optical_props_kernel!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I},
    τ::FTA2D,
    islw::Bool,
    source::Union{AbstractSource{FT,FTA2D},Nothing},
    ::Val{nlay},
    ::Val{ncol},
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    nlay,
    ncol,
}
    glay, gcol = @index(Global, NTuple)  # global col & lay ids
    p_lay = as.p_lay
    p_lev = as.p_lev
    d0 = as.d0
    α = as.α

    p0 = p_lev[1, gcol]

    τ[glay, gcol] = abs(
        (α * d0[gcol] * (p_lay[glay, gcol] ./ p0) .^ α ./ p_lay[glay, gcol]) * (p_lev[glay+1, gcol] - p_lev[glay, gcol]),
    )

    if source ≠ nothing && islw # computing sources for longwave problem
        sbc = FT(Stefan())
        source.lay_source[glay, gcol] = sbc * as.t_lay[glay, gcol]^FT(4) / FT(π)   # computing lay_source
        source.lev_source_inc[glay, gcol] =
            sbc * as.t_lev[glay+1, gcol]^FT(4) / FT(π)
        source.lev_source_dec[glay, gcol] =
            sbc * as.t_lev[glay, gcol]^FT(4) / FT(π)
        if glay == 1
            source.sfc_source[gcol] = sbc * as.t_sfc[gcol]^FT(4) / FT(π)   # computing sfc_source
        end
    end

    @synchronize
end
