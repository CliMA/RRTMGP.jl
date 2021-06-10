using UnPack

@kernel function compute_col_dry_kernel!(
    ::Val{nlay},
    ::Val{ncol},
    p_lev::FTA2D,
    t_lay::FTA2D,
    col_dry::FTA2D,
    mol_m_dry::FT,
    mol_m_h2o::FT,
    avogadro::FT,
    helmert1::FT,
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
    helmert2 = FT(0.02586)  # second constant of Helmert formula
    m2_to_cm2 = FT(100 * 100) # m^2 to cm^2
    if lat isa AbstractArray
        g0 = helmert1 - helmert2 * cos(FT(2) * FT(π) * lat[gcol] / FT(180)) # acceleration due to gravity [m/s^2]
    else
        g0 = helmert1
    end

    Δp = p_lev[glay, gcol] - p_lev[glay+1, gcol]

    if vmr_h2o isa AbstractArray
        fact = FT(1) / (FT(1) + vmr_h2o[glay, gcol])
        # Get average mass of moist air per mole of moist air
        m_air = (mol_m_dry + mol_m_h2o * vmr_h2o[glay, gcol]) * fact
        # Hydrostatic equation
        col_dry[glay, gcol] = (Δp * avogadro / (m2_to_cm2 * m_air * g0)) * fact
    else
        col_dry[glay, gcol] = (Δp * avogadro / (m2_to_cm2 * mol_m_dry * g0)) # molecules/cm^2
    end
end

@kernel function compute_optical_props_kernel!(
    ::Val{nlay},
    ::Val{ncol},
    as::AtmosphericState{FT,FTA1D,FTA1DN,FTA2D,CLDP,CLDM,VMR,I},
    lkp::AbstractLookUp{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D},
    op::AbstractOpticalProps{FT,FTA2D},
    igpt::Int,
    islw::Bool,
    sf::Union{Nothing,AbstractSource{FT}},
    lkp_cld::Union{LookUpCld{I,B,FT,FTA1D,FTA2D,FTA3D,FTA4D},Nothing},
) where {
    I<:Int,
    B<:Bool,
    FT<:AbstractFloat,
    UI8<:UInt8,
    UI8A1D<:AbstractArray{UI8,1},
    IA1D<:AbstractArray{I,1},
    IA2D<:AbstractArray{I,2},
    IA3D<:AbstractArray{I,3},
    FTA1D<:AbstractArray{FT,1},
    FTA1DN<:Union{AbstractArray{FT,1},Nothing},
    FTA2D<:AbstractArray{FT,2},
    CLDP<:Union{AbstractArray{FT,2},Nothing},
    CLDM<:Union{AbstractArray{Bool,2},Nothing},
    FTA3D<:AbstractArray{FT,3},
    FTA4D<:AbstractArray{FT,4},
    VMR<:AbstractVmr{FT},
    nlay,
    ncol,
}

    glay, gcol = @index(Global, NTuple)  # global col & lay ids

    fmajor = @private FT (2, 2, 2) # (η,p,T)
    fminor = @private FT (2, 2)    # (η,T)
    col_mix = @private FT (2)
    scaling_plnck = @private FT (2)
    jη = @private I (2)
    ig = @private I (2)
    τ_maj = @private FT (1)
    τ_min = @private FT (1)
    τ_ray = @private FT (1)
    loc_in_bnd = @private I (1)

    # setting references
    @unpack Δ_t_ref, Δ_ln_p_ref, n_t_ref, t_ref, p_ref = lkp
    @unpack key_species, p_ref_tropo, kmajor, major_gpt2bnd = lkp
    @unpack bnd_lims_gpt, vmr_ref, idx_h2o = lkp
    @unpack t_lay, p_lay, t_lev, t_sfc, col_dry, vmr = as

    tropo = p_lay[glay, gcol] > p_ref_tropo ? 1 : 2

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


    pa2hpa = FT(1) / FT(100) # pascals to hectopascals
    npres = length(p_ref)
    n_η = size(kmajor, 2)
    bnd = major_gpt2bnd[igpt]

    jtemp = loc_lower(t_lay[glay, gcol], Δ_t_ref, n_t_ref, t_ref)
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

    if ig[1] == 0 && ig[2] == 0
        ig .= 2
    end

    τ_maj[1] = τ_min[1] = τ_ray[1] = FT(0)

    # compute fminor and fmajor
    @inbounds for itemp = 1:2
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
        interp3d(col_mix, fmajor, kmajor, igpt, jη, jpresst, jtemp) *
        col_dry[glay, gcol]
    #-----------------------------------------------------------------------------------
    # computing τ_minor
    loc_in_bnd[1] = igpt - (bnd_lims_gpt[1, bnd] - 1)

    τ_min[1] = FT(0)

    @inbounds for i = minor_bnd_st[bnd]:minor_bnd_st[bnd+1]-1
        imnr = idx_gases_minor[i]
        vmr_imnr = get_vmr(vmr, imnr, glay, gcol)
        if vmr_imnr > eps(FT) * 2
            scaling = vmr_imnr * col_dry[glay, gcol]

            if minor_scales_with_density[i] == 1
                scaling *= (pa2hpa * p_lay[glay, gcol] / t_lay[glay, gcol])
                sgas = idx_scaling_gas[i]
                if sgas > 0
                    dry_fact =
                        FT(1) / (FT(1) + get_vmr(vmr, idx_h2o, glay, gcol))
                    if scale_by_complement[i] == 1
                        scaling *=
                            (FT(1) - get_vmr(vmr, sgas, glay, gcol) * dry_fact)
                    else
                        scaling *= get_vmr(vmr, sgas, glay, gcol) * dry_fact
                    end
                end
            end
            k_loc = minor_gpt_sh[i] + loc_in_bnd[1]
            τ_min[1] += interp2d(fminor, kminor, k_loc, jη, jtemp) * scaling
        end
    end

    # compute τ_Rayleigh
    if hasfield(typeof(lkp), :rayl_lower)
        if tropo == 1 # lower
            krayl = lkp.rayl_lower
        else # upper
            krayl = lkp.rayl_upper
        end
        scaling =
            (get_vmr(vmr, idx_h2o, glay, gcol) + FT(1)) * col_dry[glay, gcol]
        τ_ray[1] += interp2d(fminor, krayl, igpt, jη, jtemp) * scaling
    end
    #-----------------------------------------------------------------------------

    op.τ[glay, gcol] = τ_maj[1] + τ_min[1] + τ_ray[1]

    if hasfield(typeof(op), :ssa) # single scattering albedo
        if op.τ[glay, gcol] > 2 * eps(FT)
            op.ssa[glay, gcol] = τ_ray[1] / op.τ[glay, gcol]
        else
            op.ssa[glay, gcol] = FT(0)
        end
    end

    if hasfield(typeof(op), :g) # initializing asymmetry parameter
        op.g[glay, gcol] = FT(0)
    end

    # computing Planck sources for longwave problem
    if hasfield(typeof(lkp), :totplnk) && islw # compute source functions
        @unpack planck_fraction, t_planck, n_t_plnk, totplnk = lkp
        @unpack lay_source, lev_source_inc, lev_source_dec, sfc_source = sf
        # compute Planck fraction
        scaling_plnck .= FT(1)
        p_frac = interp3d(
            scaling_plnck,
            fmajor,
            planck_fraction,
            igpt,
            jη,
            jpresst,
            jtemp,
        )
        planck_args = (t_planck, totplnk, bnd)
        # computing lay_source
        lay_source[glay, gcol] =
            interp1d(t_lay[glay, gcol], planck_args...) * p_frac
        # computing lev_source_inc
        lev_source_inc[glay, gcol] =
            interp1d(t_lev[glay+1, gcol], planck_args...) * p_frac
        # computing lev_source_dec
        lev_source_dec[glay, gcol] =
            interp1d(t_lev[glay, gcol], planck_args...) * p_frac
        if glay == 1 # computing sfc_source
            sfc_source[gcol] = interp1d(t_sfc[gcol], planck_args...) * p_frac
        end
        #-----------------------------------------------------------------------------
    end

    # cloud optics----------------------------------------------------------------
    if lkp_cld isa LookUpCld &&
       as.cld_mask isa AbstractArray &&
       op isa TwoStream
        ibnd = I(major_gpt2bnd[igpt])
        τ_cl, τ_cl_ssa, τ_cl_ssag =
            compute_cld_props(lkp_cld, as, gcol, glay, ibnd, igpt)

        if islw == false
            τ_cl, τ_cl_ssa, τ_cl_ssag = delta_scale(τ_cl, τ_cl_ssa, τ_cl_ssag)

        end
        increment!(op, τ_cl, τ_cl_ssa, τ_cl_ssag, glay, gcol, igpt)
    end
    #-----------------------------------------------------------------------------
    #@cuprintln("τ_major = $(τ_maj[1]); τ_minor = $(τ_min[1]); glay = $glay; gcol = $gcol;")
end

loc_lower(xi, Δx, n, x) =
    max(min(unsafe_trunc(Int, (xi - x[1]) / Δx) + 1, n - 1), 1)

function interp1d(xi::FT, x, y, col) where {FT<:AbstractFloat}
    Δx = x[2] - x[1]
    n = length(x)
    loc = loc_lower(xi, Δx, n, x)
    factor = (xi - x[loc]) / Δx
    return (y[loc, col] * (FT(1) - factor) + y[loc+1, col] * factor)
end

interp2d(fminor, coeff, igpt, jη, jtemp) =
    fminor[1, 1] * coeff[igpt, jη[1], jtemp] +
    fminor[2, 1] * coeff[igpt, jη[1]+1, jtemp] +
    fminor[1, 2] * coeff[igpt, jη[2], jtemp+1] +
    fminor[2, 2] * coeff[igpt, jη[2]+1, jtemp+1]

interp3d(scaling, fmajor, coeff, igpt, jη, jpresst, jtemp) =
    scaling[1] * (
        fmajor[1, 1, 1] * coeff[igpt, jη[1], jpresst-1, jtemp] +
        fmajor[2, 1, 1] * coeff[igpt, jη[1]+1, jpresst-1, jtemp] +
        fmajor[1, 2, 1] * coeff[igpt, jη[1], jpresst, jtemp] +
        fmajor[2, 2, 1] * coeff[igpt, jη[1]+1, jpresst, jtemp]
    ) +
    scaling[2] * (
        fmajor[1, 1, 2] * coeff[igpt, jη[2], jpresst-1, jtemp+1] +
        fmajor[2, 1, 2] * coeff[igpt, jη[2]+1, jpresst-1, jtemp+1] +
        fmajor[1, 2, 2] * coeff[igpt, jη[2], jpresst, jtemp+1] +
        fmajor[2, 2, 2] * coeff[igpt, jη[2]+1, jpresst, jtemp+1]
    )

function compute_cld_props(lkp_cld, as, gcol, glay, ibnd, igpt)
    use_lut = lkp_cld.use_lut
    cld_mask = as.cld_mask[glay, gcol]
    re_liq = as.cld_r_eff_liq[glay, gcol]
    re_ice = as.cld_r_eff_ice[glay, gcol]
    ice_rgh = as.ice_rgh
    cld_path_liq = as.cld_path_liq[glay, gcol]
    cld_path_ice = as.cld_path_ice[glay, gcol]
    FT = eltype(re_liq)
    τl, τl_ssa, τl_ssag = FT(0), FT(0), FT(0)
    τi, τi_ssa, τi_ssag = FT(0), FT(0), FT(0)
    if use_lut # use LUT interpolation
        @unpack lut_extliq,
        lut_ssaliq,
        lut_asyliq,
        lut_extice,
        lut_ssaice,
        lut_asyice,
        radliq_lwr,
        radliq_upr,
        radice_lwr,
        radice_upr,
        nsize_liq,
        nsize_ice = lkp_cld
        Δr_liq = (radliq_upr - radliq_lwr) / FT(nsize_liq - 1)
        Δr_ice = (radice_upr - radice_lwr) / FT(nsize_ice - 1)
        # cloud liquid particles
        if cld_path_liq > eps(FT)
            loc = Int(max(
                min(
                    unsafe_trunc(Int, (re_liq - radliq_lwr) / Δr_liq) + 1,
                    nsize_liq - 1,
                ),
                1,
            ))
            fac = (re_liq - radliq_lwr - (loc - 1) * Δr_liq) / Δr_liq
            fc1 = FT(1) - fac
            τl =
                (fc1 * lut_extliq[loc, ibnd] + fac * lut_extliq[loc+1, ibnd]) *
                cld_path_liq
            τl_ssa =
                (fc1 * lut_ssaliq[loc, ibnd] + fac * lut_ssaliq[loc+1, ibnd]) *
                τl
            τl_ssag =
                (fc1 * lut_asyliq[loc, ibnd] + fac * lut_asyliq[loc+1, ibnd]) *
                τl_ssa
        end
        # cloud ice particles
        if cld_path_ice > eps(FT)
            loc = Int(max(
                min(
                    unsafe_trunc(Int, (re_ice - radice_lwr) / Δr_ice) + 1,
                    nsize_ice - 1,
                ),
                1,
            ))
            fac = (re_ice - radice_lwr - (loc - 1) * Δr_ice) / Δr_ice
            fc1 = FT(1) - fac
            τi =
                (
                    fc1 * lut_extice[loc, ibnd, ice_rgh] +
                    fac * lut_extice[loc+1, ibnd, ice_rgh]
                ) * cld_path_ice
            τi_ssa =
                (
                    fc1 * lut_ssaice[loc, ibnd, ice_rgh] +
                    fac * lut_ssaice[loc+1, ibnd, ice_rgh]
                ) * τi
            τi_ssag =
                (
                    fc1 * lut_asyice[loc, ibnd, ice_rgh] +
                    fac * lut_asyice[loc+1, ibnd, ice_rgh]
                ) * τi_ssa
        end
    else # use pade interpolation
        @unpack pade_extliq,
        pade_ssaliq,
        pade_asyliq,
        pade_extice,
        pade_ssaice,
        pade_asyice,
        pade_sizreg_extliq,
        pade_sizreg_ssaliq,
        pade_sizreg_asyliq,
        pade_sizreg_extice,
        pade_sizreg_ssaice,
        pade_sizreg_asyice = lkp_cld
        m_ext, m_ssa_g = 3, 3
        n_ext, n_ssa_g = 3, 2
        # Finds index into size regime table
        # This works only if there are precisely three size regimes (four bounds) and it's
        # previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
        if cld_path_liq > eps(FT)
            irad = Int(min(
                floor(
                    (re_liq - pade_sizreg_extliq[2]) / pade_sizreg_extliq[3],
                ) + 2,
                3,
            ))
            τl =
                pade_eval(ibnd, re_liq, irad, m_ext, n_ext, pade_extliq) *
                cld_path_liq

            irad = Int(min(
                floor(
                    (re_liq - pade_sizreg_ssaliq[2]) / pade_sizreg_ssaliq[3],
                ) + 2,
                3,
            ))
            τl_ssa =
                (
                    FT(1) - max(
                        FT(0),
                        pade_eval(
                            ibnd,
                            re_liq,
                            irad,
                            m_ssa_g,
                            n_ssa_g,
                            pade_ssaliq,
                        ),
                    )
                ) * τl

            irad = Int(min(
                floor(
                    (re_liq - pade_sizreg_asyliq[2]) / pade_sizreg_asyliq[3],
                ) + 2,
                3,
            ))
            τl_ssag =
                pade_eval(ibnd, re_liq, irad, m_ssa_g, n_ssa_g, pade_asyliq) *
                τl_ssa
        end

        if cld_path_ice > eps(FT)
            irad = Int(min(
                floor(
                    (re_ice - pade_sizreg_extice[2]) / pade_sizreg_extice[3],
                ) + 2,
                3,
            ))

            τi =
                pade_eval(
                    ibnd,
                    re_ice,
                    irad,
                    m_ext,
                    n_ext,
                    pade_extice,
                    ice_rgh,
                ) * cld_path_ice

            irad = Int(min(
                floor(
                    (re_ice - pade_sizreg_ssaice[2]) / pade_sizreg_ssaice[3],
                ) + 2,
                3,
            ))
            τi_ssa =
                (
                    FT(1) - max(
                        FT(0),
                        pade_eval(
                            ibnd,
                            re_ice,
                            irad,
                            m_ssa_g,
                            n_ssa_g,
                            pade_ssaice,
                            ice_rgh,
                        ),
                    )
                ) * τi

            irad = Int(min(
                floor(
                    (re_ice - pade_sizreg_asyice[2]) / pade_sizreg_asyice[3],
                ) + 2,
                3,
            ))
            τi_ssag =
                pade_eval(
                    ibnd,
                    re_ice,
                    irad,
                    m_ssa_g,
                    n_ssa_g,
                    pade_asyice,
                    ice_rgh,
                ) * τi_ssa
        end
    end

    τ = τl + τi
    τ_ssa = τl_ssa + τi_ssa
    τ_ssag = (τl_ssag + τi_ssag) / max(eps(FT), τ_ssa)
    τ_ssa /= max(eps(FT), τ)

    return (τ, τ_ssa, τ_ssag)
end

function increment!(
    op::TwoStream{FT},
    τ2,
    ssa2,
    g2,
    glay,
    gcol,
    igpt,
) where {FT<:AbstractFloat}
    τ1, ssa1, g1 = op.τ[glay, gcol], op.ssa[glay, gcol], op.g[glay, gcol]
    τ = τ1 + τ2
    ssa = τ1 * ssa1 + τ2 * ssa2
    ssag = (τ1 * ssa1 * g1 + τ2 * ssa2 * g2) / max(eps(FT), ssa)
    ssa /= max(eps(FT), τ)

    op.τ[glay, gcol] = τ
    op.ssa[glay, gcol] = ssa
    op.g[glay, gcol] = ssag
    return nothing
end

function pade_eval(
    ibnd,
    re,
    irad,
    m,
    n,
    pade_coeffs,
    irgh::Union{Int,Nothing} = nothing,
)
    FT = eltype(re)
    if irgh isa Int
        coeffs = view(pade_coeffs, :, :, :, irgh)

    else
        coeffs = pade_coeffs
    end

    denom = coeffs[ibnd, irad, n+m]
    for i = (n+m-1):-1:(1+m)
        denom = coeffs[ibnd, irad, i] + re * denom
    end
    denom = FT(1) + re * denom

    numer = coeffs[ibnd, irad, m]
    for i = (m-1):-1:2
        numer = coeffs[ibnd, irad, i] + re * numer
    end
    numer = coeffs[ibnd, irad, 1] + re * numer

    return (numer / denom)
end

function delta_scale(τ, ssa, g)
    FT = typeof(τ)
    f = g * g
    wf = ssa * f
    τ_s = (FT(1) - wf) * τ
    ssa_s = (ssa - wf) / max(eps(FT), FT(1) - wf)
    g_s = (g - f) / max(eps(FT), FT(1) - f)
    return (τ_s, ssa_s, g_s)
end
# This functions calculates the optical thickness based on pressure 
# and lapse rate for a gray atmosphere. 
# See Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
# https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2

@kernel function compute_optical_props_kernel!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I},
    τ::FTA2D,
    islw::Bool,
    source::Union{AbstractSource{FT},Nothing},
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
    # setting references
    @unpack p_lay, p_lev, t_lay, t_lev, d0, α = as
    p0 = p_lev[1, gcol]

    τ[glay, gcol] = abs(
        (α * d0[gcol] * (p_lay[glay, gcol] ./ p0) .^ α ./ p_lay[glay, gcol]) * (p_lev[glay+1, gcol] - p_lev[glay, gcol]),
    )

    if source isa SourceLWNoScat || source isa SourceLW2Str # computing sources for longwave problem
        sbc = FT(Stefan())
        @unpack lay_source, lev_source_inc, lev_source_dec, sfc_source = source
        lay_source[glay, gcol] = sbc * t_lay[glay, gcol]^FT(4) / FT(π)   # computing lay_source
        lev_source_inc[glay, gcol] = sbc * t_lev[glay+1, gcol]^FT(4) / FT(π)
        lev_source_dec[glay, gcol] = sbc * t_lev[glay, gcol]^FT(4) / FT(π)
        if glay == 1
            sfc_source[gcol] = sbc * as.t_sfc[gcol]^FT(4) / FT(π)   # computing sfc_source
        end
    end
end
