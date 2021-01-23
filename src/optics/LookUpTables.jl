module LookUpTables

using DocStringExtensions
using ..Device: array_type
using Adapt

export AbstractLookUp, LookUpLW, LookUpSW

abstract type AbstractLookUp{
    I,
    FT,
    UI8,
    UI8A1D,
    IA1D,
    IA2D,
    IA3D,
    FTA1D,
    FTA2D,
    FTA3D,
    FTA4D,
} end

struct LookUpLW{
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
} <: AbstractLookUp{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D}
    n_gases::I
    n_bnd::I
    n_gpt::I
    n_atmos_layers::I
    n_t_ref::I
    n_p_ref::I
    n_t_plnk::I
    n_mix_frac::I
    n_maj_absrb::I
    n_min_absrb::I
    n_min_absrb_lower::I
    n_min_absrb_upper::I
    n_absrb_ext::I
    n_contrib_lower::I
    n_contrib_upper::I
    idx_h2o::I

    p_ref_tropo::FT
    t_ref_absrb::FT
    p_ref_absrb::FT

    Δ_t_ref::FT
    Δ_ln_p_ref::FT

    idx_gases_minor_lower::IA1D
    idx_gases_minor_upper::IA1D

    idx_scaling_gas_lower::IA1D
    idx_scaling_gas_upper::IA1D

    key_species::IA3D

    kmajor::FTA4D
    kminor_lower::FTA3D
    kminor_upper::FTA3D

    kminor_start_lower::IA1D
    kminor_start_upper::IA1D

    planck_fraction::FTA4D
    t_planck::FTA1D
    totplnk::FTA2D

    major_gpt2bnd::UI8A1D
    bnd_lims_gpt::IA2D      # TODO: Marked for removal
    bnd_lims_wn::FTA2D

    minor_lower_gpt_lims::IA2D
    minor_upper_gpt_lims::IA2D

    minor_lower_bnd::UI8A1D
    minor_upper_bnd::UI8A1D

    minor_lower_bnd_st::UI8A1D
    minor_upper_bnd_st::UI8A1D

    minor_lower_gpt_sh::IA1D
    minor_upper_gpt_sh::IA1D

    minor_lower_scales_with_density::IA1D
    minor_upper_scales_with_density::IA1D

    lower_scale_by_complement::IA1D
    upper_scale_by_complement::IA1D

    p_ref::FTA1D
    t_ref::FTA1D
    vmr_ref::FTA3D
end
Adapt.@adapt_structure LookUpLW

function LookUpLW(
    ds,
    ::Type{I},
    ::Type{FT},
    ::Type{DA},
) where {I<:Int,FT<:AbstractFloat,DA}

    UI8 = UInt8
    UI8A1D = DA{UInt8,1}
    STA = Array{String,1}
    IA1D = DA{I,1}
    IA2D = DA{I,2}
    IA3D = DA{I,3}
    FTA1D = DA{FT,1}
    FTA2D = DA{FT,2}
    FTA3D = DA{FT,3}
    FTA4D = DA{FT,4}
    DSTAI = Dict{String,I}

    n_bnd = I(ds.dim["bnd"])
    n_gpt = I(ds.dim["gpt"])
    n_atmos_layers = I(ds.dim["atmos_layer"])
    n_t_ref = I(ds.dim["temperature"])
    n_p_ref = I(ds.dim["pressure"])
    n_t_plnk = I(ds.dim["temperature_Planck"])
    n_mix_frac = I(ds.dim["mixing_fraction"])
    n_maj_absrb = I(ds.dim["absorber"])
    n_min_absrb = I(ds.dim["minor_absorber"])
    n_min_absrb_lower = I(ds.dim["minor_absorber_intervals_lower"])
    n_min_absrb_upper = I(ds.dim["minor_absorber_intervals_upper"])

    n_absrb_ext = I(ds.dim["absorber_ext"])
    n_contrib_lower = I(ds.dim["contributors_lower"])
    n_contrib_upper = I(ds.dim["contributors_upper"])

    p_ref_tropo = FT(ds["press_ref_trop"][:][1])
    t_ref_absrb = FT(ds["absorption_coefficient_ref_T"][:][1])
    p_ref_absrb = FT(ds["absorption_coefficient_ref_P"][:][1])

    gases_major = STA(undef, n_maj_absrb)
    gases_minor = STA(undef, n_min_absrb)
    id_minor = STA(undef, n_min_absrb)
    gases_minor_lower = STA(undef, n_min_absrb_lower)
    gases_minor_upper = STA(undef, n_min_absrb_upper)
    scaling_gas_lower = STA(undef, n_min_absrb_lower)
    scaling_gas_upper = STA(undef, n_min_absrb_upper)
    idx_gases = DSTAI()

    idx_gases_minor_lower = zeros(I, n_min_absrb_lower)
    idx_gases_minor_upper = zeros(I, n_min_absrb_upper)

    idx_scaling_gas_lower = zeros(I, n_min_absrb_lower)
    idx_scaling_gas_upper = zeros(I, n_min_absrb_upper)

    for igas = 1:n_maj_absrb
        gases_major[igas] = strip(String(ds["gas_names"][:, igas]))
        idx_gases[gases_major[igas]] = I(igas)
    end

    idx_h2o = idx_gases["h2o"]
    idx_gases["h2o_frgn"] = idx_h2o # water vapor - foreign
    idx_gases["h2o_self"] = idx_h2o # water vapor - self-continua


    n_gases = n_maj_absrb

    for igas = 1:n_min_absrb
        gases_minor[igas] = strip(String(ds["gas_minor"][:, igas]))
        id_minor[igas] = strip(String(ds["identifier_minor"][:, igas]))
    end

    for igas = 1:n_min_absrb_lower
        gases_minor_lower[igas] =
            strip(String(ds["minor_gases_lower"][:, igas]))
        scaling_gas_lower[igas] =
            strip(String(ds["scaling_gas_lower"][:, igas]))
        if ~isempty(gases_minor_lower[igas])
            idx_gases_minor_lower[igas] = idx_gases[gases_minor_lower[igas]]
        end
        if ~isempty(scaling_gas_lower[igas])
            idx_scaling_gas_lower[igas] = idx_gases[scaling_gas_lower[igas]]
        end
    end
    idx_gases_minor_lower = DA(idx_gases_minor_lower)
    idx_scaling_gas_lower = DA(idx_scaling_gas_lower)

    for igas = 1:n_min_absrb_upper
        gases_minor_upper[igas] =
            strip(String(ds["minor_gases_upper"][:, igas]))
        scaling_gas_upper[igas] =
            strip(String(ds["scaling_gas_upper"][:, igas]))
        if ~isempty(gases_minor_upper[igas])
            idx_gases_minor_upper[igas] = idx_gases[gases_minor_upper[igas]]
        end
        if ~isempty(scaling_gas_upper[igas])
            idx_scaling_gas_upper[igas] = idx_gases[scaling_gas_upper[igas]]
        end
    end
    idx_gases_minor_upper = DA(idx_gases_minor_upper)
    idx_scaling_gas_upper = DA(idx_scaling_gas_upper)

    key_species = IA3D(ds["key_species"][:])

    kmajor = FTA4D(ds["kmajor"][:])
    kminor_lower = FTA3D(ds["kminor_lower"][:])
    kminor_upper = FTA3D(ds["kminor_upper"][:])
    kminor_start_lower = IA1D(ds["kminor_start_lower"][:])
    kminor_start_upper = IA1D(ds["kminor_start_upper"][:])

    planck_fraction = FTA4D(ds["plank_fraction"][:])
    t_planck = FTA1D(ds["temperature_Planck"][:])

    totplnk = FTA2D(ds["totplnk"][:])
    bnd_lims_gpt = IA2D(ds["bnd_limits_gpt"][:])
    bnd_lims_wn = FTA2D(ds["bnd_limits_wavenumber"][:])
    #-----------------------
    major_gpt2bnd = Array{UI8,1}(undef, n_gpt)
    for i = 1:n_bnd
        major_gpt2bnd[bnd_lims_gpt[1, i]:bnd_lims_gpt[2, i]] .= UI8(i)
    end
    #-----------------------
    minor_lower_gpt_lims = IA2D(ds["minor_limits_gpt_lower"][:])
    minor_upper_gpt_lims = IA2D(ds["minor_limits_gpt_upper"][:])
    #-----------------------
    minor_lower_bnd = zeros(UI8, n_min_absrb_lower)
    minor_upper_bnd = zeros(UI8, n_min_absrb_upper)

    minor_lower_bnd_st = Array{UI8,1}(undef, n_bnd + 1)
    minor_upper_bnd_st = Array{UI8,1}(undef, n_bnd + 1)

    minor_lower_gpt_sh = Array{Int,1}(undef, n_min_absrb_lower)
    minor_upper_gpt_sh = Array{Int,1}(undef, n_min_absrb_upper)

    minor_lower_gpt_sh[1] = 0
    for i = 1:n_min_absrb_lower
        minor_lower_bnd[i] = major_gpt2bnd[minor_lower_gpt_lims[1, i]]
        if i > 1
            minor_lower_gpt_sh[i] =
                minor_lower_gpt_sh[i-1] + minor_lower_gpt_lims[2, i-1] -
                minor_lower_gpt_lims[1, i-1] + 1
        end
    end

    minor_upper_gpt_sh[1] = 0
    for i = 1:n_min_absrb_upper
        minor_upper_bnd[i] = major_gpt2bnd[minor_upper_gpt_lims[1, i]]
        if i > 1
            minor_upper_gpt_sh[i] =
                minor_upper_gpt_sh[i-1] + minor_upper_gpt_lims[2, i-1] -
                minor_upper_gpt_lims[1, i-1] + 1
        end
    end
    #-----------------------
    minor_lower_bnd_st[1] = 1
    minor_upper_bnd_st[1] = 1

    for ibnd = 2:n_bnd+1
        loc_low = findlast(isequal(UI8(ibnd - 1)), minor_lower_bnd)
        loc_upp = findlast(isequal(UI8(ibnd - 1)), minor_upper_bnd)
        if isnothing(loc_low)
            minor_lower_bnd_st[ibnd] = minor_lower_bnd_st[ibnd-1]
        else
            minor_lower_bnd_st[ibnd] = UI8(loc_low + 1)
        end

        if isnothing(loc_upp)
            minor_upper_bnd_st[ibnd] = minor_upper_bnd_st[ibnd-1]
        else
            minor_upper_bnd_st[ibnd] = UI8(loc_upp + 1)
        end
    end
    #-----------------------
    major_gpt2bnd = DA(major_gpt2bnd)
    minor_lower_bnd = DA(minor_lower_bnd)
    minor_upper_bnd = DA(minor_upper_bnd)
    minor_lower_bnd_st = DA(minor_lower_bnd_st)
    minor_upper_bnd_st = DA(minor_upper_bnd_st)
    minor_lower_gpt_sh = DA(minor_lower_gpt_sh)
    minor_upper_gpt_sh = DA(minor_upper_gpt_sh)
    #-----------------------

    minor_lower_scales_with_density =
        IA1D(ds["minor_scales_with_density_lower"][:])
    minor_upper_scales_with_density =
        IA1D(ds["minor_scales_with_density_upper"][:])

    lower_scale_by_complement = IA1D(ds["scale_by_complement_lower"][:])
    upper_scale_by_complement = IA1D(ds["scale_by_complement_upper"][:])

    p_ref = FTA1D(ds["press_ref"][:])
    t_ref = FTA1D(ds["temp_ref"][:])
    vmr_ref = FTA3D(ds["vmr_ref"][:])

    Δ_t_ref = t_ref[2] - t_ref[1]
    Δ_ln_p_ref = log(p_ref[1]) - log(p_ref[2])

    return (
        LookUpLW{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D}(
            n_gases,
            n_bnd,
            n_gpt,
            n_atmos_layers,
            n_t_ref,
            n_p_ref,
            n_t_plnk,
            n_mix_frac,
            n_maj_absrb,
            n_min_absrb,
            n_min_absrb_lower,
            n_min_absrb_upper,
            n_absrb_ext,
            n_contrib_lower,
            n_contrib_upper,
            idx_h2o,
            p_ref_tropo,
            t_ref_absrb,
            p_ref_absrb,
            Δ_t_ref,
            Δ_ln_p_ref,
            idx_gases_minor_lower,
            idx_gases_minor_upper,
            idx_scaling_gas_lower,
            idx_scaling_gas_upper,
            key_species,
            kmajor,
            kminor_lower,
            kminor_upper,
            kminor_start_lower,
            kminor_start_upper,
            planck_fraction,
            t_planck,
            totplnk,
            major_gpt2bnd,
            bnd_lims_gpt,
            bnd_lims_wn,
            minor_lower_gpt_lims,
            minor_upper_gpt_lims,
            minor_lower_bnd,
            minor_upper_bnd,
            minor_lower_bnd_st,
            minor_upper_bnd_st,
            minor_lower_gpt_sh,
            minor_upper_gpt_sh,
            minor_lower_scales_with_density,
            minor_lower_scales_with_density,
            lower_scale_by_complement,
            upper_scale_by_complement,
            p_ref,
            t_ref,
            vmr_ref,
        ),
        idx_gases,
    )
end

struct LookUpSW{
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
} <: AbstractLookUp{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D}
    n_gases::I
    n_bnd::I
    n_gpt::I
    n_atmos_layers::I
    n_t_ref::I
    n_p_ref::I
    n_mix_frac::I
    n_maj_absrb::I
    n_min_absrb::I
    n_min_absrb_lower::I
    n_min_absrb_upper::I
    n_absrb_ext::I
    n_contrib_lower::I
    n_contrib_upper::I
    idx_h2o::I

    p_ref_tropo::FT
    t_ref_absrb::FT
    p_ref_absrb::FT

    Δ_t_ref::FT
    Δ_ln_p_ref::FT

    idx_gases_minor_lower::IA1D
    idx_gases_minor_upper::IA1D

    idx_scaling_gas_lower::IA1D
    idx_scaling_gas_upper::IA1D

    key_species::IA3D

    kmajor::FTA4D
    kminor_lower::FTA3D
    kminor_upper::FTA3D

    kminor_start_lower::IA1D
    kminor_start_upper::IA1D

    major_gpt2bnd::UI8A1D
    bnd_lims_gpt::IA2D
    bnd_lims_wn::FTA2D

    minor_lower_gpt_lims::IA2D
    minor_upper_gpt_lims::IA2D

    minor_lower_bnd::UI8A1D # added---------------------------------
    minor_upper_bnd::UI8A1D

    minor_lower_bnd_st::UI8A1D
    minor_upper_bnd_st::UI8A1D

    minor_lower_gpt_sh::IA1D
    minor_upper_gpt_sh::IA1D #--------------------------------------

    minor_lower_scales_with_density::IA1D
    minor_upper_scales_with_density::IA1D

    lower_scale_by_complement::IA1D
    upper_scale_by_complement::IA1D

    p_ref::FTA1D
    t_ref::FTA1D
    vmr_ref::FTA3D

    rayl_lower::FTA3D
    rayl_upper::FTA3D
    solar_src::FTA1D
end
Adapt.@adapt_structure LookUpSW


function LookUpSW(
    ds,
    ::Type{I},
    ::Type{FT},
    ::Type{DA},
) where {I<:Int,FT<:AbstractFloat,DA}

    UI8 = UInt8
    UI8A1D = DA{UInt8,1}
    STA = Array{String,1}
    IA1D = DA{I,1}
    IA2D = DA{I,2}
    IA3D = DA{I,3}
    FTA1D = DA{FT,1}
    FTA2D = DA{FT,2}
    FTA3D = DA{FT,3}
    FTA4D = DA{FT,4}
    DSTAI = Dict{String,I}

    n_bnd = I(ds.dim["bnd"])
    n_gpt = I(ds.dim["gpt"])
    n_atmos_layers = I(ds.dim["atmos_layer"])
    n_t_ref = I(ds.dim["temperature"])
    n_p_ref = I(ds.dim["pressure"])
    n_mix_frac = I(ds.dim["mixing_fraction"])

    n_maj_absrb = I(ds.dim["absorber"])
    n_min_absrb = I(ds.dim["minor_absorber"])

    n_min_absrb_lower = I(ds.dim["minor_absorber_intervals_lower"])
    n_min_absrb_upper = I(ds.dim["minor_absorber_intervals_upper"])
    n_absrb_ext = I(ds.dim["absorber_ext"])

    n_contrib_lower = I(ds.dim["contributors_lower"])
    n_contrib_upper = I(ds.dim["contributors_upper"])

    p_ref_tropo = FT(ds["press_ref_trop"][:][1])
    t_ref_absrb = FT(ds["absorption_coefficient_ref_T"][:][1])
    p_ref_absrb = FT(ds["absorption_coefficient_ref_P"][:][1])

    gases_major = STA(undef, n_maj_absrb)
    gases_minor = STA(undef, n_min_absrb)
    id_minor = STA(undef, n_min_absrb)
    gases_minor_lower = STA(undef, n_min_absrb_lower)
    gases_minor_upper = STA(undef, n_min_absrb_upper)
    scaling_gas_lower = STA(undef, n_min_absrb_lower)
    scaling_gas_upper = STA(undef, n_min_absrb_upper)
    idx_gases = DSTAI()

    idx_gases_minor_lower = zeros(I, n_min_absrb_lower)
    idx_gases_minor_upper = zeros(I, n_min_absrb_upper)

    idx_scaling_gas_lower = zeros(I, n_min_absrb_lower)
    idx_scaling_gas_upper = zeros(I, n_min_absrb_upper)

    for igas = 1:n_maj_absrb
        gases_major[igas] = strip(String(ds["gas_names"][:, igas]))
        idx_gases[gases_major[igas]] = I(igas)
    end

    idx_h2o = idx_gases["h2o"]
    idx_gases["h2o_frgn"] = idx_h2o # water vapor - foreign
    idx_gases["h2o_self"] = idx_h2o # water vapor - self-continua

    n_gases = n_maj_absrb

    for igas = 1:n_min_absrb
        gases_minor[igas] = strip(String(ds["gas_minor"][:, igas]))
        id_minor[igas] = strip(String(ds["identifier_minor"][:, igas]))
    end

    for igas = 1:n_min_absrb_lower
        gases_minor_lower[igas] =
            strip(String(ds["minor_gases_lower"][:, igas]))
        scaling_gas_lower[igas] =
            strip(String(ds["scaling_gas_lower"][:, igas]))
        if ~isempty(gases_minor_lower[igas])
            idx_gases_minor_lower[igas] = idx_gases[gases_minor_lower[igas]]
        end
        if ~isempty(scaling_gas_lower[igas])
            idx_scaling_gas_lower[igas] = idx_gases[scaling_gas_lower[igas]]
        end
    end

    for igas = 1:n_min_absrb_upper
        gases_minor_upper[igas] =
            strip(String(ds["minor_gases_upper"][:, igas]))
        scaling_gas_upper[igas] =
            strip(String(ds["scaling_gas_upper"][:, igas]))
        if ~isempty(gases_minor_upper[igas])
            idx_gases_minor_upper[igas] = idx_gases[gases_minor_upper[igas]]
        end
        if ~isempty(scaling_gas_upper[igas])
            idx_scaling_gas_upper[igas] = idx_gases[scaling_gas_upper[igas]]
        end
    end

    key_species = IA3D(ds["key_species"][:])

    kmajor = FTA4D(ds["kmajor"][:])
    kminor_lower = FTA3D(ds["kminor_lower"][:])
    kminor_upper = FTA3D(ds["kminor_upper"][:])
    kminor_start_lower = IA1D(ds["kminor_start_lower"][:])
    kminor_start_upper = IA1D(ds["kminor_start_upper"][:])

    bnd_lims_gpt = IA2D(ds["bnd_limits_gpt"][:])
    bnd_lims_wn = FTA2D(ds["bnd_limits_wavenumber"][:])
    #-----------------------
    major_gpt2bnd = Array{UI8,1}(undef, n_gpt)
    for i = 1:n_bnd
        major_gpt2bnd[bnd_lims_gpt[1, i]:bnd_lims_gpt[2, i]] .= UI8(i)
    end
    #-----------------------
    minor_lower_gpt_lims = IA2D(ds["minor_limits_gpt_lower"][:])
    minor_upper_gpt_lims = IA2D(ds["minor_limits_gpt_upper"][:])
    #-----------------------
    minor_lower_bnd = zeros(UI8, n_min_absrb_lower)
    minor_upper_bnd = zeros(UI8, n_min_absrb_upper)

    minor_lower_bnd_st = Array{UI8,1}(undef, n_bnd + 1)
    minor_upper_bnd_st = Array{UI8,1}(undef, n_bnd + 1)

    minor_lower_gpt_sh = Array{Int,1}(undef, n_min_absrb_lower)
    minor_upper_gpt_sh = Array{Int,1}(undef, n_min_absrb_upper)

    minor_lower_gpt_sh[1] = 0
    for i = 1:n_min_absrb_lower
        minor_lower_bnd[i] = major_gpt2bnd[minor_lower_gpt_lims[1, i]]
        if i > 1
            minor_lower_gpt_sh[i] =
                minor_lower_gpt_sh[i-1] + minor_lower_gpt_lims[2, i-1] -
                minor_lower_gpt_lims[1, i-1] + 1
        end
    end

    minor_upper_gpt_sh[1] = 0
    for i = 1:n_min_absrb_upper
        minor_upper_bnd[i] = major_gpt2bnd[minor_upper_gpt_lims[1, i]]
        if i > 1
            minor_upper_gpt_sh[i] =
                minor_upper_gpt_sh[i-1] + minor_upper_gpt_lims[2, i-1] -
                minor_upper_gpt_lims[1, i-1] + 1
        end
    end
    #-----------------------
    minor_lower_bnd_st[1] = 1
    minor_upper_bnd_st[1] = 1

    for ibnd = 2:n_bnd+1
        loc_low = findlast(isequal(UI8(ibnd - 1)), minor_lower_bnd)
        loc_upp = findlast(isequal(UI8(ibnd - 1)), minor_upper_bnd)
        if isnothing(loc_low)
            minor_lower_bnd_st[ibnd] = minor_lower_bnd_st[ibnd-1]
        else
            minor_lower_bnd_st[ibnd] = UI8(loc_low + 1)
        end

        if isnothing(loc_upp)
            minor_upper_bnd_st[ibnd] = minor_upper_bnd_st[ibnd-1]
        else
            minor_upper_bnd_st[ibnd] = UI8(loc_upp + 1)
        end
    end
    #------------------------
    major_gpt2bnd = DA(major_gpt2bnd)
    minor_lower_bnd = DA(minor_lower_bnd)
    minor_upper_bnd = DA(minor_upper_bnd)
    minor_lower_bnd_st = DA(minor_lower_bnd_st)
    minor_upper_bnd_st = DA(minor_upper_bnd_st)
    minor_lower_gpt_sh = DA(minor_lower_gpt_sh)
    minor_upper_gpt_sh = DA(minor_upper_gpt_sh)
    #-----------------------    
    minor_lower_scales_with_density =
        IA1D(ds["minor_scales_with_density_lower"][:])
    minor_upper_scales_with_density =
        IA1D(ds["minor_scales_with_density_upper"][:])

    lower_scale_by_complement = IA1D(ds["scale_by_complement_lower"][:])
    upper_scale_by_complement = IA1D(ds["scale_by_complement_upper"][:])

    p_ref = FTA1D(ds["press_ref"][:])
    t_ref = FTA1D(ds["temp_ref"][:])
    vmr_ref = FTA3D(ds["vmr_ref"][:])

    Δ_t_ref = t_ref[2] - t_ref[1]
    Δ_ln_p_ref = log(p_ref[1]) - log(p_ref[2])

    rayl_lower = FTA3D(ds["rayl_lower"][:])
    rayl_upper = FTA3D(ds["rayl_upper"][:])
    solar_src = FTA1D(ds["solar_source"][:])

    return (
        LookUpSW{I,FT,UI8,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D}(
            n_gases,
            n_bnd,
            n_gpt,
            n_atmos_layers,
            n_t_ref,
            n_p_ref,
            n_mix_frac,
            n_maj_absrb,
            n_min_absrb,
            n_min_absrb_lower,
            n_min_absrb_upper,
            n_absrb_ext,
            n_contrib_lower,
            n_contrib_upper,
            idx_h2o,
            p_ref_tropo,
            t_ref_absrb,
            p_ref_absrb,
            Δ_t_ref,
            Δ_ln_p_ref,
            idx_gases_minor_lower,
            idx_gases_minor_upper,
            idx_scaling_gas_lower,
            idx_scaling_gas_upper,
            key_species,
            kmajor,
            kminor_lower,
            kminor_upper,
            kminor_start_lower,
            kminor_start_upper,
            major_gpt2bnd,
            bnd_lims_gpt,
            bnd_lims_wn,
            minor_lower_gpt_lims,
            minor_upper_gpt_lims,
            minor_lower_bnd,
            minor_upper_bnd,
            minor_lower_bnd_st,
            minor_upper_bnd_st,
            minor_lower_gpt_sh,
            minor_upper_gpt_sh,
            minor_lower_scales_with_density,
            minor_lower_scales_with_density,
            lower_scale_by_complement,
            upper_scale_by_complement,
            p_ref,
            t_ref,
            vmr_ref,
            rayl_lower,
            rayl_upper,
            solar_src,
        ),
        idx_gases,
    )
end

end
