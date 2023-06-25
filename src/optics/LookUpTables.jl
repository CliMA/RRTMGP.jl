module LookUpTables

using DocStringExtensions
using Adapt

export AbstractLookUp, LookUpLW, LookUpSW, LookUpCld

"""
    AbstractLookUp{FT}

Abstract lookup table for longwave and shortwave problems.
"""
abstract type AbstractLookUp{FT} end

"""
    LookUpLW{FT,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D} <: 
        AbstractLookUp{FT}

Longwave lookup tables, used to compute optical properties. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LookUpLW{
    FT <: AbstractFloat,
    UI8A1D <: AbstractArray{UInt8, 1},
    IA1D <: AbstractArray{Int, 1},
    IA2D <: AbstractArray{Int, 2},
    IA3D <: AbstractArray{Int, 3},
    FTA1D <: AbstractArray{FT, 1},
    FTA2D <: AbstractArray{FT, 2},
    FTA3D <: AbstractArray{FT, 3},
    FTA4D <: AbstractArray{FT, 4},
} <: AbstractLookUp{FT}
    "number of gases used in the lookup table"
    n_gases::Int
    "number of longwave bands"
    n_bnd::Int
    "number of `g-points`"
    n_gpt::Int
    "number of atmospheric layers (=2, lower and upper atmospheres)"
    n_atmos_layers::Int
    "number of reference temperatures for absorption coefficient lookup table"
    n_t_ref::Int
    "number of reference pressures for absorption lookup table"
    n_p_ref::Int
    "number of reference binary mixing fractions, for absorption coefficient lookup table"
    n_η::Int
    "number of reference temperatures, for Planck source calculations"
    n_t_plnk::Int
    "number of major absorbing gases"
    n_maj_absrb::Int
    "number of minor absorbing gases"
    n_min_absrb::Int
    "number of minor absorbers in lower atmosphere"
    n_min_absrb_lower::Int
    "number of minor absorbers in upper atmosphere"
    n_min_absrb_upper::Int
    n_absrb_ext::Int # not used
    "number of minor contributors in the lower atmosphere"
    n_contrib_lower::Int
    "number of minor contributors in the upper atmosphere"
    n_contrib_upper::Int
    "vmr array index for h2o"
    idx_h2o::Int
    "Reference pressure separating upper and lower atmosphere"
    p_ref_tropo::FT
    "Reference temperature"
    t_ref_absrb::FT
    "Reference pressure"
    p_ref_absrb::FT
    "minimum pressure supported by RRTMGP lookup tables"
    p_ref_min::FT
    "Δt for reference temperature values (Δt is constant)"
    Δ_t_ref::FT
    "Δ for log of reference pressure values (Δp is constant)"
    Δ_ln_p_ref::FT
    "indices for minor gases contributing to absorption in the lower atmosphere `(n_min_absrb_lower)`"
    idx_gases_minor_lower::IA1D
    "indices for minor gases contributing to absorption in the upper atmosphere `(n_min_absrb_upper)`"
    idx_gases_minor_upper::IA1D
    "indices for scaling gases in the lower atmosphere `(n_min_absrb_lower)`"
    idx_scaling_gas_lower::IA1D
    "indices for scaling gases in the upper atmosphere `(n_min_absrb_upper)`"
    idx_scaling_gas_upper::IA1D
    "major absorbing species in each band `(2, n_atmos_layers, n_bnd)`"
    key_species::IA3D
    "major absorption coefficient `(n_gpt, n_η, n_p_ref, n_t_ref)`"
    kmajor::FTA4D
    "minor absorption coefficient in lower atmosphere `(n_contrib_lower, n_η, n_t_ref)`"
    kminor_lower::FTA3D
    "minor absorption coefficient in upper atmosphere `(n_contrib_upper, n_η, n_t_ref)`"
    kminor_upper::FTA3D
    kminor_start_lower::IA1D # not currently used
    kminor_start_upper::IA1D # not currently used
    "Planck fraction `(n_gpt, n_η, n_p_ref, n_t_ref)`"
    planck_fraction::FTA4D
    "reference temperatures for Planck source calculations `(n_t_plnk)`"
    t_planck::FTA1D
    "total Planck source for each band `(n_t_plnk, n_bnd)`"
    totplnk::FTA2D
    "map from `g-point` to band"
    major_gpt2bnd::UI8A1D
    "starting and ending `g-point` for each band `(2, n_bnd)`"
    bnd_lims_gpt::IA2D
    "starting and ending wavenumber for each band `(2, n_bnd)`"
    bnd_lims_wn::FTA2D
    "`g-point` limits for minor contributors in lower atmosphere `(2, n_contrib_lower)`"
    minor_lower_gpt_lims::IA2D
    "`g-point` limits for minor contributors in upper atmosphere `(2, n_contrib_upper)`"
    minor_upper_gpt_lims::IA2D
    "band number for minor contributor in the lower atmosphere `(n_contrib_lower)`"
    minor_lower_bnd::UI8A1D
    "band number for minor contributor in the upper atmosphere `(n_contrib_upper)`"
    minor_upper_bnd::UI8A1D
    "starting index to `idx_gases_minor_lower` for each band `(n_bnd)`"
    minor_lower_bnd_st::UI8A1D
    "starting index to `idx_gases_minor_upper` for each band `(n_bnd)`"
    minor_upper_bnd_st::UI8A1D
    "shift in `kminor_lower` for each band `(n_min_absrb_lower)`"
    minor_lower_gpt_sh::IA1D
    "shift in `kminor_upper` for each band `(n_min_absrb_upper)`"
    minor_upper_gpt_sh::IA1D
    "minor gas (lower atmosphere) scales with density? `(n_min_absrb_lower)`"
    minor_lower_scales_with_density::IA1D
    "minor gas (upper atmosphere) scales with density? `(n_min_absrb_upper)`"
    minor_upper_scales_with_density::IA1D
    "minor gas (lower atmosphere) scales by compliment `(n_min_absrb_lower)`"
    lower_scale_by_complement::IA1D
    "minor gas (upper atmosphere) scales by compliment `(n_min_absrb_upper)`"
    upper_scale_by_complement::IA1D
    "reference pressures used by the lookup table `(n_p_ref)`"
    p_ref::FTA1D
    "reference temperatures used by the lookup table `(n_t_ref)`"
    t_ref::FTA1D
    "reference volume mixing ratios used by the lookup table `(2, n_gases, n_t_ref)`"
    vmr_ref::FTA3D
end
Adapt.@adapt_structure LookUpLW

function LookUpLW(ds, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}

    UI8 = UInt8
    UI8A1D = DA{UInt8, 1}
    STA = Array{String, 1}
    IA1D = DA{Int, 1}
    IA2D = DA{Int, 2}
    IA3D = DA{Int, 3}
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    FTA3D = DA{FT, 3}
    FTA4D = DA{FT, 4}
    DSTAI = Dict{String, Int}

    n_bnd = Int(ds.dim["bnd"])
    n_gpt = Int(ds.dim["gpt"])
    n_atmos_layers = Int(ds.dim["atmos_layer"])
    n_t_ref = Int(ds.dim["temperature"])
    n_p_ref = Int(ds.dim["pressure"])
    n_t_plnk = Int(ds.dim["temperature_Planck"])
    n_maj_absrb = Int(ds.dim["absorber"])
    n_min_absrb = Int(ds.dim["minor_absorber"])
    n_min_absrb_lower = Int(ds.dim["minor_absorber_intervals_lower"])
    n_min_absrb_upper = Int(ds.dim["minor_absorber_intervals_upper"])

    n_absrb_ext = Int(ds.dim["absorber_ext"])
    n_contrib_lower = Int(ds.dim["contributors_lower"])
    n_contrib_upper = Int(ds.dim["contributors_upper"])

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

    idx_gases_minor_lower = zeros(Int, n_min_absrb_lower)
    idx_gases_minor_upper = zeros(Int, n_min_absrb_upper)

    idx_scaling_gas_lower = zeros(Int, n_min_absrb_lower)
    idx_scaling_gas_upper = zeros(Int, n_min_absrb_upper)

    for igas in 1:n_maj_absrb
        gases_major[igas] = strip(String(ds["gas_names"][:, igas]))
        idx_gases[gases_major[igas]] = Int(igas)
    end

    idx_h2o = idx_gases["h2o"]
    idx_gases["h2o_frgn"] = idx_h2o # water vapor - foreign
    idx_gases["h2o_self"] = idx_h2o # water vapor - self-continua


    n_gases = n_maj_absrb

    for igas in 1:n_min_absrb
        gases_minor[igas] = strip(String(ds["gas_minor"][:, igas]))
        id_minor[igas] = strip(String(ds["identifier_minor"][:, igas]))
    end

    for igas in 1:n_min_absrb_lower
        gases_minor_lower[igas] = strip(String(ds["minor_gases_lower"][:, igas]))
        scaling_gas_lower[igas] = strip(String(ds["scaling_gas_lower"][:, igas]))
        if ~isempty(gases_minor_lower[igas])
            idx_gases_minor_lower[igas] = idx_gases[gases_minor_lower[igas]]
        end
        if ~isempty(scaling_gas_lower[igas])
            idx_scaling_gas_lower[igas] = idx_gases[scaling_gas_lower[igas]]
        end
    end
    idx_gases_minor_lower = DA(idx_gases_minor_lower)
    idx_scaling_gas_lower = DA(idx_scaling_gas_lower)

    for igas in 1:n_min_absrb_upper
        gases_minor_upper[igas] = strip(String(ds["minor_gases_upper"][:, igas]))
        scaling_gas_upper[igas] = strip(String(ds["scaling_gas_upper"][:, igas]))
        if ~isempty(gases_minor_upper[igas])
            idx_gases_minor_upper[igas] = idx_gases[gases_minor_upper[igas]]
        end
        if ~isempty(scaling_gas_upper[igas])
            idx_scaling_gas_upper[igas] = idx_gases[scaling_gas_upper[igas]]
        end
    end
    idx_gases_minor_upper = DA(idx_gases_minor_upper)
    idx_scaling_gas_upper = DA(idx_scaling_gas_upper)

    key_species = ds["key_species"][:]
    for j in 1:size(key_species, 3)
        for i in 1:size(key_species, 2)
            if key_species[1, i, j] == 0 && key_species[2, i, j] == 0
                key_species[1:2, i, j] .= 2
            end
        end
    end

    key_species = IA3D(key_species)

    kmajor = FTA4D(ds["kmajor"][:])
    kminor_lower = FTA3D(ds["kminor_lower"][:])
    kminor_upper = FTA3D(ds["kminor_upper"][:])
    kminor_start_lower = IA1D(ds["kminor_start_lower"][:])
    kminor_start_upper = IA1D(ds["kminor_start_upper"][:])

    planck_fraction = FTA4D(ds["plank_fraction"][:])
    t_planck = FTA1D(ds["temperature_Planck"][:])

    totplnk = FTA2D(ds["totplnk"][:])
    bnd_lims_gpt = Array{Int, 2}(ds["bnd_limits_gpt"][:])
    bnd_lims_wn = FTA2D(ds["bnd_limits_wavenumber"][:])
    #-----------------------
    major_gpt2bnd = Array{UI8, 1}(undef, n_gpt)
    for i in 1:n_bnd
        major_gpt2bnd[bnd_lims_gpt[1, i]:bnd_lims_gpt[2, i]] .= UI8(i)
    end
    #-----------------------
    bnd_lims_gpt = IA2D(bnd_lims_gpt)
    minor_lower_gpt_lims = Array{Int, 2}(ds["minor_limits_gpt_lower"][:])
    minor_upper_gpt_lims = Array{Int, 2}(ds["minor_limits_gpt_upper"][:])
    #-----------------------
    minor_lower_bnd = zeros(UI8, n_min_absrb_lower)
    minor_upper_bnd = zeros(UI8, n_min_absrb_upper)

    minor_lower_bnd_st = Array{UI8, 1}(undef, n_bnd + 1)
    minor_upper_bnd_st = Array{UI8, 1}(undef, n_bnd + 1)

    minor_lower_gpt_sh = Array{Int, 1}(undef, n_min_absrb_lower)
    minor_upper_gpt_sh = Array{Int, 1}(undef, n_min_absrb_upper)

    minor_lower_gpt_sh[1] = 0
    for i in 1:n_min_absrb_lower
        minor_lower_bnd[i] = major_gpt2bnd[minor_lower_gpt_lims[1, i]]
        if i > 1
            minor_lower_gpt_sh[i] =
                minor_lower_gpt_sh[i - 1] + minor_lower_gpt_lims[2, i - 1] - minor_lower_gpt_lims[1, i - 1] + 1
        end
    end

    minor_upper_gpt_sh[1] = 0
    for i in 1:n_min_absrb_upper
        minor_upper_bnd[i] = major_gpt2bnd[minor_upper_gpt_lims[1, i]]
        if i > 1
            minor_upper_gpt_sh[i] =
                minor_upper_gpt_sh[i - 1] + minor_upper_gpt_lims[2, i - 1] - minor_upper_gpt_lims[1, i - 1] + 1
        end
    end
    #-----------------------
    minor_lower_bnd_st[1] = 1
    minor_upper_bnd_st[1] = 1

    for ibnd in 2:(n_bnd + 1)
        loc_low = findlast(isequal(UI8(ibnd - 1)), minor_lower_bnd)
        loc_upp = findlast(isequal(UI8(ibnd - 1)), minor_upper_bnd)
        if isnothing(loc_low)
            minor_lower_bnd_st[ibnd] = minor_lower_bnd_st[ibnd - 1]
        else
            minor_lower_bnd_st[ibnd] = UI8(loc_low + 1)
        end

        if isnothing(loc_upp)
            minor_upper_bnd_st[ibnd] = minor_upper_bnd_st[ibnd - 1]
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
    minor_lower_gpt_lims = IA2D(minor_lower_gpt_lims)
    minor_upper_gpt_lims = IA2D(minor_upper_gpt_lims)
    #-----------------------

    minor_lower_scales_with_density = IA1D(ds["minor_scales_with_density_lower"][:])
    minor_upper_scales_with_density = IA1D(ds["minor_scales_with_density_upper"][:])

    lower_scale_by_complement = IA1D(ds["scale_by_complement_lower"][:])
    upper_scale_by_complement = IA1D(ds["scale_by_complement_upper"][:])

    p_ref = Array{FT, 1}(ds["press_ref"][:])
    t_ref = Array{FT, 1}(ds["temp_ref"][:])

    p_ref_min = minimum(p_ref)

    Δ_t_ref = t_ref[2] - t_ref[1]
    Δ_ln_p_ref = log(p_ref[1]) - log(p_ref[2])

    p_ref = FTA1D(p_ref)
    t_ref = FTA1D(t_ref)
    vmr_ref = FTA3D(ds["vmr_ref"][:])

    n_η = size(kmajor, 2)

    return (
        LookUpLW{FT, UI8A1D, IA1D, IA2D, IA3D, FTA1D, FTA2D, FTA3D, FTA4D}(
            n_gases,
            n_bnd,
            n_gpt,
            n_atmos_layers,
            n_t_ref,
            n_p_ref,
            n_η,
            n_t_plnk,
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
            p_ref_min,
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
            minor_upper_scales_with_density,
            lower_scale_by_complement,
            upper_scale_by_complement,
            p_ref,
            t_ref,
            vmr_ref,
        ),
        idx_gases,
    )
end

"""
    LookUpSW{FT,UI8A1D,IA1D,IA2D,IA3D,FTA1D,FTA2D,FTA3D,FTA4D} <:
        AbstractLookUp{FT}

Shortwave lookup tables, used to compute optical properties. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LookUpSW{
    FT <: AbstractFloat,
    UI8A1D <: AbstractArray{UInt8, 1},
    IA1D <: AbstractArray{Int, 1},
    IA2D <: AbstractArray{Int, 2},
    IA3D <: AbstractArray{Int, 3},
    FTA1D <: AbstractArray{FT, 1},
    FTA2D <: AbstractArray{FT, 2},
    FTA3D <: AbstractArray{FT, 3},
    FTA4D <: AbstractArray{FT, 4},
} <: AbstractLookUp{FT}
    "number of gases used in the lookup table"
    n_gases::Int
    "number of shortwave bands"
    n_bnd::Int
    "number of `g-points`"
    n_gpt::Int
    "number of atmospheric layers (=2, lower and upper atmospheres)"
    n_atmos_layers::Int
    "number of reference temperatures for absorption coefficient lookup table"
    n_t_ref::Int
    "number of reference pressures for absorption lookup table"
    n_p_ref::Int
    "number of reference binary mixing fractions, for absorption coefficient lookup table"
    n_η::Int
    "number of major absorbing gases"
    n_maj_absrb::Int
    "number of minor absorbing gases"
    n_min_absrb::Int
    "number of minor absorbers in lower atmosphere"
    n_min_absrb_lower::Int
    "number of minor absorbers in upper atmosphere"
    n_min_absrb_upper::Int
    n_absrb_ext::Int # not used
    "number of minor contributors in the lower atmosphere"
    n_contrib_lower::Int
    "number of minor contributors in the upper atmosphere"
    n_contrib_upper::Int
    "vmr array index for h2o"
    idx_h2o::Int
    "Reference pressure separating upper and lower atmosphere"
    p_ref_tropo::FT
    "Reference temperature"
    t_ref_absrb::FT
    "Reference pressure"
    p_ref_absrb::FT
    "minimum pressure supported by RRTMGP lookup tables"
    p_ref_min::FT
    "Δt for reference temperature values (Δt is constant, array is linearly equispaced)"
    Δ_t_ref::FT
    "Δ for log of reference pressure values (Δ is constant, equispaced on logarithmic scale)"
    Δ_ln_p_ref::FT
    "total solar irradiation"
    solar_src_tot::FT
    "indices for minor gases contributing to absorption in the lower atmosphere `(n_min_absrb_lower)`"
    idx_gases_minor_lower::IA1D
    "indices for minor gases contributing to absorption in the upper atmosphere `(n_min_absrb_upper)`"
    idx_gases_minor_upper::IA1D
    "indices for scaling gases in the lower atmosphere `(n_min_absrb_lower)`"
    idx_scaling_gas_lower::IA1D
    "indices for scaling gases in the upper atmosphere `(n_min_absrb_upper)`"
    idx_scaling_gas_upper::IA1D
    "major absorbing species in each band `(2, n_atmos_layers, n_bnd)`"
    key_species::IA3D
    "major absorption coefficient `(n_gpt, n_η, n_p_ref, n_t_ref)`"
    kmajor::FTA4D
    "minor absorption coefficient in lower atmosphere `(n_contrib_lower, n_η, n_t_ref)`"
    kminor_lower::FTA3D
    "minor absorption coefficient in upper atmosphere `(n_contrib_upper, n_η, n_t_ref)`"
    kminor_upper::FTA3D
    kminor_start_lower::IA1D # not currently used
    kminor_start_upper::IA1D # not currently used
    "map from g-point to band"
    major_gpt2bnd::UI8A1D
    "starting and ending `g-point` for each band `(2, n_bnd)`"
    bnd_lims_gpt::IA2D
    "starting and ending wavenumber for each band `(2, n_bnd)`"
    bnd_lims_wn::FTA2D
    "`g-point` limits for minor contributors in lower atmosphere `(2, n_contrib_lower)`"
    minor_lower_gpt_lims::IA2D
    "`g-point` limits for minor contributors in upper atmosphere `(2, n_contrib_upper)`"
    minor_upper_gpt_lims::IA2D
    "band number for minor contributor in the lower atmosphere `(n_contrib_lower)`"
    minor_lower_bnd::UI8A1D
    "band number for minor contributor in the upper atmosphere `(n_contrib_upper)`"
    minor_upper_bnd::UI8A1D
    "starting index to `idx_gases_minor_lower` for each band"
    minor_lower_bnd_st::UI8A1D
    "starting index to `idx_gases_minor_upper` for each band"
    minor_upper_bnd_st::UI8A1D
    "shift in `kminor_lower` for each band `(n_min_absrb_lower)`"
    minor_lower_gpt_sh::IA1D
    "shift in `kminor_upper` for each band `(n_min_absrb_upper)`"
    minor_upper_gpt_sh::IA1D
    "minor gas (lower atmosphere) scales with density? `(n_min_absrb_lower)`"
    minor_lower_scales_with_density::IA1D
    "minor gas (upper atmosphere) scales with density? `(n_min_absrb_upper)`"
    minor_upper_scales_with_density::IA1D
    "minor gas (lower atmosphere) scales by compliment `(n_min_absrb_lower)`"
    lower_scale_by_complement::IA1D
    "minor gas (upper atmosphere) scales by compliment `(n_min_absrb_upper)`"
    upper_scale_by_complement::IA1D
    "reference pressures used by the lookup table `(n_p_ref)`"
    p_ref::FTA1D
    "reference temperatures used by the lookup table `(n_t_ref)`"
    t_ref::FTA1D
    "reference volume mixing ratios used by the lookup table `(2, n_gases, n_t_ref)`"
    vmr_ref::FTA3D
    "Rayleigh absorption coefficient for lower atmosphere `(n_gpt, n_η, n_t_ref)`"
    rayl_lower::FTA3D
    "Rayleigh absorption coefficient for upper atmosphere `(n_gpt, n_η, n_t_ref)`"
    rayl_upper::FTA3D
    "relative solar source contribution from each `g-point` `(n_gpt)`"
    solar_src_scaled::FTA1D
end
Adapt.@adapt_structure LookUpSW


function LookUpSW(ds, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}

    UI8 = UInt8
    UI8A1D = DA{UInt8, 1}
    STA = Array{String, 1}
    IA1D = DA{Int, 1}
    IA2D = DA{Int, 2}
    IA3D = DA{Int, 3}
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    FTA3D = DA{FT, 3}
    FTA4D = DA{FT, 4}
    DSTAI = Dict{String, Int}

    n_bnd = Int(ds.dim["bnd"])
    n_gpt = Int(ds.dim["gpt"])
    n_atmos_layers = Int(ds.dim["atmos_layer"])
    n_t_ref = Int(ds.dim["temperature"])
    n_p_ref = Int(ds.dim["pressure"])

    n_maj_absrb = Int(ds.dim["absorber"])
    n_min_absrb = Int(ds.dim["minor_absorber"])

    n_min_absrb_lower = Int(ds.dim["minor_absorber_intervals_lower"])
    n_min_absrb_upper = Int(ds.dim["minor_absorber_intervals_upper"])
    n_absrb_ext = Int(ds.dim["absorber_ext"])

    n_contrib_lower = Int(ds.dim["contributors_lower"])
    n_contrib_upper = Int(ds.dim["contributors_upper"])

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

    idx_gases_minor_lower = zeros(Int, n_min_absrb_lower)
    idx_gases_minor_upper = zeros(Int, n_min_absrb_upper)

    idx_scaling_gas_lower = zeros(Int, n_min_absrb_lower)
    idx_scaling_gas_upper = zeros(Int, n_min_absrb_upper)

    for igas in 1:n_maj_absrb
        gases_major[igas] = strip(String(ds["gas_names"][:, igas]))
        idx_gases[gases_major[igas]] = Int(igas)
    end

    idx_h2o = idx_gases["h2o"]
    idx_gases["h2o_frgn"] = idx_h2o # water vapor - foreign
    idx_gases["h2o_self"] = idx_h2o # water vapor - self-continua

    n_gases = n_maj_absrb

    for igas in 1:n_min_absrb
        gases_minor[igas] = strip(String(ds["gas_minor"][:, igas]))
        id_minor[igas] = strip(String(ds["identifier_minor"][:, igas]))
    end

    for igas in 1:n_min_absrb_lower
        gases_minor_lower[igas] = strip(String(ds["minor_gases_lower"][:, igas]))
        scaling_gas_lower[igas] = strip(String(ds["scaling_gas_lower"][:, igas]))
        if ~isempty(gases_minor_lower[igas])
            idx_gases_minor_lower[igas] = idx_gases[gases_minor_lower[igas]]
        end
        if ~isempty(scaling_gas_lower[igas])
            idx_scaling_gas_lower[igas] = idx_gases[scaling_gas_lower[igas]]
        end
    end

    for igas in 1:n_min_absrb_upper
        gases_minor_upper[igas] = strip(String(ds["minor_gases_upper"][:, igas]))
        scaling_gas_upper[igas] = strip(String(ds["scaling_gas_upper"][:, igas]))
        if ~isempty(gases_minor_upper[igas])
            idx_gases_minor_upper[igas] = idx_gases[gases_minor_upper[igas]]
        end
        if ~isempty(scaling_gas_upper[igas])
            idx_scaling_gas_upper[igas] = idx_gases[scaling_gas_upper[igas]]
        end
    end

    key_species = ds["key_species"][:]
    for j in 1:size(key_species, 3)
        for i in 1:size(key_species, 2)
            if key_species[1, i, j] == 0 && key_species[2, i, j] == 0
                key_species[1:2, i, j] .= 2
            end
        end
    end

    key_species = IA3D(key_species)

    kmajor = FTA4D(ds["kmajor"][:])
    kminor_lower = FTA3D(ds["kminor_lower"][:])
    kminor_upper = FTA3D(ds["kminor_upper"][:])
    kminor_start_lower = IA1D(ds["kminor_start_lower"][:])
    kminor_start_upper = IA1D(ds["kminor_start_upper"][:])

    bnd_lims_gpt = Array{Int, 2}(ds["bnd_limits_gpt"][:])
    bnd_lims_wn = FTA2D(ds["bnd_limits_wavenumber"][:])
    #-----------------------
    major_gpt2bnd = Array{UI8, 1}(undef, n_gpt)
    for i in 1:n_bnd
        major_gpt2bnd[bnd_lims_gpt[1, i]:bnd_lims_gpt[2, i]] .= UI8(i)
    end
    #-----------------------
    bnd_lims_gpt = IA2D(bnd_lims_gpt)
    minor_lower_gpt_lims = Array{Int, 2}(ds["minor_limits_gpt_lower"][:])
    minor_upper_gpt_lims = Array{Int, 2}(ds["minor_limits_gpt_upper"][:])
    #-----------------------
    minor_lower_bnd = zeros(UI8, n_min_absrb_lower)
    minor_upper_bnd = zeros(UI8, n_min_absrb_upper)

    minor_lower_bnd_st = Array{UI8, 1}(undef, n_bnd + 1)
    minor_upper_bnd_st = Array{UI8, 1}(undef, n_bnd + 1)

    minor_lower_gpt_sh = Array{Int, 1}(undef, n_min_absrb_lower)
    minor_upper_gpt_sh = Array{Int, 1}(undef, n_min_absrb_upper)

    minor_lower_gpt_sh[1] = 0
    for i in 1:n_min_absrb_lower
        minor_lower_bnd[i] = major_gpt2bnd[minor_lower_gpt_lims[1, i]]
        if i > 1
            minor_lower_gpt_sh[i] =
                minor_lower_gpt_sh[i - 1] + minor_lower_gpt_lims[2, i - 1] - minor_lower_gpt_lims[1, i - 1] + 1
        end
    end

    minor_upper_gpt_sh[1] = 0
    for i in 1:n_min_absrb_upper
        minor_upper_bnd[i] = major_gpt2bnd[minor_upper_gpt_lims[1, i]]
        if i > 1
            minor_upper_gpt_sh[i] =
                minor_upper_gpt_sh[i - 1] + minor_upper_gpt_lims[2, i - 1] - minor_upper_gpt_lims[1, i - 1] + 1
        end
    end
    #-----------------------
    minor_lower_bnd_st[1] = 1
    minor_upper_bnd_st[1] = 1

    for ibnd in 2:(n_bnd + 1)
        loc_low = findlast(isequal(UI8(ibnd - 1)), minor_lower_bnd)
        loc_upp = findlast(isequal(UI8(ibnd - 1)), minor_upper_bnd)
        if isnothing(loc_low)
            minor_lower_bnd_st[ibnd] = minor_lower_bnd_st[ibnd - 1]
        else
            minor_lower_bnd_st[ibnd] = UI8(loc_low + 1)
        end

        if isnothing(loc_upp)
            minor_upper_bnd_st[ibnd] = minor_upper_bnd_st[ibnd - 1]
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
    minor_lower_gpt_lims = IA2D(minor_lower_gpt_lims)
    minor_upper_gpt_lims = IA2D(minor_upper_gpt_lims)
    #-----------------------    
    minor_lower_scales_with_density = IA1D(ds["minor_scales_with_density_lower"][:])
    minor_upper_scales_with_density = IA1D(ds["minor_scales_with_density_upper"][:])
    lower_scale_by_complement = IA1D(ds["scale_by_complement_lower"][:])
    upper_scale_by_complement = IA1D(ds["scale_by_complement_upper"][:])

    p_ref = Array{FT, 1}(ds["press_ref"][:])
    t_ref = Array{FT, 1}(ds["temp_ref"][:])

    p_ref_min = minimum(p_ref)

    Δ_t_ref = t_ref[2] - t_ref[1]
    Δ_ln_p_ref = log(p_ref[1]) - log(p_ref[2])

    p_ref = FTA1D(p_ref)
    t_ref = FTA1D(t_ref)
    vmr_ref = FTA3D(ds["vmr_ref"][:])

    n_η = size(kmajor, 2)

    rayl_lower = FTA3D(ds["rayl_lower"][:])
    rayl_upper = FTA3D(ds["rayl_upper"][:])

    solar_src = ds["solar_source"][:]
    solar_src_tot = FT(sum(solar_src))
    solar_src_scaled = FTA1D(solar_src ./ solar_src_tot)

    return (
        LookUpSW{FT, UI8A1D, IA1D, IA2D, IA3D, FTA1D, FTA2D, FTA3D, FTA4D}(
            n_gases,
            n_bnd,
            n_gpt,
            n_atmos_layers,
            n_t_ref,
            n_p_ref,
            n_η,
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
            p_ref_min,
            Δ_t_ref,
            Δ_ln_p_ref,
            solar_src_tot,
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
            minor_upper_scales_with_density,
            lower_scale_by_complement,
            upper_scale_by_complement,
            p_ref,
            t_ref,
            vmr_ref,
            rayl_lower,
            rayl_upper,
            solar_src_scaled,
        ),
        idx_gases,
    )
end

"""
    LookUpCld{FT,FTA1D,FTA2D,FTA3D,FTA4D}

Lookup table for cloud optics.

This struct stores the tables and Pade coefficients for determing extinction coeffient, 
single-scattering albedo, and asymmetry parameter g as a function of effective radius.
We compute the optical depth tau (=exintinction coeff * condensed water path)
and the products tau*ssa and tau*ssa*g for liquid and ice cloud separately.
These are used to determine the optical properties of ice and water cloud together.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LookUpCld{B, FT, FTA1D, FTA2D, FTA3D, FTA4D} <: AbstractLookUp{FT}
    "number of bands"
    nband::Int
    "number of ice roughness types"
    nrghice::Int
    "number of liquid particle sizes"
    nsize_liq::Int
    "number of ice particle sizes"
    nsize_ice::Int
    "number of size regimes"
    nsizereg::Int
    "number of extinction coefficients for pade approximation"
    ncoeff_ext::Int
    "number of ssa/g coefficients for pade approximation"
    ncoeff_ssa_g::Int
    "number of size regime boundaries for pade interpolation"
    nbound::Int
    "pair = 2"
    pair::Int
    "liquid particle size lower bound for LUT interpolation"
    radliq_lwr::FT
    "liquid particle size upper bound for LUT interpolation"
    radliq_upr::FT
    "factor for calculating LUT interpolation for liquid particle"
    radliq_fac::FT
    "ice particle size lower bound for LUT interpolation"
    radice_lwr::FT
    "ice particle size upper bound for LUT interpolation"
    radice_upr::FT
    "factor for calculating LUT interpolation for ice particle"
    radice_fac::FT
    "LUT liquid extinction coefficient (`nsize_liq, nbnd`) m²/g"
    lut_extliq::FTA2D
    "LUT liquid single scattering albedo (`nsize_liq, nbnd`)"
    lut_ssaliq::FTA2D
    "LUT liquid asymmetry parameter (`nsize_liq, nbnd`)"
    lut_asyliq::FTA2D
    "LUT ice extinction coefficient (`nsize_ice, nband, nrghice`) m²/g"
    lut_extice::FTA3D
    "LUT ice single scattering albedo (`nsize_ice, nband, nrghice`)"
    lut_ssaice::FTA3D
    "LUT ice asymmetry parameter (`nsize_ice, nband, nrghice`)"
    lut_asyice::FTA3D
    "pade coefficients for liquid extinction (`nband, nsizereg, ncoeff_ext`)"
    pade_extliq::FTA3D
    "pade coefficients for liquid single scattening albedo (`nband, nsizereg, ncoeff_ssa_g`)"
    pade_ssaliq::FTA3D
    "pade coefficients for liquid asymmetry paramter (`nband, nsizereg, ncoeff_ssa_g`)"
    pade_asyliq::FTA3D
    "pade coefficients for ice extinction (`nband, nsizereg, ncoeff_ext, nrghice`)"
    pade_extice::FTA4D
    "pade coefficients for ice single scattening albedo (`nband, nsizereg, ncoeff_ssa_g, nrghice`)"
    pade_ssaice::FTA4D
    "pade coefficients for ice asymmetry paramter (`nband, nsizereg, ncoeff_ssa_g, nrghice`)"
    pade_asyice::FTA4D
    "pade size regime boundaries for liquid extinction coefficient for pade interpolation (`nbound`) μm"
    pade_sizreg_extliq::FTA1D
    "pade size regime boundaries for liquid single scattering albedo for pade interpolation (`nbound`) μm"
    pade_sizreg_ssaliq::FTA1D
    "pade size regime boundaries for liquid asymmetry parameter for pade interpolation (`nbound`) μm"
    pade_sizreg_asyliq::FTA1D
    "pade size regime boundaries for ice extinction coefficient for pade interpolation (`nbound`) μm"
    pade_sizreg_extice::FTA1D
    "pade size regime boundaries for ice single scattering albedo for pade interpolation (`nbound`) μm"
    pade_sizreg_ssaice::FTA1D
    "pade size regime boundaries for ice asymmetry parameter for pade interpolation (`nbound`) μm"
    pade_sizreg_asyice::FTA1D
    "beginning and ending wavenumber for each band (`2, nband`) cm⁻¹"
    bnd_lims_wn::FTA2D
    "use LUT (default) or pade interpolation"
    use_lut::B
end
Adapt.@adapt_structure LookUpCld

function LookUpCld(ds, ::Type{FT}, ::Type{DA}, use_lut::Bool = true) where {FT <: AbstractFloat, DA}

    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    FTA3D = DA{FT, 3}
    FTA4D = DA{FT, 4}

    nband = Int(ds.dim["nband"])
    nrghice = Int(ds.dim["nrghice"])
    nsize_liq = Int(ds.dim["nsize_liq"])
    nsize_ice = Int(ds.dim["nsize_ice"])
    nsizereg = Int(ds.dim["nsizereg"])
    ncoeff_ext = Int(ds.dim["ncoeff_ext"])
    ncoeff_ssa_g = Int(ds.dim["ncoeff_ssa_g"])
    nbound = Int(ds.dim["nbound"])
    pair = Int(ds.dim["pair"])

    @assert nsizereg == 3 # RRTMGP pade approximation assumes exactly (3) size regimes

    radliq_lwr = FT(ds["radliq_lwr"][:])
    radliq_upr = FT(ds["radliq_upr"][:])
    radliq_fac = FT(ds["radliq_fac"][:])

    radice_lwr = FT(ds["radice_lwr"][:])
    radice_upr = FT(ds["radice_upr"][:])
    radice_fac = FT(ds["radice_fac"][:])

    lut_extliq = FTA2D(ds["lut_extliq"][:])
    lut_ssaliq = FTA2D(ds["lut_ssaliq"][:])
    lut_asyliq = FTA2D(ds["lut_asyliq"][:])

    lut_extice = FTA3D(ds["lut_extice"][:])
    lut_ssaice = FTA3D(ds["lut_ssaice"][:])
    lut_asyice = FTA3D(ds["lut_asyice"][:])

    pade_extliq = FTA3D(ds["pade_extliq"][:])
    pade_ssaliq = FTA3D(ds["pade_ssaliq"][:])
    pade_asyliq = FTA3D(ds["pade_asyliq"][:])

    pade_extice = FTA4D(ds["pade_extice"][:])
    pade_ssaice = FTA4D(ds["pade_ssaice"][:])
    pade_asyice = FTA4D(ds["pade_asyice"][:])

    pade_sizreg_extliq = FTA1D(ds["pade_sizreg_extliq"][:])
    pade_sizreg_ssaliq = FTA1D(ds["pade_sizreg_ssaliq"][:])
    pade_sizreg_asyliq = FTA1D(ds["pade_sizreg_asyliq"][:])

    pade_sizreg_extice = FTA1D(ds["pade_sizreg_extice"][:])
    pade_sizreg_ssaice = FTA1D(ds["pade_sizreg_ssaice"][:])
    pade_sizreg_asyice = FTA1D(ds["pade_sizreg_asyice"][:])

    bnd_lims_wn = FTA2D(ds["bnd_limits_wavenumber"][:])

    return (LookUpCld{Bool, FT, FTA1D, FTA2D, FTA3D, FTA4D}(
        nband,
        nrghice,
        nsize_liq,
        nsize_ice,
        nsizereg,
        ncoeff_ext,
        ncoeff_ssa_g,
        nbound,
        pair,
        radliq_lwr,
        radliq_upr,
        radliq_fac,
        radice_lwr,
        radice_upr,
        radice_fac,
        lut_extliq,
        lut_ssaliq,
        lut_asyliq,
        lut_extice,
        lut_ssaice,
        lut_asyice,
        pade_extliq,
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
        pade_sizreg_asyice,
        bnd_lims_wn,
        use_lut,
    ))

end

end
