module LookUpTables

using DocStringExtensions
using Adapt

export AbstractLookUp, LookUpLW, LookUpSW, LookUpCld, PadeCld, LookUpMinor, ReferencePoints, LookUpPlanck, BandData

"""
    AbstractLookUp

Abstract lookup table for longwave and shortwave problems.
"""
abstract type AbstractLookUp end

"""
    LookUpMinor{O, G, K}

Lookup table for computing optical properties of minor gases.
"""
struct LookUpMinor{O, G, K}
    "starting index to `idx_gases_minor_lower/upper` for each band `(n_bnd + 1)`"
    bnd_st::O
    "starting index in `kminor_lower/upper` for each g-point `(n_gpt + 1)`"
    gpt_st::O
    "contains indices for minor gases contributing to absorption in the lower/upper atmosphere,
     indices for scaling gases,
     minor gas scales with density,
     minor gas scales by complement `(4, n_min_absrb_lower/upper)`"
    gasdata::G
    "minor absorption coefficient in lower/upper atmosphere `(n_η, n_t_ref, n_contrib_lower/upper)`"
    kminor::K
end
Adapt.@adapt_structure LookUpMinor

# number of minor absorbors in the lower/upper atmosphere
@inline get_n_min_absrb(lkp::LookUpMinor) = size(lkp.idx_gases, 2)
# number of minor contributors in the lower/upper atmosphere
@inline get_n_contrib(lkp::LookUpMinor) = size(lkp.kminor, 3)
# get bounds for each g-point
@inline getbounds(lkp::LookUpMinor, ibnd, igpt) =
    @inbounds (lkp.bnd_st[ibnd], lkp.gpt_st[igpt], lkp.gpt_st[igpt + 1] - lkp.gpt_st[igpt])
# returns minor gas index, scaling gas index, minor gas scales with density and minor gas scales by complement at the given `idx`.
@inline get_minor_gas_data(lkp::LookUpMinor, idx) = @inbounds view(lkp.gasdata, :, idx)

"""
    ReferencePoints{RD1, RD3}

Log of reference pressures, reference temperatures, and volume mixing ratios used by the lookup table
"""
struct ReferencePoints{RD1, RD3}
    "log of reference pressures used by the lookup table `(n_p_ref)`"
    ln_p_ref::RD1
    "reference temperatures used by the lookup table `(n_t_ref)`"
    t_ref::RD1
    "reference volume mixing ratios used by the lookup table `(2, n_gases, n_t_ref)`"
    vmr_ref::RD3
end
Adapt.@adapt_structure ReferencePoints

"""
    LookUpPlanck{PD1, PD2, PD4}

Look up data for Planck source calculations.
"""
struct LookUpPlanck{PD1, PD2, PD4}
    "Planck fraction `(n_η, n_p_ref, n_t_ref, n_gpt)`"
    planck_fraction::PD4
    "reference temperatures for Planck source calculations `(n_t_plnk)`"
    t_planck::PD1
    "total Planck source for each band `(n_t_plnk, n_bnd)`"
    tot_planck::PD2
end
Adapt.@adapt_structure LookUpPlanck

"""
    BandData{BNDI1, BNDI2, BNDD2}

Band/g-point data for the lookup table.
"""
struct BandData{BNDI1, BNDI2, BNDD2}
    "map from `g-point` to band"
    major_gpt2bnd::BNDI1
    "starting and ending `g-point` for each band `(2, n_bnd)`"
    bnd_lims_gpt::BNDI2
    "starting and ending wavenumber for each band `(2, n_bnd)`"
    bnd_lims_wn::BNDD2
end
Adapt.@adapt_structure BandData

"""
    LookUpLW{FT, IA3D, FTA4D, BND, P, R, LMNR} <: AbstractLookUp

Longwave lookup tables, used to compute optical properties. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LookUpLW{FT, IA3D, FTA4D, BND, P, R, LMNR} <: AbstractLookUp
    "vmr array index for h2o"
    idx_h2o::Int
    "Reference pressure separating upper and lower atmosphere"
    p_ref_tropo::FT
    "minimum pressure supported by RRTMGP lookup tables"
    p_ref_min::FT
    "major absorbing species in each band `(2, n_atmos_layers, n_bnd)`"
    key_species::IA3D
    "major absorption coefficient `(n_η, n_p_ref + 1, n_t_ref, n_gpt)`"
    kmajor::FTA4D
    "lookup data for Planck source calculations"
    planck::P
    "band data"
    band_data::BND
    "reference temperatures, pressures and volume mixing ratios"
    ref_points::R
    "lookup data for minor gases in the lower atmosphere"
    minor_lower::LMNR
    "lookup data for minor gases in the upper atmosphere"
    minor_upper::LMNR
end
Adapt.@adapt_structure LookUpLW

function LookUpLW(ds, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}
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

    n_absrb_ext = Int(ds.dim["absorber_ext"]) # not currently used
    n_contrib_lower = Int(ds.dim["contributors_lower"])
    n_contrib_upper = Int(ds.dim["contributors_upper"])

    p_ref_tropo = FT(Array(ds["press_ref_trop"])[1])
    t_ref_absrb = FT(Array(ds["absorption_coefficient_ref_T"])[1]) # not currently used
    p_ref_absrb = FT(Array(ds["absorption_coefficient_ref_P"])[1]) # not currently used

    gases_major = STA(undef, n_maj_absrb)
    gases_minor = STA(undef, n_min_absrb)
    id_minor = STA(undef, n_min_absrb)
    gases_minor_lower = STA(undef, n_min_absrb_lower)
    gases_minor_upper = STA(undef, n_min_absrb_upper)
    scaling_gas_lower = STA(undef, n_min_absrb_lower)
    scaling_gas_upper = STA(undef, n_min_absrb_upper)
    idx_gases = DSTAI()

    idx_gases_minor_lower = zeros(Int, 1, n_min_absrb_lower)
    idx_gases_minor_upper = zeros(Int, 1, n_min_absrb_upper)

    idx_scaling_gas_lower = zeros(Int, 1, n_min_absrb_lower)
    idx_scaling_gas_upper = zeros(Int, 1, n_min_absrb_upper)

    @inbounds for igas in 1:n_maj_absrb
        gases_major[igas] = strip(String(ds["gas_names"][:, igas]))
        idx_gases[gases_major[igas]] = Int(igas)
    end

    idx_h2o = idx_gases["h2o"]
    idx_gases["h2o_frgn"] = idx_h2o # water vapor - foreign
    idx_gases["h2o_self"] = idx_h2o # water vapor - self-continua


    n_gases = n_maj_absrb

    @inbounds for igas in 1:n_min_absrb
        gases_minor[igas] = strip(String(ds["gas_minor"][:, igas]))
        id_minor[igas] = strip(String(ds["identifier_minor"][:, igas]))
    end

    @inbounds for igas in 1:n_min_absrb_lower
        gases_minor_lower[igas] = strip(String(ds["minor_gases_lower"][:, igas]))
        scaling_gas_lower[igas] = strip(String(ds["scaling_gas_lower"][:, igas]))
        if ~isempty(gases_minor_lower[igas])
            idx_gases_minor_lower[igas] = idx_gases[gases_minor_lower[igas]]
        end
        if ~isempty(scaling_gas_lower[igas])
            idx_scaling_gas_lower[igas] = idx_gases[scaling_gas_lower[igas]]
        end
    end
    idx_gases_minor_lower = idx_gases_minor_lower
    idx_scaling_gas_lower = idx_scaling_gas_lower

    @inbounds for igas in 1:n_min_absrb_upper
        gases_minor_upper[igas] = strip(String(ds["minor_gases_upper"][:, igas]))
        scaling_gas_upper[igas] = strip(String(ds["scaling_gas_upper"][:, igas]))
        if ~isempty(gases_minor_upper[igas])
            idx_gases_minor_upper[igas] = idx_gases[gases_minor_upper[igas]]
        end
        if ~isempty(scaling_gas_upper[igas])
            idx_scaling_gas_upper[igas] = idx_gases[scaling_gas_upper[igas]]
        end
    end
    idx_gases_minor_upper = idx_gases_minor_upper
    idx_scaling_gas_upper = idx_scaling_gas_upper

    key_species = Array(ds["key_species"])
    @inbounds for j in 1:size(key_species, 3)
        for i in 1:size(key_species, 2)
            if key_species[1, i, j] == 0 && key_species[2, i, j] == 0
                key_species[1:2, i, j] .= 2
            end
        end
    end

    key_species = IA3D(key_species)

    kmajor = FTA4D(permutedims(Array(ds["kmajor"]), [2, 3, 4, 1]))
    kminor_start_lower = IA1D(Array(ds["kminor_start_lower"])) # not currently used
    kminor_start_upper = IA1D(Array(ds["kminor_start_upper"])) # not currently used
    planck_fraction = FTA4D(permutedims(Array(ds["plank_fraction"]), [2, 3, 4, 1]))
    t_planck = FTA1D(Array(ds["temperature_Planck"]))

    tot_planck = FTA2D(Array(ds["totplnk"]))

    planck = LookUpPlanck(planck_fraction, t_planck, tot_planck)

    bnd_lims_gpt = Array{Int, 2}(Array(ds["bnd_limits_gpt"]))
    bnd_lims_wn = FTA2D(Array(ds["bnd_limits_wavenumber"]))
    #-----------------------
    major_gpt2bnd = Array{Int, 1}(undef, n_gpt)
    @inbounds for i in 1:n_bnd
        major_gpt2bnd[bnd_lims_gpt[1, i]:bnd_lims_gpt[2, i]] .= i
    end
    #-----------------------
    #`g-point` limits for minor contributors in lower/upper atmosphere `(2, n_min_absrb_lower/upper)`
    minor_lower_gpt_lims = Array{Int, 2}(Array(ds["minor_limits_gpt_lower"]))
    minor_upper_gpt_lims = Array{Int, 2}(Array(ds["minor_limits_gpt_upper"]))
    #-----------------------
    #band number for minor contributor in the lower/upper atmosphere `(n_min_absrb_lower)`
    minor_lower_bnd = zeros(Int, n_min_absrb_lower)
    minor_upper_bnd = zeros(Int, n_min_absrb_upper)

    minor_lower_bnd_st = Array{Int, 1}(undef, n_bnd + 1)
    minor_upper_bnd_st = Array{Int, 1}(undef, n_bnd + 1)

    minor_lower_gpt_sh = Array{Int, 1}(undef, n_min_absrb_lower)
    minor_upper_gpt_sh = Array{Int, 1}(undef, n_min_absrb_upper)

    minor_lower_gpt_st = ones(Int, n_gpt + 1)
    minor_upper_gpt_st = ones(Int, n_gpt + 1)

    minor_lower_gpt_sh[1] = 0
    @inbounds for i in 1:n_min_absrb_lower
        minor_lower_bnd[i] = major_gpt2bnd[minor_lower_gpt_lims[1, i]]
        if i > 1
            minor_lower_gpt_sh[i] =
                minor_lower_gpt_sh[i - 1] + minor_lower_gpt_lims[2, i - 1] - minor_lower_gpt_lims[1, i - 1] + 1
        end
    end

    minor_upper_gpt_sh[1] = 0
    @inbounds for i in 1:n_min_absrb_upper
        minor_upper_bnd[i] = major_gpt2bnd[minor_upper_gpt_lims[1, i]]
        if i > 1
            minor_upper_gpt_sh[i] =
                minor_upper_gpt_sh[i - 1] + minor_upper_gpt_lims[2, i - 1] - minor_upper_gpt_lims[1, i - 1] + 1
        end
    end
    #-----------------------
    minor_lower_bnd_st[1] = 1
    minor_upper_bnd_st[1] = 1

    @inbounds for ibnd in 2:(n_bnd + 1)
        loc_low = findlast(isequal(ibnd - 1), minor_lower_bnd)
        loc_upp = findlast(isequal(ibnd - 1), minor_upper_bnd)
        if isnothing(loc_low)
            minor_lower_bnd_st[ibnd] = minor_lower_bnd_st[ibnd - 1]
        else
            minor_lower_bnd_st[ibnd] = loc_low + 1
        end

        if isnothing(loc_upp)
            minor_upper_bnd_st[ibnd] = minor_upper_bnd_st[ibnd - 1]
        else
            minor_upper_bnd_st[ibnd] = loc_upp + 1
        end
    end
    # reorder kminor
    reorder_minor_lower = zeros(Int, n_contrib_lower)
    reorder_minor_upper = zeros(Int, n_contrib_upper)
    loc_lower, loc_upper = 1, 1
    for ibnd in 1:n_bnd
        nminorgases_lower = minor_lower_bnd_st[ibnd + 1] - minor_lower_bnd_st[ibnd]
        nminorgases_upper = minor_upper_bnd_st[ibnd + 1] - minor_upper_bnd_st[ibnd]
        for (loc_in_bnd, igpt) in enumerate(bnd_lims_gpt[1, ibnd]:bnd_lims_gpt[2, ibnd])
            minor_lower_gpt_st[igpt + 1] = minor_lower_gpt_st[igpt] + nminorgases_lower
            minor_upper_gpt_st[igpt + 1] = minor_upper_gpt_st[igpt] + nminorgases_upper
            for i in minor_lower_bnd_st[ibnd]:(minor_lower_bnd_st[ibnd + 1] - 1)
                reorder_minor_lower[loc_lower] = minor_lower_gpt_sh[i] + loc_in_bnd
                loc_lower += 1
            end
            for i in minor_upper_bnd_st[ibnd]:(minor_upper_bnd_st[ibnd + 1] - 1)
                reorder_minor_upper[loc_upper] = minor_upper_gpt_sh[i] + loc_in_bnd
                loc_upper += 1
            end
        end
    end

    kminor_lower = FTA3D(permutedims(Array(ds["kminor_lower"]), [2, 3, 1])[:, :, reorder_minor_lower])
    kminor_upper = FTA3D(permutedims(Array(ds["kminor_upper"]), [2, 3, 1])[:, :, reorder_minor_upper])
    #-----------------------
    major_gpt2bnd = DA(major_gpt2bnd)
    minor_lower_bnd = DA(minor_lower_bnd)
    minor_upper_bnd = DA(minor_upper_bnd)
    minor_lower_gpt_sh = DA(minor_lower_gpt_sh)
    minor_upper_gpt_sh = DA(minor_upper_gpt_sh)
    minor_lower_gpt_lims = IA2D(minor_lower_gpt_lims)
    minor_upper_gpt_lims = IA2D(minor_upper_gpt_lims)
    bnd_lims_gpt = IA2D(bnd_lims_gpt)
    #-----------------------
    band_data = BandData(major_gpt2bnd, bnd_lims_gpt, bnd_lims_wn)

    minor_lower_scales_with_density = reshape(Array(ds["minor_scales_with_density_lower"]), 1, :)
    minor_upper_scales_with_density = reshape(Array(ds["minor_scales_with_density_upper"]), 1, :)

    lower_scale_by_complement = reshape(Array(ds["scale_by_complement_lower"]), 1, :)
    upper_scale_by_complement = reshape(Array(ds["scale_by_complement_upper"]), 1, :)

    p_ref = Array{FT, 1}(Array(ds["press_ref"]))
    t_ref = Array{FT, 1}(Array(ds["temp_ref"]))

    p_ref_min = minimum(p_ref)

    Δ_t_ref = t_ref[2] - t_ref[1]
    Δ_ln_p_ref = log(p_ref[1]) - log(p_ref[2])

    ln_p_ref = FTA1D(log.(p_ref))
    t_ref = FTA1D(t_ref)
    vmr_ref = FTA3D(Array(ds["vmr_ref"]))

    n_η = size(kmajor, 1)
    minor_lower = LookUpMinor(
        DA(minor_lower_bnd_st),
        DA(minor_lower_gpt_st),
        DA(
            vcat(
                idx_gases_minor_lower,
                idx_scaling_gas_lower,
                minor_lower_scales_with_density,
                lower_scale_by_complement,
            ),
        ),
        kminor_lower,
    )

    minor_upper = LookUpMinor(
        DA(minor_upper_bnd_st),
        DA(minor_upper_gpt_st),
        DA(
            vcat(
                idx_gases_minor_upper,
                idx_scaling_gas_upper,
                minor_upper_scales_with_density,
                upper_scale_by_complement,
            ),
        ),
        kminor_upper,
    )
    ref_points = ReferencePoints(ln_p_ref, t_ref, vmr_ref)
    return (
        LookUpLW{FT, IA3D, FTA4D, typeof(band_data), typeof(planck), typeof(ref_points), typeof(minor_lower)}(
            idx_h2o,
            p_ref_tropo,
            p_ref_min,
            key_species,
            kmajor,
            planck,
            band_data,
            ref_points,
            minor_lower,
            minor_upper,
        ),
        idx_gases,
    )
end

# number of minor absorbors in the lower/upper atmosphere
@inline get_n_min_absrb_lower(lkp::LookUpLW) = get_n_min_absrb(lkp.minor_lower)
@inline get_n_min_absrb_upper(lkp::LookUpLW) = get_n_min_absrb(lkp.minor_upper)
# number of minor contributors in the lower/upper atmosphere
@inline get_n_contrib_lower(lkp::LookUpLW) = get_n_contrib(lkp.minor_lower)
@inline get_n_contrib_upper(lkp::LookUpLW) = get_n_contrib(lkp.minor_upper)
# number of bands in the longwave lookup table
@inline get_n_bnd(lkp::LookUpLW) = size(lkp.key_species, 3)
# number of atmospheric layers (=2, lower and upper atmospheres)
@inline get_n_atmos_layers(lkp::LookUpLW) = size(lkp.key_species, 2)
# number of gases used in the lookup table
@inline get_n_gases(lkp::LookUpLW) = size(lkp.ref_points.vmr_ref, 2)

@inline get_n_η(lkp::LookUpLW) = size(lkp.kmajor, 1)

"""
    LookUpSW{FT, IA3D, FTA1D, FTA3D, FTA4D, BND, R, LMNR} <: AbstractLookUp

Shortwave lookup tables, used to compute optical properties. 

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LookUpSW{FT, IA3D, FTA1D, FTA3D, FTA4D, BND, R, LMNR} <: AbstractLookUp
    "vmr array index for h2o"
    idx_h2o::Int
    "Reference pressure separating upper and lower atmosphere"
    p_ref_tropo::FT
    "minimum pressure supported by RRTMGP lookup tables"
    p_ref_min::FT
    "total solar irradiation"
    solar_src_tot::FT
    "major absorbing species in each band `(2, n_atmos_layers, n_bnd)`"
    key_species::IA3D
    "major absorption coefficient `(n_η, n_p_ref + 1, n_t_ref, n_gpt)`"
    kmajor::FTA4D
    "band data"
    band_data::BND
    "reference temperatures, pressures and volume mixing ratios"
    ref_points::R
    "Rayleigh absorption coefficient for lower atmosphere `(n_η, n_t_ref, n_gpt)`"
    rayl_lower::FTA3D
    "Rayleigh absorption coefficient for upper atmosphere `(n_η, n_t_ref, n_gpt)`"
    rayl_upper::FTA3D
    "relative solar source contribution from each `g-point` `(n_gpt)`"
    solar_src_scaled::FTA1D
    "lookup data for minor gases in the lower atmosphere"
    minor_lower::LMNR
    "lookup data for minor gases in the upper atmosphere"
    minor_upper::LMNR
end
Adapt.@adapt_structure LookUpSW


function LookUpSW(ds, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}
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
    n_absrb_ext = Int(ds.dim["absorber_ext"]) # not currently used

    n_contrib_lower = Int(ds.dim["contributors_lower"])
    n_contrib_upper = Int(ds.dim["contributors_upper"])

    p_ref_tropo = FT(Array(ds["press_ref_trop"])[1])
    t_ref_absrb = FT(Array(ds["absorption_coefficient_ref_T"])[1]) # not currently used
    p_ref_absrb = FT(Array(ds["absorption_coefficient_ref_P"])[1]) # not currently used

    gases_major = STA(undef, n_maj_absrb)
    gases_minor = STA(undef, n_min_absrb)
    id_minor = STA(undef, n_min_absrb)
    gases_minor_lower = STA(undef, n_min_absrb_lower)
    gases_minor_upper = STA(undef, n_min_absrb_upper)
    scaling_gas_lower = STA(undef, n_min_absrb_lower)
    scaling_gas_upper = STA(undef, n_min_absrb_upper)
    idx_gases = DSTAI()

    idx_gases_minor_lower = zeros(Int, 1, n_min_absrb_lower)
    idx_gases_minor_upper = zeros(Int, 1, n_min_absrb_upper)

    idx_scaling_gas_lower = zeros(Int, 1, n_min_absrb_lower)
    idx_scaling_gas_upper = zeros(Int, 1, n_min_absrb_upper)

    @inbounds for igas in 1:n_maj_absrb
        gases_major[igas] = strip(String(ds["gas_names"][:, igas]))
        idx_gases[gases_major[igas]] = Int(igas)
    end

    idx_h2o = idx_gases["h2o"]
    idx_gases["h2o_frgn"] = idx_h2o # water vapor - foreign
    idx_gases["h2o_self"] = idx_h2o # water vapor - self-continua

    n_gases = n_maj_absrb

    @inbounds for igas in 1:n_min_absrb
        gases_minor[igas] = strip(String(ds["gas_minor"][:, igas]))
        id_minor[igas] = strip(String(ds["identifier_minor"][:, igas]))
    end

    @inbounds for igas in 1:n_min_absrb_lower
        gases_minor_lower[igas] = strip(String(ds["minor_gases_lower"][:, igas]))
        scaling_gas_lower[igas] = strip(String(ds["scaling_gas_lower"][:, igas]))
        if ~isempty(gases_minor_lower[igas])
            idx_gases_minor_lower[igas] = idx_gases[gases_minor_lower[igas]]
        end
        if ~isempty(scaling_gas_lower[igas])
            idx_scaling_gas_lower[igas] = idx_gases[scaling_gas_lower[igas]]
        end
    end

    @inbounds for igas in 1:n_min_absrb_upper
        gases_minor_upper[igas] = strip(String(ds["minor_gases_upper"][:, igas]))
        scaling_gas_upper[igas] = strip(String(ds["scaling_gas_upper"][:, igas]))
        if ~isempty(gases_minor_upper[igas])
            idx_gases_minor_upper[igas] = idx_gases[gases_minor_upper[igas]]
        end
        if ~isempty(scaling_gas_upper[igas])
            idx_scaling_gas_upper[igas] = idx_gases[scaling_gas_upper[igas]]
        end
    end

    key_species = Array(ds["key_species"])
    @inbounds for j in 1:size(key_species, 3)
        for i in 1:size(key_species, 2)
            if key_species[1, i, j] == 0 && key_species[2, i, j] == 0
                key_species[1:2, i, j] .= 2
            end
        end
    end

    key_species = IA3D(key_species)

    kmajor = FTA4D(permutedims(Array(ds["kmajor"]), [2, 3, 4, 1]))
    kminor_start_lower = IA1D(Array(ds["kminor_start_lower"])) # not currently used
    kminor_start_upper = IA1D(Array(ds["kminor_start_upper"])) # not currently used

    bnd_lims_gpt = Array{Int, 2}(Array(ds["bnd_limits_gpt"]))
    bnd_lims_wn = FTA2D(Array(ds["bnd_limits_wavenumber"]))
    #-----------------------
    major_gpt2bnd = Array{Int, 1}(undef, n_gpt)
    @inbounds for i in 1:n_bnd
        major_gpt2bnd[bnd_lims_gpt[1, i]:bnd_lims_gpt[2, i]] .= i
    end
    #-----------------------
    minor_lower_gpt_lims = Array{Int, 2}(Array(ds["minor_limits_gpt_lower"]))
    minor_upper_gpt_lims = Array{Int, 2}(Array(ds["minor_limits_gpt_upper"]))
    #-----------------------
    minor_lower_bnd = zeros(Int, n_min_absrb_lower)
    minor_upper_bnd = zeros(Int, n_min_absrb_upper)

    minor_lower_bnd_st = Array{Int, 1}(undef, n_bnd + 1)
    minor_upper_bnd_st = Array{Int, 1}(undef, n_bnd + 1)

    minor_lower_gpt_sh = Array{Int, 1}(undef, n_min_absrb_lower)
    minor_upper_gpt_sh = Array{Int, 1}(undef, n_min_absrb_upper)

    minor_lower_gpt_st = ones(Int, n_gpt + 1)
    minor_upper_gpt_st = ones(Int, n_gpt + 1)

    minor_lower_gpt_sh[1] = 0
    @inbounds for i in 1:n_min_absrb_lower
        minor_lower_bnd[i] = major_gpt2bnd[minor_lower_gpt_lims[1, i]]
        if i > 1
            minor_lower_gpt_sh[i] =
                minor_lower_gpt_sh[i - 1] + minor_lower_gpt_lims[2, i - 1] - minor_lower_gpt_lims[1, i - 1] + 1
        end
    end

    minor_upper_gpt_sh[1] = 0
    @inbounds for i in 1:n_min_absrb_upper
        minor_upper_bnd[i] = major_gpt2bnd[minor_upper_gpt_lims[1, i]]
        if i > 1
            minor_upper_gpt_sh[i] =
                minor_upper_gpt_sh[i - 1] + minor_upper_gpt_lims[2, i - 1] - minor_upper_gpt_lims[1, i - 1] + 1
        end
    end
    #-----------------------
    minor_lower_bnd_st[1] = 1
    minor_upper_bnd_st[1] = 1

    @inbounds for ibnd in 2:(n_bnd + 1)
        loc_low = findlast(isequal(ibnd - 1), minor_lower_bnd)
        loc_upp = findlast(isequal(ibnd - 1), minor_upper_bnd)
        if isnothing(loc_low)
            minor_lower_bnd_st[ibnd] = minor_lower_bnd_st[ibnd - 1]
        else
            minor_lower_bnd_st[ibnd] = loc_low + 1
        end

        if isnothing(loc_upp)
            minor_upper_bnd_st[ibnd] = minor_upper_bnd_st[ibnd - 1]
        else
            minor_upper_bnd_st[ibnd] = loc_upp + 1
        end
    end
    # reorder kminor
    reorder_minor_lower = zeros(Int, n_contrib_lower)
    reorder_minor_upper = zeros(Int, n_contrib_upper)
    loc_lower, loc_upper = 1, 1
    for ibnd in 1:n_bnd
        nminorgases_lower = minor_lower_bnd_st[ibnd + 1] - minor_lower_bnd_st[ibnd]
        nminorgases_upper = minor_upper_bnd_st[ibnd + 1] - minor_upper_bnd_st[ibnd]
        for (loc_in_bnd, igpt) in enumerate(bnd_lims_gpt[1, ibnd]:bnd_lims_gpt[2, ibnd])
            minor_lower_gpt_st[igpt + 1] = minor_lower_gpt_st[igpt] + nminorgases_lower
            minor_upper_gpt_st[igpt + 1] = minor_upper_gpt_st[igpt] + nminorgases_upper
            for i in minor_lower_bnd_st[ibnd]:(minor_lower_bnd_st[ibnd + 1] - 1)
                reorder_minor_lower[loc_lower] = minor_lower_gpt_sh[i] + loc_in_bnd
                loc_lower += 1
            end
            for i in minor_upper_bnd_st[ibnd]:(minor_upper_bnd_st[ibnd + 1] - 1)
                reorder_minor_upper[loc_upper] = minor_upper_gpt_sh[i] + loc_in_bnd
                loc_upper += 1
            end
        end
    end

    kminor_lower = FTA3D(permutedims(Array(ds["kminor_lower"]), [2, 3, 1])[:, :, reorder_minor_lower])
    kminor_upper = FTA3D(permutedims(Array(ds["kminor_upper"]), [2, 3, 1])[:, :, reorder_minor_upper])
    #------------------------
    major_gpt2bnd = DA(major_gpt2bnd)
    minor_lower_bnd = DA(minor_lower_bnd)
    minor_upper_bnd = DA(minor_upper_bnd)
    minor_lower_gpt_sh = DA(minor_lower_gpt_sh)
    minor_upper_gpt_sh = DA(minor_upper_gpt_sh)
    minor_lower_gpt_lims = IA2D(minor_lower_gpt_lims)
    minor_upper_gpt_lims = IA2D(minor_upper_gpt_lims)
    bnd_lims_gpt = IA2D(bnd_lims_gpt)
    #-----------------------    
    band_data = BandData(major_gpt2bnd, bnd_lims_gpt, bnd_lims_wn)

    minor_lower_scales_with_density = reshape(Array(ds["minor_scales_with_density_lower"]), 1, :)
    minor_upper_scales_with_density = reshape(Array(ds["minor_scales_with_density_upper"]), 1, :)
    lower_scale_by_complement = reshape(Array(ds["scale_by_complement_lower"]), 1, :)
    upper_scale_by_complement = reshape(Array(ds["scale_by_complement_upper"]), 1, :)

    p_ref = Array{FT, 1}(Array(ds["press_ref"]))
    t_ref = Array{FT, 1}(Array(ds["temp_ref"]))

    p_ref_min = minimum(p_ref)

    Δ_t_ref = t_ref[2] - t_ref[1]
    Δ_ln_p_ref = log(p_ref[1]) - log(p_ref[2])

    ln_p_ref = FTA1D(log.(p_ref))
    t_ref = FTA1D(t_ref)
    vmr_ref = FTA3D(Array(ds["vmr_ref"]))

    n_η = size(kmajor, 1)

    rayl_lower = FTA3D(permutedims(Array(ds["rayl_lower"]), [2, 3, 1]))
    rayl_upper = FTA3D(permutedims(Array(ds["rayl_upper"]), [2, 3, 1]))

    solar_src = Array(ds["solar_source"])
    solar_src_tot = FT(sum(solar_src))
    solar_src_scaled = FTA1D(solar_src ./ solar_src_tot)

    minor_lower = LookUpMinor(
        DA(minor_lower_bnd_st),
        DA(minor_lower_gpt_st),
        DA(
            vcat(
                idx_gases_minor_lower,
                idx_scaling_gas_lower,
                minor_lower_scales_with_density,
                lower_scale_by_complement,
            ),
        ),
        kminor_lower,
    )

    minor_upper = LookUpMinor(
        DA(minor_upper_bnd_st),
        DA(minor_upper_gpt_st),
        DA(
            vcat(
                idx_gases_minor_upper,
                idx_scaling_gas_upper,
                minor_upper_scales_with_density,
                upper_scale_by_complement,
            ),
        ),
        kminor_upper,
    )
    ref_points = ReferencePoints(ln_p_ref, t_ref, vmr_ref)

    return (
        LookUpSW{FT, IA3D, FTA1D, FTA3D, FTA4D, typeof(band_data), typeof(ref_points), typeof(minor_lower)}(
            idx_h2o,
            p_ref_tropo,
            p_ref_min,
            solar_src_tot,
            key_species,
            kmajor,
            band_data,
            ref_points,
            rayl_lower,
            rayl_upper,
            solar_src_scaled,
            minor_lower,
            minor_upper,
        ),
        idx_gases,
    )
end

# number of minor absorbors in the lower/upper atmosphere
@inline get_n_min_absrb_lower(lkp::LookUpSW) = get_n_min_absrb(lkp.minor_lower)
@inline get_n_min_absrb_upper(lkp::LookUpSW) = get_n_min_absrb(lkp.minor_upper)
# number of minor contributors in the lower/upper atmosphere
@inline get_n_contrib_lower(lkp::LookUpSW) = get_n_contrib(lkp.minor_lower)
@inline get_n_contrib_upper(lkp::LookUpSW) = get_n_contrib(lkp.minor_upper)
# number of bands in the shortwave lookup table
@inline get_n_bnd(lkp::LookUpSW) = size(lkp.key_species, 3)
# number of atmospheric layers (=2, lower and upper atmospheres)
@inline get_n_atmos_layers(lkp::LookUpSW) = size(lkp.key_species, 2)
# number of gases used in the lookup table
@inline get_n_gases(lkp::LookUpSW) = size(lkp.ref_points.vmr_ref, 2)

@inline get_n_η(lkp::LookUpSW) = size(lkp.kmajor, 1)
"""
    LookUpCld{D, B, L, I, W} <: AbstractLookUp

Lookup table for cloud optics.

This struct stores the lookup tables for determing extinction coeffient, 
single-scattering albedo, and asymmetry parameter g as a function of effective radius.
We compute the optical depth tau (=exintinction coeff * condensed water path)
and the products tau*ssa and tau*ssa*g for liquid and ice cloud separately.
These are used to determine the optical properties of ice and water cloud together.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LookUpCld{D, B, L, I, W} <: AbstractLookUp
    "`nband`, `nrghice`, `nsize_liq`, `nsize_ice`, `pair`"
    dims::D
    "particle size lower and upper bounds and factor for LUT interpolation for liquid and ice particles"
    bounds::B
    "liquid extinction coefficient, single scattering albedo and symmetry paramter `(3*nsize_liq, nbnd)`"
    liqdata::L
    "ice extinction coefficient, single scattering albedo and symmetry paramter `(3*nsize_ice, nbnd, nrghice)`"
    icedata::I
    "beginning and ending wavenumber for each band (`2, nband`) cm⁻¹"
    bnd_lims_wn::W
end
Adapt.@adapt_structure LookUpCld

function LookUpCld(ds, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}
    dims = DA([
        Int(ds.dim["nband"]),
        Int(ds.dim["nrghice"]),
        Int(ds.dim["nsize_liq"]),
        Int(ds.dim["nsize_ice"]),
        Int(ds.dim["pair"]),
    ])
    bounds = DA([
        # liquid particle size lower bound for LUT interpolation
        FT(Array(ds["radliq_lwr"])),
        # liquid particle size upper bound for LUT interpolation
        FT(Array(ds["radliq_upr"])),
        # factor for calculating LUT interpolation for liquid particle
        FT(Array(ds["radliq_fac"])),
        # ice particle size lower bound for LUT interpolation
        FT(Array(ds["radice_lwr"])),
        # ice particle size upper bound for LUT interpolation
        FT(Array(ds["radice_upr"])),
        # factor for calculating LUT interpolation for ice particle
        FT(Array(ds["radice_fac"])),
    ])
    liqdata = DA(vcat(Array(ds["lut_extliq"]), Array(ds["lut_ssaliq"]), Array(ds["lut_asyliq"])))
    icedata = DA(vcat(Array(ds["lut_extice"]), Array(ds["lut_ssaice"]), Array(ds["lut_asyice"])))
    bnd_lims_wn = DA(Array(ds["bnd_limits_wavenumber"]))
    return LookUpCld(dims, bounds, liqdata, icedata, bnd_lims_wn)
end

@inline get_dims(lkp::LookUpCld) = lkp.dims
# number of bands (cloud lookup table dimension)
@inline get_nband(lkp::LookUpCld) = @inbounds lkp.dims[1]
# number of ice roughness types (cloud lookup table dimension)
@inline get_nrghice(lkp::LookUpCld) = @inbounds lkp.dims[2]
# number of liquid particle sizes (cloud lookup table dimension)
@inline get_nsize_liq(lkp::LookUpCld) = @inbounds lkp.dims[3]
# number of ice particle sizes (cloud lookup table dimension)
@inline get_nsize_ice(lkp::LookUpCld) = @inbounds lkp.dims[4]
# pair = 2 (cloud lookup table dimension)
@inline get_pair(::LookUpCld) = @inbounds lkp.dims[5]

@inline function getview_liqdata(lkp::LookUpCld, ibnd)
    n = get_nsize_liq(lkp)
    liqdata = lkp.liqdata
    return @inbounds (
        # LUT liquid extinction coefficient (`nsize_liq, nbnd`) m²/g
        view(liqdata, 1:n, ibnd),
        # LUT liquid single scattering albedo (`nsize_liq, nbnd`)
        view(liqdata, (n + 1):(n * 2), ibnd),
        # LUT liquid asymmetry parameter (`nsize_liq, nbnd`)
        view(liqdata, (n * 2 + 1):(n * 3), ibnd),
    )
end

@inline function getview_icedata(lkp::LookUpCld, ibnd, ice_rgh)
    n = get_nsize_ice(lkp)
    icedata = lkp.icedata
    return @inbounds (
        # LUT ice extinction coefficient (`nsize_ice, nband, nrghice`) m²/g
        view(icedata, 1:n, ibnd, ice_rgh),
        # LUT ice single scattering albedo (`nsize_ice, nband, nrghice`)
        view(icedata, (n + 1):(n * 2), ibnd, ice_rgh),
        # LUT ice asymmetry parameter (`nsize_ice, nband, nrghice`)
        view(icedata, (n * 2 + 1):(n * 3), ibnd, ice_rgh),
    )
end

"""
    PadeCld{FTA1D, FTA2D, FTA3D, FTA4D} <: AbstractLookUp

Pade coefficients for cloud optics.

This struct stores the Pade coefficients for determing extinction coeffient, 
single-scattering albedo, and asymmetry parameter g as a function of effective radius.
We compute the optical depth tau (=exintinction coeff * condensed water path)
and the products tau*ssa and tau*ssa*g for liquid and ice cloud separately.
These are used to determine the optical properties of ice and water cloud together.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct PadeCld{FTA1D, FTA2D, FTA3D, FTA4D} <: AbstractLookUp
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
end
Adapt.@adapt_structure PadeCld

# number of bands for Pade approximation
get_nband(lkp::PadeCld) = size(lkp.pade_extliq, 1)
# number of size regimes for Pade approximation
get_nsizereg(lkp::PadeCld) = size(lkp.pade_extliq, 2)
# number of extinction coefficients for pade approximation
get_ncoeff_ext(lkp::PadeCld) = size(lkp.pade_extliq, 3)
# number of ssa/g coefficients for pade approximation
get_ncoeff_ssa_g(lkp::PadeCld) = size(lkp.pade_ssaliq, 3)
# number of size regime boundaries for pade interpolation
get_nbound(lkp::PadeCld) = size(lkp.pade_sizreg_extliq, 1)

function PadeCld(ds, ::Type{FT}, ::Type{DA}) where {FT <: AbstractFloat, DA}
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    FTA3D = DA{FT, 3}
    FTA4D = DA{FT, 4}

    nsizereg = Int(ds.dim["nsizereg"])
    @assert nsizereg == 3 # RRTMGP pade approximation assumes exactly (3) size regimes

    pade_extliq = FTA3D(Array(ds["pade_extliq"]))
    pade_ssaliq = FTA3D(Array(ds["pade_ssaliq"]))
    pade_asyliq = FTA3D(Array(ds["pade_asyliq"]))

    pade_extice = FTA4D(Array(ds["pade_extice"]))
    pade_ssaice = FTA4D(Array(ds["pade_ssaice"]))
    pade_asyice = FTA4D(Array(ds["pade_asyice"]))

    pade_sizreg_extliq = FTA1D(Array(ds["pade_sizreg_extliq"]))
    pade_sizreg_ssaliq = FTA1D(Array(ds["pade_sizreg_ssaliq"]))
    pade_sizreg_asyliq = FTA1D(Array(ds["pade_sizreg_asyliq"]))

    pade_sizreg_extice = FTA1D(Array(ds["pade_sizreg_extice"]))
    pade_sizreg_ssaice = FTA1D(Array(ds["pade_sizreg_ssaice"]))
    pade_sizreg_asyice = FTA1D(Array(ds["pade_sizreg_asyice"]))

    bnd_lims_wn = FTA2D(Array(ds["bnd_limits_wavenumber"]))

    return (PadeCld{FTA1D, FTA2D, FTA3D, FTA4D}(
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
    ))
end

end
