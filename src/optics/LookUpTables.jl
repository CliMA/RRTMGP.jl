module LookUpTables

using DocStringExtensions
using Adapt

export AbstractLookUp,
    LookUpLW, LookUpSW, LookUpCld, LookUpAerosolMerra, LookUpMinor, ReferencePoints, LookUpPlanck, BandData

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
    LookUpAerosolMerra{D, D1, D2, D3, D4, W} <: AbstractLookUp

Merra lookup table for aersols. 

This struct stores the lookup tables for determing extinction coeffient,
single-scattering albedo, and asymmetry parameter g as a function of aerosol
particle size, relative humidity and band. Data is provided for dust, sea salt,
sulfate, black carbon (hydrophobic and hydrophilic) and organic carbon 
(hydrophobic and hydrophilic).


# Fields
$(DocStringExtensions.FIELDS)
"""
struct LookUpAerosolMerra{D, D1, D2, D3, D4, W} <: AbstractLookUp
    "`nband`, `nval`, `nbin`, `nrh`, `pair`"
    dims::D
    "beginning and ending limit for each MERRA aerosol size bin (microns)"
    size_bin_limits::D2
    "relative humidity levels for MERRA hydrophilic aerosols"
    rh_levels::D1
    "dust `(nval, nbin, nband)`"
    dust::D3
    "sea salt `(nval, nrh, nbin, nband)`"
    sea_salt::D4
    "sulfate `(nval, nrh, nband)`"
    sulfate::D3
    "black carbon - hydrophilic `(nval, nhr, nband)`"
    black_carbon_rh::D3
    "black carbon - hydrophobic `(nval, nband)`"
    black_carbon::D2
    "organic carbon - hydrophilic `(nval, nhr, nband)`"
    organic_carbon_rh::D3
    "organic carbon - hydrophobic `(nval, nband)`"
    organic_carbon::D2
    "beginning and ending wavenumber for each band (`2, nband`) cm⁻¹"
    bnd_lims_wn::W
    "Band number index corresponding to 550 nm"
    iband_550nm::Int
end
Adapt.@adapt_structure LookUpAerosolMerra

end
