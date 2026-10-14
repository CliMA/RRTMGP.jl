module LookUpTables

using Adapt

export AbstractLookUp,
    LookUpLW,
    LookUpSW,
    LookUpCld,
    LookUpAerosolMerra,
    LookUpMinor,
    ReferencePoints,
    LookUpPlanck,
    BandData

"""
    AbstractLookUp

Abstract lookup table for longwave and shortwave problems.
"""
abstract type AbstractLookUp end

"""
    LookUpMinor{O, G, K}

Lookup table for computing optical properties of minor gases.

# Fields
- `bnd_st`: Starting index to `idx_gases_minor_lower/upper` for each band `(n_bnd + 1)`.
- `gpt_st`: Starting index in `kminor_lower/upper` for each g-point `(n_gpt + 1)`.
- `gasdata`: Indices for minor gases contributing to absorption in the lower/upper
  atmosphere, indices for scaling gases, whether the minor gas scales with density, and
  whether the minor gas scales by complement; `(4, n_min_absrb_lower/upper)`.
- `kminor`: Minor absorption coefficient in the lower/upper atmosphere
  `(n_η, n_t_ref, n_contrib_lower/upper)`.
"""
struct LookUpMinor{O, G, K}
    bnd_st::O
    gpt_st::O
    gasdata::G
    kminor::K
end
Adapt.@adapt_structure LookUpMinor

# number of minor absorbors in the lower/upper atmosphere
@inline get_n_min_absrb(lkp::LookUpMinor) = size(lkp.gasdata, 2)
# number of minor contributors in the lower/upper atmosphere
@inline get_n_contrib(lkp::LookUpMinor) = size(lkp.kminor, 3)
# get bounds for each g-point
@inline getbounds(lkp::LookUpMinor, ibnd, igpt) = @inbounds (
    lkp.bnd_st[ibnd],
    lkp.gpt_st[igpt],
    lkp.gpt_st[igpt + 1] - lkp.gpt_st[igpt],
)
# returns minor gas index, scaling gas index, minor gas scales with density and minor gas scales by complement at the given `idx`.
@inline get_minor_gas_data(lkp::LookUpMinor, idx) =
    @inbounds view(lkp.gasdata, :, idx)

"""
    ReferencePoints{RD1, RD3}

Log of reference pressures, reference temperatures, and volume mixing ratios used by the
lookup table.

# Fields
- `ln_p_ref`: Log of reference pressures used by the lookup table `(n_p_ref)`.
- `t_ref`: Reference temperatures used by the lookup table `(n_t_ref)`.
- `vmr_ref`: Reference volume mixing ratios used by the lookup table
  `(2, n_gases, n_t_ref)`.
"""
struct ReferencePoints{RD1, RD3}
    ln_p_ref::RD1
    t_ref::RD1
    vmr_ref::RD3
end
Adapt.@adapt_structure ReferencePoints

"""
    LookUpPlanck{PD1, PD2, PD4}

Lookup data for Planck source calculations.

# Fields
- `planck_fraction`: Planck fraction `(n_η, n_p_ref, n_t_ref, n_gpt)`.
- `t_planck`: Reference temperatures for Planck source calculations `(n_t_plnk)`.
- `tot_planck`: Total Planck source for each band `(n_t_plnk, n_bnd)`.
"""
struct LookUpPlanck{PD1, PD2, PD4}
    planck_fraction::PD4
    t_planck::PD1
    tot_planck::PD2
end
Adapt.@adapt_structure LookUpPlanck

"""
    BandData{BNDI1, BNDI2, BNDD2}

Band/g-point data for the lookup table.

# Fields
- `major_gpt2bnd`: Map from g-point to band.
- `bnd_lims_gpt`: Starting and ending g-point for each band `(2, n_bnd)`.
- `bnd_lims_wn`: Starting and ending wavenumber for each band `(2, n_bnd)`.
"""
struct BandData{BNDI1, BNDI2, BNDD2}
    major_gpt2bnd::BNDI1
    bnd_lims_gpt::BNDI2
    bnd_lims_wn::BNDD2
end
Adapt.@adapt_structure BandData

"""
    LookUpLW{FT, IA3D, FTA4D, BND, P, R, LMNR} <: AbstractLookUp

Longwave lookup tables, used to compute optical properties.

# Fields
- `idx_h2o`: Vmr array index for H₂O.
- `p_ref_tropo`: Reference pressure separating the upper and lower atmosphere.
- `p_ref_min`: Minimum pressure supported by RRTMGP lookup tables.
- `t_ref_min`: Minimum temperature supported by RRTMGP lookup tables (first reference temperature).
- `t_ref_max`: Maximum temperature supported by RRTMGP lookup tables (last reference temperature).
- `key_species`: Major absorbing species in each band `(2, n_atmos_layers, n_bnd)`.
- `kmajor`: Major absorption coefficient `(n_η, n_p_ref + 1, n_t_ref, n_gpt)`.
- `planck`: Lookup data for Planck source calculations.
- `band_data`: Band data.
- `ref_points`: Reference temperatures, pressures and volume mixing ratios.
- `minor_lower`: Lookup data for minor gases in the lower atmosphere.
- `minor_upper`: Lookup data for minor gases in the upper atmosphere.
"""
struct LookUpLW{FT, IA3D, FTA4D, BND, P, R, LMNR} <: AbstractLookUp
    idx_h2o::Int
    p_ref_tropo::FT
    p_ref_min::FT
    t_ref_min::FT
    t_ref_max::FT
    key_species::IA3D
    kmajor::FTA4D
    planck::P
    band_data::BND
    ref_points::R
    minor_lower::LMNR
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
- `idx_h2o`: Vmr array index for H₂O.
- `p_ref_tropo`: Reference pressure separating the upper and lower atmosphere.
- `p_ref_min`: Minimum pressure supported by RRTMGP lookup tables.
- `t_ref_min`: Minimum temperature supported by RRTMGP lookup tables (first reference temperature).
- `t_ref_max`: Maximum temperature supported by RRTMGP lookup tables (last reference temperature).
- `solar_src_tot`: Total solar irradiation.
- `key_species`: Major absorbing species in each band `(2, n_atmos_layers, n_bnd)`.
- `kmajor`: Major absorption coefficient `(n_η, n_p_ref + 1, n_t_ref, n_gpt)`.
- `band_data`: Band data.
- `ref_points`: Reference temperatures, pressures and volume mixing ratios.
- `rayl_lower`: Rayleigh absorption coefficient for the lower atmosphere
  `(n_η, n_t_ref, n_gpt)`.
- `rayl_upper`: Rayleigh absorption coefficient for the upper atmosphere
  `(n_η, n_t_ref, n_gpt)`.
- `solar_src_scaled`: Relative solar source contribution from each g-point `(n_gpt)`.
- `minor_lower`: Lookup data for minor gases in the lower atmosphere.
- `minor_upper`: Lookup data for minor gases in the upper atmosphere.
"""
struct LookUpSW{FT, IA3D, FTA1D, FTA3D, FTA4D, BND, R, LMNR} <: AbstractLookUp
    idx_h2o::Int
    p_ref_tropo::FT
    p_ref_min::FT
    t_ref_min::FT
    t_ref_max::FT
    solar_src_tot::FT
    key_species::IA3D
    kmajor::FTA4D
    band_data::BND
    ref_points::R
    rayl_lower::FTA3D
    rayl_upper::FTA3D
    solar_src_scaled::FTA1D
    minor_lower::LMNR
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
    LookUpCld{D, B, L, I, R, S, RB, SB, W} <: AbstractLookUp

Lookup table for cloud (and, optionally, precipitation) optics.

Stores extinction coefficient, single-scattering albedo, and asymmetry parameter
g as a function of effective radius, from which we form τ (= ext * water path)
and the products τ*ssa and τ*ssa*g for each condensate species. Liquid and ice
are always present; rain and snow are present when `raindata`/`snowdata` are not
`nothing`.

# Fields
- `dims`: Dimensions `nband`, `nrghice`, `nsize_liq`, `nsize_ice`, `pair`.
- `bounds`: Liquid/ice particle size lower and upper bounds and LUT interpolation
  factor: `(radliq_lwr, radliq_upr, radice_lwr, radice_upr)`.
- `liqdata`: Liquid ext, ssa, asymmetry `(3*nsize_liq, nbnd)`.
- `icedata`: Ice ext, ssa, asymmetry `(3*nsize_ice, nbnd, nrghice)`.
- `raindata`: Rain ext, ssa, asymmetry `(3*nsize_rai, nbnd)`, or `nothing`.
- `snowdata`: Snow ext, ssa, asymmetry `(3*nsize_sno, nbnd)`, or `nothing`.
- `rai_bounds`: `(radrai_lwr, radrai_upr)`, or `nothing`.
- `sno_bounds`: `(radsno_lwr, radsno_upr)`, or `nothing`.
- `bnd_lims_wn`: Beginning and ending wavenumber for each band [cm⁻¹] `(2, nband)`.
"""
struct LookUpCld{D, B, L, I, R, S, RB, SB, W} <: AbstractLookUp
    dims::D
    bounds::B
    liqdata::L
    icedata::I
    raindata::R
    snowdata::S
    rai_bounds::RB
    sno_bounds::SB
    bnd_lims_wn::W
end
Adapt.@adapt_structure LookUpCld

# Backward-compatible constructor: cloud-only (no precip) tables.
LookUpCld(dims, bounds, liqdata, icedata, bnd_lims_wn) =
    LookUpCld(dims, bounds, liqdata, icedata, nothing, nothing, nothing, nothing, bnd_lims_wn)

# True when rain/snow optics are available.
@inline has_precip_optics(lkp::LookUpCld) = !isnothing(lkp.raindata)

@inline get_dims(lkp::LookUpCld) = lkp.dims
@inline get_nband(lkp::LookUpCld) = @inbounds lkp.dims[1]
@inline get_nrghice(lkp::LookUpCld) = @inbounds lkp.dims[2]
@inline get_nsize_liq(lkp::LookUpCld) = @inbounds lkp.dims[3]
@inline get_nsize_ice(lkp::LookUpCld) = @inbounds lkp.dims[4]
@inline get_pair(lkp::LookUpCld) = @inbounds lkp.dims[5]
# Derived from the data arrays (ext/ssa/asy packed along dim 1), so no dims change.
@inline get_nsize_rai(lkp::LookUpCld) = size(lkp.raindata, 1) ÷ 3
@inline get_nsize_sno(lkp::LookUpCld) = size(lkp.snowdata, 1) ÷ 3

@inline function getview_liqdata(lkp::LookUpCld, ibnd)
    n = get_nsize_liq(lkp)
    liqdata = lkp.liqdata
    return @inbounds (
        view(liqdata, 1:n, ibnd),
        view(liqdata, (n + 1):(n * 2), ibnd),
        view(liqdata, (n * 2 + 1):(n * 3), ibnd),
    )
end

@inline function getview_icedata(lkp::LookUpCld, ibnd, ice_rgh)
    n = get_nsize_ice(lkp)
    icedata = lkp.icedata
    return @inbounds (
        view(icedata, 1:n, ibnd, ice_rgh),
        view(icedata, (n + 1):(n * 2), ibnd, ice_rgh),
        view(icedata, (n * 2 + 1):(n * 3), ibnd, ice_rgh),
    )
end

@inline function getview_raindata(lkp::LookUpCld, ibnd)
    n = get_nsize_rai(lkp)
    raindata = lkp.raindata
    return @inbounds (
        view(raindata, 1:n, ibnd),
        view(raindata, (n + 1):(n * 2), ibnd),
        view(raindata, (n * 2 + 1):(n * 3), ibnd),
    )
end

@inline function getview_snowdata(lkp::LookUpCld, ibnd)
    n = get_nsize_sno(lkp)
    snowdata = lkp.snowdata
    return @inbounds (
        view(snowdata, 1:n, ibnd),
        view(snowdata, (n + 1):(n * 2), ibnd),
        view(snowdata, (n * 2 + 1):(n * 3), ibnd),
    )
end
"""
    with_ice_snow_liquid_rain(lkp::LookUpCld; ice_rgh = 1) -> LookUpCld

STOPGAP / TEST ONLY. Returns a copy of `lkp` whose precip optics are aliased to
existing cloud tables: rain → liquid table (bounds = liquid's), snow → ice table
at roughness `ice_rgh` (bounds = ice's). Because rain/snow effective radii exceed
these bounds, every precip layer clamps to the largest tabulated particle — rain
gets liquid-at-max-radius optics, snow gets ice-at-max-radius optics. This
validates the precip plumbing end-to-end; it is NOT physically meaningful (no
precip-sized scattering, no size dependence). Replace with real rain/snow LUTs.
"""
function with_ice_snow_liquid_rain(lkp::LookUpCld; ice_rgh::Integer = 1)
    LookUpTables.has_precip_optics(lkp) && @warn "with_ice_snow_liquid_rain: table already has precip optics; overwriting with aliased stopgap data."

    DA = typeof(lkp.liqdata)
    BT = eltype(lkp.bounds)
    precip_bounds(lwr, upr) = DA(BT[BT(lwr), BT(upr)])

    radliq_lwr, radliq_upr, radice_lwr, radice_upr = lkp.bounds

    # rain → liquid table (already (3*nsize_liq, nbnd), no roughness dim)
    raindata = DA(copy(lkp.liqdata))

    # snow → ice table at one roughness level; icedata is
    # (3*nsize_ice, nbnd, nrghice), flatten to (3*nsize_ice, nbnd)
    snowdata = DA(lkp.icedata[:, :, ice_rgh])

    return LookUpCld(
        lkp.dims,
        lkp.bounds,
        lkp.liqdata,
        lkp.icedata,
        raindata,
        snowdata,
        precip_bounds(radliq_lwr, radliq_upr),   # rai_bounds = liquid bounds
        precip_bounds(radice_lwr, radice_upr),   # sno_bounds = ice bounds
        lkp.bnd_lims_wn,
    )
end
"""
    LookUpAerosolMerra{D, D1, D2, D3, D4, W} <: AbstractLookUp

MERRA lookup table for aerosols.

This struct stores the lookup tables for determining the extinction coefficient,
single-scattering albedo, and asymmetry parameter g as a function of aerosol
particle size, relative humidity and band. Data is provided for dust, sea salt,
sulfate, black carbon (hydrophobic and hydrophilic) and organic carbon
(hydrophobic and hydrophilic).

# Fields
- `dims`: Dimensions `nband`, `nval`, `nbin`, `nrh`, `pair`.
- `size_bin_limits`: Beginning and ending limit for each MERRA aerosol size bin
  [microns].
- `rh_levels`: Relative humidity levels for MERRA hydrophilic aerosols.
- `dust`: Dust `(nval, nbin, nband)`.
- `sea_salt`: Sea salt `(nval, nrh, nbin, nband)`.
- `sulfate`: Sulfate `(nval, nrh, nband)`.
- `black_carbon_rh`: Black carbon, hydrophilic `(nval, nrh, nband)`.
- `black_carbon`: Black carbon, hydrophobic `(nval, nband)`.
- `organic_carbon_rh`: Organic carbon, hydrophilic `(nval, nrh, nband)`.
- `organic_carbon`: Organic carbon, hydrophobic `(nval, nband)`.
- `bnd_lims_wn`: Beginning and ending wavenumber for each band [cm⁻¹] `(2, nband)`.
- `iband_550nm`: Band number index corresponding to 550 nm.
"""
struct LookUpAerosolMerra{D, D1, D2, D3, D4, W} <: AbstractLookUp
    dims::D
    size_bin_limits::D2
    rh_levels::D1
    dust::D3
    sea_salt::D4
    sulfate::D3
    black_carbon_rh::D3
    black_carbon::D2
    organic_carbon_rh::D3
    organic_carbon::D2
    bnd_lims_wn::W
    iband_550nm::Int
end
Adapt.@adapt_structure LookUpAerosolMerra

end
