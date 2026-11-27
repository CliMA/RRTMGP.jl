"""
    add_aerosol_optics_1scalar!(
        τ,
        aod_ext,
        aod_sca,
        aero_mask,
        aero_size,
        aero_mass,
        rel_hum,
        lkp_aero,
        ibnd,
        iband_550nm
    )

Compute the `OneScalar` aerosol optics properties and add them
to the existing `OneScalar` optics properties.
"""
@inline function add_aerosol_optics_1scalar!(
    τ,
    aod_ext,
    aod_sca,
    aero_mask,
    aero_size,
    aero_mass,
    rel_hum,
    lkp_aero,
    ibnd,
    iband_550nm,
)
    FT = eltype(τ)
    @inbounds begin
        collect_aod = !isnothing(aod_ext) && (ibnd == iband_550nm)
        nlay = length(τ)
        collect_aod && (aod_ext[] = 0)
        collect_aod && (aod_sca[] = 0)
        for glay in 1:nlay
            if aero_mask[glay]
                τ_aero, τ_ssa_aero, τ_ssag_aero = compute_lookup_aerosol(
                    lkp_aero,
                    ibnd,
                    aero_mass,
                    aero_size,
                    rel_hum[glay],
                    glay,
                )
                τ[glay] += (τ_aero - τ_ssa_aero)
                if collect_aod
                    aod_ext[] += τ_aero
                    aod_sca[] += τ_ssa_aero
                end
            end
        end
    end
    return nothing
end

"""
    add_aerosol_optics_2stream!(
        τ,
        ssa,
        g,
        aod_ext,
        aod_sca,
        aero_mask,
        aero_size,
        aero_mass,
        rel_hum,
        lkp_aero,
        ibnd,
        iband_550nm;
        delta_scaling = false,
    )

Compute the `TwoStream` aerosol optics properties and add them
to the existing `TwoStream` optics properties:

 - `aod_ext` total aerosol optical depth
 - `aod_sca` scattering component of aerosol optical depth
"""
@inline function add_aerosol_optics_2stream!(
    τ,
    ssa,
    g,
    aod_ext,
    aod_sca,
    aero_mask,
    aero_size,
    aero_mass,
    rel_hum,
    lkp_aero,
    ibnd,
    iband_550nm;
    delta_scaling = false,
)
    @inbounds begin
        collect_aod = !isnothing(aod_ext) && (ibnd == iband_550nm)
        nlay = length(τ)
        FT = eltype(τ)
        collect_aod && (aod_ext[] = 0)
        collect_aod && (aod_sca[] = 0)
        for glay in 1:nlay
            if aero_mask[glay]
                (τ_aero, ssa_aero, g_aero), (τ_ext, τ_sca) =
                    aerosol_optics_2stream(
                        lkp_aero,
                        ibnd,
                        aero_mass,
                        aero_size,
                        rel_hum[glay],
                        glay,
                        delta_scaling,
                    )
                if collect_aod
                    aod_ext[] += τ_ext
                    aod_sca[] += τ_sca
                end
                τ[glay], ssa[glay], g[glay] = increment_2stream(
                    τ[glay],
                    ssa[glay],
                    g[glay],
                    τ_aero,
                    ssa_aero,
                    g_aero,
                )
            end
        end
    end
    return nothing
end

# The `TwoStream` aerosol optics of layer `glay` in band `ibnd`, delta-scaled for
# the shortwave, and the layer's optical depth and scattering optical depth.
@inline function aerosol_optics_2stream(
    lkp_aero,
    ibnd,
    aero_mass,
    aero_size,
    rel_hum,
    glay,
    delta_scaling,
)
    return aerosol_optics_2stream(
        lkp_aero,
        ibnd,
        aero_mass,
        glay,
        aerosol_searches(lkp_aero, aero_mass, aero_size, rel_hum, glay),
        delta_scaling,
    )
end
@inline function aerosol_optics_2stream(
    lkp_aero,
    ibnd,
    aero_mass,
    glay,
    (loc, factor, bins)::Tuple,
    delta_scaling,
)
    FT = typeof(factor)
    τ_aero, τ_ssa_aero, τ_ssag_aero = compute_lookup_aerosol(
        lkp_aero,
        ibnd,
        aero_mass,
        glay,
        loc,
        factor,
        bins,
    )
    g_aero = τ_ssag_aero / max(eps(FT), τ_ssa_aero)
    ssa_aero = τ_ssa_aero / max(eps(FT), τ_aero)
    optics =
        delta_scaling ? delta_scale(τ_aero, ssa_aero, g_aero) :
        (τ_aero, ssa_aero, g_aero)
    return optics, (τ_aero, τ_ssa_aero)
end

"""
    compute_lookup_aerosol(lkp_aero, ibnd::Integer, aero_mass, aero_size, rh::FT, glay) where {FT}

Compute cumulative aerosol optical properties for various aerosol particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), aerosol mass (`aeromass`), aerosol size (`aerosize`), and relative humidity (`rh`).
"""
function compute_lookup_aerosol(
    lkp_aero,
    ibnd::Integer,
    aero_mass,
    aero_size,
    rh::FT,
    glay,
) where {FT}
    return compute_lookup_aerosol(
        lkp_aero,
        ibnd,
        aero_mass,
        glay,
        aerosol_searches(lkp_aero, aero_mass, aero_size, rh, glay)...,
    )
end

# The dust and sea salt species, whose optics depend on their size bin
const DUST_SPECIES = (1, 8, 9, 10, 11)
const SEA_SALT_SPECIES = (2, 12, 13, 14, 15)

"""
    aerosol_searches(lkp_aero, aero_mass, aero_size, rh, glay)

The parts of the aerosol optics of layer `glay` that do not depend on the band:
the humidity interpolation location and factor, which every humidity-dependent
species shares, and the size bins of the dust and sea salt species with mass,
three bits each, dust first.
"""
@inline function aerosol_searches(lkp_aero, aero_mass, aero_size, rh, glay)
    FT = typeof(rh)
    @inbounds loc, factor = interp1d_loc_factor(rh, lkp_aero.rh_levels)
    bins = UInt32(0)
    for (k, i_aero_size) in enumerate((DUST_SPECIES..., SEA_SALT_SPECIES...))
        @inbounds if aero_mass[i_aero_size, glay] > FT(0)
            bin = locate_merra_size_bin(
                lkp_aero.size_bin_limits,
                aero_size[i_aero_size, glay],
            )
            bins |= UInt32(bin) << (3 * (k - 1))
        end
    end
    return loc, factor, bins
end

# the size bin of the `k`-th of the dust and sea salt species
@inline size_bin(bins, k) = Int((bins >> (3 * (k - 1))) & 0x7)

"""
    compute_lookup_aerosol(lkp_aero, ibnd, aero_mass, glay, loc, factor, bins)

The aerosol optics of layer `glay` in band `ibnd`, from the searches of
[`aerosol_searches`](@ref).
"""
function compute_lookup_aerosol(
    lkp_aero,
    ibnd::Integer,
    aero_mass,
    glay,
    loc,
    factor::FT,
    bins::UInt32,
) where {FT}
    τ_cum, τ_ssa_cum, τ_ssag_cum = FT(0), FT(0), FT(0)

    for (k, i_aero_size) in enumerate(DUST_SPECIES)
        if aero_mass[i_aero_size, glay] > FT(0)
            τ, τ_ssa, τ_ssag = dust_props(
                lkp_aero,
                ibnd,
                aero_mass[i_aero_size, glay],
                size_bin(bins, k),
            )
            τ_cum += τ
            τ_ssa_cum += τ_ssa
            τ_ssag_cum += τ_ssag
        end
    end

    for (k, i_aero_size) in enumerate(SEA_SALT_SPECIES)
        if aero_mass[i_aero_size, glay] > FT(0)
            τ, τ_ssa, τ_ssag = sea_salt_props(
                lkp_aero,
                ibnd,
                aero_mass[i_aero_size, glay],
                size_bin(bins, length(DUST_SPECIES) + k),
                loc,
                factor,
            )
            τ_cum += τ
            τ_ssa_cum += τ_ssa
            τ_ssag_cum += τ_ssag
        end
    end

    if aero_mass[3, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_sulfate_props(
            lkp_aero,
            ibnd,
            aero_mass[3, glay],
            loc,
            factor,
        )
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[4, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_black_carbon_props(
            lkp_aero,
            ibnd,
            aero_mass[4, glay],
            loc,
            factor,
        )
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[5, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_black_carbon_props(
            lkp_aero,
            ibnd,
            aero_mass[5, glay],
        )
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[6, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_organic_carbon_props(
            lkp_aero,
            ibnd,
            aero_mass[6, glay],
            loc,
            factor,
        )
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[7, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_organic_carbon_props(
            lkp_aero,
            ibnd,
            aero_mass[7, glay],
        )
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    return τ_cum, τ_ssa_cum, τ_ssag_cum
end

"""
    compute_lookup_dust_props(lkp_aero, ibnd, aeromass::FT, aerosize::FT) where {FT}

Compute aerosol optical properties for dust particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), aerosol mass (`aeromass`) and aerosol size (`aerosize`).
"""
function compute_lookup_dust_props(
    lkp_aero,
    ibnd,
    aeromass::FT,
    aerosize::FT,
) where {FT}
    bin = locate_merra_size_bin(lkp_aero.size_bin_limits, aerosize)
    return dust_props(lkp_aero, ibnd, aeromass, bin)
end

# the dust optics of size bin `bin` in band `ibnd`
@inline function dust_props(lkp_aero, ibnd, aeromass, bin)
    (; dust) = lkp_aero
    @inbounds begin
        τ = aeromass * dust[1, bin, ibnd]
        τ_ssa = τ * dust[2, bin, ibnd]
        τ_ssag = τ_ssa * dust[3, bin, ibnd]
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_sea_salt_props(lkp_aero, ibnd, aeromass::FT, aerosize::FT, loc, factor::FT) where {FT}

Compute aerosol optical properties for sea salt particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), aerosol mass (`aeromass`), aerosol size (`aerosize`), and the humidity
interpolation weights (`loc`, `factor`) from `interp1d_loc_factor`.
"""
function compute_lookup_sea_salt_props(
    lkp_aero,
    ibnd,
    aeromass::FT,
    aerosize::FT,
    loc,
    factor::FT,
) where {FT}
    bin = locate_merra_size_bin(lkp_aero.size_bin_limits, aerosize)
    return sea_salt_props(lkp_aero, ibnd, aeromass, bin, loc, factor)
end

# the sea salt optics of size bin `bin` in band `ibnd`, at the humidity
# interpolation location `loc` and factor `factor`
@inline function sea_salt_props(lkp_aero, ibnd, aeromass, bin, loc, factor)
    FT = typeof(aeromass)
    (; sea_salt) = lkp_aero
    @inbounds begin
        τ =
            aeromass * (
                sea_salt[1, loc, bin, ibnd] * (FT(1) - factor) +
                sea_salt[1, loc + 1, bin, ibnd] * factor
            )
        τ_ssa =
            τ * (
                sea_salt[2, loc, bin, ibnd] * (FT(1) - factor) +
                sea_salt[2, loc + 1, bin, ibnd] * factor
            )
        τ_ssag =
            τ_ssa * (
                sea_salt[3, loc, bin, ibnd] * (FT(1) - factor) +
                sea_salt[3, loc + 1, bin, ibnd] * factor
            )
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_sulfate_props(lkp_aero, ibnd, aeromass::FT, loc, factor::FT) where {FT}

Compute aerosol optical properties for sulfate particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), aerosol mass (`aeromass`), and the humidity interpolation weights
(`loc`, `factor`) from `interp1d_loc_factor`.
"""
function compute_lookup_sulfate_props(
    lkp_aero,
    ibnd,
    aeromass::FT,
    loc,
    factor::FT,
) where {FT}
    (; sulfate) = lkp_aero
    @inbounds begin
        τ =
            aeromass * (
                sulfate[1, loc, ibnd] * (FT(1) - factor) +
                sulfate[1, loc + 1, ibnd] * factor
            )
        τ_ssa =
            τ * (
                sulfate[2, loc, ibnd] * (FT(1) - factor) +
                sulfate[2, loc + 1, ibnd] * factor
            )
        τ_ssag =
            τ_ssa * (
                sulfate[3, loc, ibnd] * (FT(1) - factor) +
                sulfate[3, loc + 1, ibnd] * factor
            )
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_black_carbon_props(lkp_aero, ibnd, aeromass::FT) where {FT}

Compute aerosol optical properties for hydrophobic black carbon particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), and aerosol mass (`aeromass`).
"""
function compute_lookup_black_carbon_props(
    lkp_aero,
    ibnd,
    aeromass::FT,
) where {FT}
    (; black_carbon) = lkp_aero
    @inbounds begin
        τ = aeromass * black_carbon[1, ibnd]
        τ_ssa = τ * black_carbon[2, ibnd]
        τ_ssag = τ_ssa * black_carbon[3, ibnd]
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_black_carbon_props(lkp_aero, ibnd, aeromass::FT, loc, factor::FT) where {FT}

Compute aerosol optical properties for hydrophilic black carbon particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), aerosol mass (`aeromass`), and the humidity interpolation weights
(`loc`, `factor`) from `interp1d_loc_factor`.
"""
function compute_lookup_black_carbon_props(
    lkp_aero,
    ibnd,
    aeromass::FT,
    loc,
    factor::FT,
) where {FT}
    (; black_carbon_rh) = lkp_aero
    @inbounds begin
        τ =
            aeromass * (
                black_carbon_rh[1, loc, ibnd] * (FT(1) - factor) +
                black_carbon_rh[1, loc + 1, ibnd] * factor
            )
        τ_ssa =
            τ * (
                black_carbon_rh[2, loc, ibnd] * (FT(1) - factor) +
                black_carbon_rh[2, loc + 1, ibnd] * factor
            )
        τ_ssag =
            τ_ssa * (
                black_carbon_rh[3, loc, ibnd] * (FT(1) - factor) +
                black_carbon_rh[3, loc + 1, ibnd] * factor
            )
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_organic_carbon_props(lkp_aero, ibnd, aeromass)

Compute aerosol optical properties for hydrophobic organic carbon particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), and aerosol mass (`aeromass`).
"""
function compute_lookup_organic_carbon_props(lkp_aero, ibnd, aeromass)
    (; organic_carbon) = lkp_aero
    @inbounds begin
        τ = aeromass * organic_carbon[1, ibnd]
        τ_ssa = τ * organic_carbon[2, ibnd]
        τ_ssag = τ_ssa * organic_carbon[3, ibnd]
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_organic_carbon_props(lkp_aero, ibnd, aeromass::FT, loc, factor::FT) where {FT}

Compute aerosol optical properties for hydrophilic organic carbon particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), aerosol mass (`aeromass`), and the humidity interpolation weights
(`loc`, `factor`) from `interp1d_loc_factor`.
"""
function compute_lookup_organic_carbon_props(
    lkp_aero,
    ibnd,
    aeromass::FT,
    loc,
    factor::FT,
) where {FT}
    (; organic_carbon_rh) = lkp_aero
    @inbounds begin
        τ =
            aeromass * (
                organic_carbon_rh[1, loc, ibnd] * (FT(1) - factor) +
                organic_carbon_rh[1, loc + 1, ibnd] * factor
            )
        τ_ssa =
            τ * (
                organic_carbon_rh[2, loc, ibnd] * (FT(1) - factor) +
                organic_carbon_rh[2, loc + 1, ibnd] * factor
            )
        τ_ssag =
            τ_ssa * (
                organic_carbon_rh[3, loc, ibnd] * (FT(1) - factor) +
                organic_carbon_rh[3, loc + 1, ibnd] * factor
            )
    end
    return τ, τ_ssa, τ_ssag
end

"""
    locate_merra_size_bin(size_bin_limits, aerosize)

Locate merra bin number for a given aerosol size (`aerosize`).
"""
function locate_merra_size_bin(size_bin_limits, aerosize)
    nbins = size(size_bin_limits, 2)
    bin = 1
    @inbounds begin
        for ibin in 1:nbins
            if size_bin_limits[1, ibin] ≤ aerosize ≤ size_bin_limits[2, ibin]
                bin = ibin
                break
            end
            bin = nbins
        end
    end
    return bin
end

function _aero_mask(aero_mass_layer::AbstractArray{FT, 1}) where {FT}
    aero_mask = false
    for aeromass in aero_mass_layer
        if aeromass > FT(0)
            aero_mask = true
            break
        end
    end
    return aero_mask
end

function compute_aero_mask!(
    aeromask::AbstractArray{B, 1},
    aeromass::AbstractArray{FT, 2},
) where {B, FT}
    nlay = length(aeromask)
    naerosols = size(aeromass, 1)
    @inbounds begin
        for ilay in 1:nlay
            mask = false
            for iaero in 1:naerosols
                if aeromass[iaero, ilay] > zero(FT)
                    mask = true
                    break
                end
            end
            aeromask[ilay] = mask
        end
    end
    return nothing
end
