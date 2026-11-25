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
                τ_aero, τ_ssa_aero, τ_ssag_aero = compute_lookup_aerosol(
                    lkp_aero,
                    ibnd,
                    aero_mass,
                    aero_size,
                    rel_hum[glay],
                    glay,
                )
                g_aero = τ_ssag_aero / max(eps(FT), τ_ssa_aero)
                ssa_aero = τ_ssa_aero / max(eps(FT), τ_aero)
                if collect_aod
                    aod_ext[] += τ_aero
                    aod_sca[] += τ_ssa_aero
                end
                if delta_scaling # delta scaling is applied for shortwave problem
                    τ_aero, ssa_aero, g_aero =
                        delta_scale(τ_aero, ssa_aero, g_aero)
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

"""
    compute_aerosol_band_optics!(
        τ_aero, ssa_aero, g_aero, aod_ext, aod_sca, aero_mask, aero_size,
        aero_mass, rel_hum, lkp_aero, ibnd, iband_550nm; delta_scaling = false,
    )

Compute one column's `TwoStream` aerosol optics for band `ibnd` into the
per-layer stores `τ_aero`, `ssa_aero` and `g_aero`, without adding them to any
gas optics.

Aerosol optics depend on the band, the column and the layer but not on the
g-point, so a solve that keeps them per band computes them once and adds them to
each of the band's g-points with [`add_aerosol_band_optics!`](@ref). The values
stored are the ones [`add_aerosol_optics_2stream!`](@ref) hands to
`increment_2stream`, so the two routes give the same optics bit for bit.

A layer the mask leaves out is not written: the reader tests the same mask.
"""
@inline function compute_aerosol_band_optics!(
    τ_store,
    ssa_store,
    g_store,
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
        nlay = length(τ_store)
        FT = eltype(τ_store)
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
                g_aero = τ_ssag_aero / max(eps(FT), τ_ssa_aero)
                ssa_aero = τ_ssa_aero / max(eps(FT), τ_aero)
                if collect_aod
                    aod_ext[] += τ_aero
                    aod_sca[] += τ_ssa_aero
                end
                if delta_scaling # delta scaling is applied for shortwave problem
                    τ_aero, ssa_aero, g_aero =
                        delta_scale(τ_aero, ssa_aero, g_aero)
                end
                τ_store[glay], ssa_store[glay], g_store[glay] =
                    τ_aero, ssa_aero, g_aero
            end
        end
    end
    return nothing
end

"""
    add_aerosol_band_optics!(τ, ssa, g, aero_mask, τ_aero, ssa_aero, g_aero)

Add one column's stored band aerosol optics, as
[`compute_aerosol_band_optics!`](@ref) left them, to the `TwoStream` optics
`τ`, `ssa` and `g` of one g-point.
"""
@inline function add_aerosol_band_optics!(
    τ,
    ssa,
    g,
    aero_mask,
    τ_store,
    ssa_store,
    g_store,
)
    @inbounds for glay in 1:length(τ)
        if aero_mask[glay]
            τ[glay], ssa[glay], g[glay] = increment_2stream(
                τ[glay],
                ssa[glay],
                g[glay],
                τ_store[glay],
                ssa_store[glay],
                g_store[glay],
            )
        end
    end
    return nothing
end

"""
    compute_lookup_aerosol(lkp_aero, ibnd::Int, aero_mass, aero_size, rh::FT, glay) where {FT}

Compute cumulative aerosol optical properties for various aerosol particles given the aerosol lookup table (`lkp_aero`),
band number (`ibnd`), aerosol mass (`aeromass`), aerosol size (`aerosize`), and relative humidity (`rh`).
"""
function compute_lookup_aerosol(
    lkp_aero,
    ibnd::Int,
    aero_mass,
    aero_size,
    rh::FT,
    glay,
) where {FT}
    τ_cum, τ_ssa_cum, τ_ssag_cum = FT(0), FT(0), FT(0)
    # Every humidity-dependent species interpolates the same `rh` over the same
    # levels, so the search runs once here rather than once per species.
    @inbounds loc, factor = interp1d_loc_factor(rh, lkp_aero.rh_levels)

    for i_aero_size in (1, 8, 9, 10, 11)
        if aero_mass[i_aero_size, glay] > FT(0)
            τ, τ_ssa, τ_ssag = compute_lookup_dust_props(
                lkp_aero,
                ibnd,
                aero_mass[i_aero_size, glay],
                aero_size[i_aero_size, glay],
            )
            τ_cum += τ
            τ_ssa_cum += τ_ssa
            τ_ssag_cum += τ_ssag
        end
    end

    for i_aero_size in (2, 12, 13, 14, 15)
        if aero_mass[i_aero_size, glay] > FT(0)
            τ, τ_ssa, τ_ssag = compute_lookup_sea_salt_props(
                lkp_aero,
                ibnd,
                aero_mass[i_aero_size, glay],
                aero_size[i_aero_size, glay],
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
    (; size_bin_limits, dust) = lkp_aero
    bin = locate_merra_size_bin(size_bin_limits, aerosize)
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
    (; size_bin_limits, sea_salt) = lkp_aero
    @inbounds begin
        bin = locate_merra_size_bin(size_bin_limits, aerosize)
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
