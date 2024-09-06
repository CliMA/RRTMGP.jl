"""
    add_aerosol_optics_1scalar!(
        τ,
        aero_mask,
        aero_size,
        aero_mass,
        rel_hum,
        lkp_aero,
        ibnd,
    )

This function computes the `OneScalar` aerosol optics properties and adds them
to the exising `OneScalar` optics properties.
"""
@inline function add_aerosol_optics_1scalar!(τ, aero_mask, aero_size, aero_mass, rel_hum, lkp_aero, ibnd)
    FT = eltype(τ)
    @inbounds begin
        nlay = length(τ)
        for glay in 1:nlay
            if aero_mask[glay]
                τ_aero, τ_ssa_aero, τ_ssag_aero =
                    compute_lookup_aerosol(lkp_aero, ibnd, aero_mass, aero_size, rel_hum[glay], glay)
                τ[glay] += (τ_aero - τ_ssa_aero)
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
        aero_mask,
        aero_size,
        aero_mass,
        rel_hum,
        lkp_aero,
        ibnd;
        delta_scaling = false,
    )

This function computes the `TwoStream` aerosol optics properties and adds them
to the exising `TwoStream` optics properties.
"""
@inline function add_aerosol_optics_2stream!(
    τ,
    ssa,
    g,
    aero_mask,
    aero_size,
    aero_mass,
    rel_hum,
    lkp_aero,
    ibnd;
    delta_scaling = false,
)
    @inbounds begin
        nlay = length(τ)
        FT = eltype(τ)
        for glay in 1:nlay
            if aero_mask[glay]
                τ_aero, τ_ssa_aero, τ_ssag_aero =
                    compute_lookup_aerosol(lkp_aero, ibnd, aero_mass, aero_size, rel_hum[glay], glay)
                g_aero = τ_ssag_aero / max(eps(FT), τ_ssa_aero)
                ssa_aero = τ_ssa_aero / max(eps(FT), τ_aero)
                if delta_scaling # delta scaling is applied for shortwave problem
                    τ_aero, ssa_aero, g_aero = delta_scale(τ_aero, ssa_aero, g_aero)
                end
                τ[glay], ssa[glay], g[glay] = increment_2stream(τ[glay], ssa[glay], g[glay], τ_aero, ssa_aero, g_aero)
            end
        end
    end
    return nothing
end

"""
    compute_lookup_aerosol(lkp_aero, ibnd::Int, aero_mass, aero_size, rh::FT, glay) where {FT}

Compute cumulative aerosol optical properties for various aerosol particles given the aerosol lookup table (`lkp_aero`), 
band number (`ibnd`), aerosol mass (`aeromass`), aerosol size (`aerosize`), and relative humidity (`rh`).
"""
function compute_lookup_aerosol(lkp_aero, ibnd::Int, aero_mass, aero_size, rh::FT, glay) where {FT}
    τ_cum, τ_ssa_cum, τ_ssag_cum = FT(0), FT(0), FT(0)

    if aero_mass[1, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_dust_props(lkp_aero, ibnd, aero_mass[1, glay], aero_size[1, glay])
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[2, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_sea_salt_props(lkp_aero, ibnd, aero_mass[2, glay], aero_size[2, glay], rh)
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[3, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_sulfate_props(lkp_aero, ibnd, aero_mass[3, glay], rh)
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[4, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_black_carbon_props(lkp_aero, ibnd, aero_mass[4, glay], rh)
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[5, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_black_carbon_props(lkp_aero, ibnd, aero_mass[5, glay])
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[6, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_organic_carbon_props(lkp_aero, ibnd, aero_mass[6, glay], rh)
        τ_cum += τ
        τ_ssa_cum += τ_ssa
        τ_ssag_cum += τ_ssag
    end

    if aero_mass[7, glay] > FT(0)
        τ, τ_ssa, τ_ssag = compute_lookup_organic_carbon_props(lkp_aero, ibnd, aero_mass[7, glay])
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
function compute_lookup_dust_props(lkp_aero, ibnd, aeromass::FT, aerosize::FT) where {FT}
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
    compute_lookup_sea_salt_props(lkp_aero, ibnd, aeromass::FT, aerosize::FT, rh::FT) where {FT}

Compute aerosol optical properties for sea salt particles given the aerosol lookup table (`lkp_aero`), 
band number (`ibnd`), aerosol mass (`aeromass`), aerosol size (`aerosize`), and relative humidity (`rh`).
"""
function compute_lookup_sea_salt_props(lkp_aero, ibnd, aeromass::FT, aerosize::FT, rh::FT) where {FT}
    (; size_bin_limits, rh_levels, sea_salt) = lkp_aero
    @inbounds begin
        bin = locate_merra_size_bin(size_bin_limits, aerosize)
        loc, factor = interp1d_loc_factor(rh, rh_levels)
        τ = aeromass * (sea_salt[1, loc, bin, ibnd] * (FT(1) - factor) + sea_salt[1, loc + 1, bin, ibnd] * factor)
        τ_ssa = τ * (sea_salt[2, loc, bin, ibnd] * (FT(1) - factor) + sea_salt[2, loc + 1, bin, ibnd] * factor)
        τ_ssag = τ_ssa * (sea_salt[3, loc, bin, ibnd] * (FT(1) - factor) + sea_salt[3, loc + 1, bin, ibnd] * factor)
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_sulfate_props(lkp_aero, ibnd, aeromass::FT, rh::FT) where {FT}

Compute aerosol optical properties for sulfate particles given the aerosol lookup table (`lkp_aero`), 
band number (`ibnd`), aerosol mass (`aeromass`), and relative humidity (`rh`).
"""
function compute_lookup_sulfate_props(lkp_aero, ibnd, aeromass::FT, rh::FT) where {FT}
    (; rh_levels, sulfate) = lkp_aero
    @inbounds begin
        loc, factor = interp1d_loc_factor(rh, rh_levels)
        τ = aeromass * (sulfate[1, loc, ibnd] * (FT(1) - factor) + sulfate[1, loc + 1, ibnd] * factor)
        τ_ssa = τ * (sulfate[2, loc, ibnd] * (FT(1) - factor) + sulfate[2, loc + 1, ibnd] * factor)
        τ_ssag = τ_ssa * (sulfate[3, loc, ibnd] * (FT(1) - factor) + sulfate[3, loc + 1, ibnd] * factor)
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_black_carbon_props(lkp_aero, ibnd, aeromass::FT) where {FT}

Compute aerosol optical properties for hydrophobic black carbon particles given the aerosol lookup table (`lkp_aero`), 
band number (`ibnd`), and aerosol mass (`aeromass`).
"""
function compute_lookup_black_carbon_props(lkp_aero, ibnd, aeromass::FT) where {FT}
    (; black_carbon) = lkp_aero
    @inbounds begin
        τ = aeromass * black_carbon[1, ibnd]
        τ_ssa = τ * black_carbon[2, ibnd]
        τ_ssag = τ_ssa * black_carbon[3, ibnd]
    end
    return τ, τ_ssa, τ_ssag
end

"""
    compute_lookup_black_carbon_props(lkp_aero, ibnd, aeromass::FT, rh::FT) where {FT}

Compute aerosol optical properties for hydrophilic black carbon particles given the aerosol lookup table (`lkp_aero`), 
band number (`ibnd`), aerosol mass (`aeromass`), and relative humidity (`rh`).
"""
function compute_lookup_black_carbon_props(lkp_aero, ibnd, aeromass::FT, rh::FT) where {FT}
    (; rh_levels, black_carbon_rh) = lkp_aero
    @inbounds begin
        loc, factor = interp1d_loc_factor(rh, rh_levels)
        τ = aeromass * (black_carbon_rh[1, loc, ibnd] * (FT(1) - factor) + black_carbon_rh[1, loc + 1, ibnd] * factor)
        τ_ssa = τ * (black_carbon_rh[2, loc, ibnd] * (FT(1) - factor) + black_carbon_rh[2, loc + 1, ibnd] * factor)
        τ_ssag = τ_ssa * (black_carbon_rh[3, loc, ibnd] * (FT(1) - factor) + black_carbon_rh[3, loc + 1, ibnd] * factor)
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
    compute_lookup_organic_carbon_props(lkp_aero, ibnd, aeromass::FT, rh::FT) where {FT}

Compute aerosol optical properties for hydrophilic organic carbon particles given the aerosol lookup table (`lkp_aero`), 
band number (`ibnd`), aerosol mass (`aeromass`), and relative humidity (`rh`).
"""
function compute_lookup_organic_carbon_props(lkp_aero, ibnd, aeromass::FT, rh::FT) where {FT}
    (; rh_levels, organic_carbon_rh) = lkp_aero
    @inbounds begin
        loc, factor = interp1d_loc_factor(rh, rh_levels)
        τ =
            aeromass *
            (organic_carbon_rh[1, loc, ibnd] * (FT(1) - factor) + organic_carbon_rh[1, loc + 1, ibnd] * factor)
        τ_ssa = τ * (organic_carbon_rh[2, loc, ibnd] * (FT(1) - factor) + organic_carbon_rh[2, loc + 1, ibnd] * factor)
        τ_ssag =
            τ_ssa * (organic_carbon_rh[3, loc, ibnd] * (FT(1) - factor) + organic_carbon_rh[3, loc + 1, ibnd] * factor)
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

function compute_aero_mask!(aeromask::AbstractArray{B, 1}, aeromass::AbstractArray{FT, 2}) where {B, FT}
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
