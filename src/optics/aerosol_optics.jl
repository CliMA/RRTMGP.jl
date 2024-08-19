"""
    add_aerosol_optics_1scalar!(
        τ,
        aero_type,
        aero_size,
        aero_mass,
        rel_hum,
        lkp_aero::LookUpAerosolMerra,
        ibnd,
    )

This function computes the `OneScalar` aerosol optics properties and adds them
to the exising `OneScalar` optics properties.
"""
@inline function add_aerosol_optics_1scalar!(
    τ,
    aero_type,
    aero_size,
    aero_mass,
    rel_hum,
    lkp_aero::LookUpAerosolMerra,
    ibnd,
)
    FT = eltype(τ)
    @inbounds begin
        nlay = length(τ)
        for glay in 1:nlay
            aero_type_lay = aero_type[glay]
            if aero_type_lay ≠ 0
                τ_aero, τ_ssa_aero, τ_ssag_aero = compute_lookup_aerosol(
                    lkp_aero,
                    ibnd,
                    aero_type_lay,
                    aero_mass[glay],
                    aero_size[glay],
                    rel_hum[glay],
                )
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
        aero_type,
        aero_size,
        aero_mass,
        rel_hum,
        lkp_aero::LookUpAerosolMerra,
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
    aero_type,
    aero_size,
    aero_mass,
    rel_hum,
    lkp_aero::LookUpAerosolMerra,
    ibnd;
    delta_scaling = false,
)
    @inbounds begin
        nlay = length(τ)
        FT = eltype(τ)
        for glay in 1:nlay
            aero_type_lay = aero_type[glay]
            if aero_type_lay > 0
                τ_aero, τ_ssa_aero, τ_ssag_aero = compute_lookup_aerosol(
                    lkp_aero,
                    ibnd,
                    aero_type_lay,
                    aero_mass[glay],
                    aero_size[glay],
                    rel_hum[glay],
                )
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

function compute_lookup_aerosol(lkp_aero, ibnd::Int, aerotype::Int, aeromass::FT, aerosize::FT, rh::FT) where {FT}
    cargs = (lkp_aero, ibnd, aeromass) # common args
    return aerotype == 1 ? compute_lookup_dust_props(cargs..., aerosize) :
           (
        aerotype == 2 ? compute_lookup_sea_salt_props(cargs..., aerosize, rh) :
        (
            aerotype == 3 ? compute_lookup_sulfate_props(cargs..., rh) :
            (
                aerotype == 4 ? compute_lookup_black_carbon_props(cargs...) :
                (
                    aerotype == 5 ? compute_lookup_black_carbon_props(cargs..., rh) :
                    (
                        aerotype == 6 ? compute_lookup_organic_carbon_props(cargs...) :
                        (aerotype == 7 ? compute_lookup_organic_carbon_props(cargs..., rh) : (FT(0), FT(0), FT(0)))
                    )
                )
            )
        )
    )
end


function compute_lookup_dust_props(lkp_aero, ibnd::Int, aeromass::FT, aerosize::FT) where {FT}
    (; size_bin_limits, dust) = lkp_aero
    bin = locate_merra_size_bin(size_bin_limits, aerosize)
    τ = aeromass * dust[1, bin, ibnd]
    τ_ssa = τ * dust[2, bin, ibnd]
    τ_ssag = τ_ssa * dust[3, bin, ibnd]
    return τ, τ_ssa, τ_ssag
end

function compute_lookup_sea_salt_props(lkp_aero, ibnd::Int, aeromass::FT, aerosize::FT, rh::FT) where {FT}
    return (FT(0), FT(0), FT(0))
end

function compute_lookup_sulfate_props(lkp_aero, ibnd::Int, aeromass::FT, rh::FT) where {FT}
    (; rh_levels, sulfate) = lkp_aero
    @inbounds begin
        loc, factor = interp1d_loc_factor(rh, rh_levels)
        τ = aeromass * (sulfate[1, loc, ibnd] * (FT(1) - factor) + sulfate[1, loc + 1, ibnd] * factor)
        τ_ssa = τ * (sulfate[2, loc, ibnd] * (FT(1) - factor) + sulfate[2, loc + 1, ibnd] * factor)
        τ_ssag = τ_ssa * (sulfate[3, loc, ibnd] * (FT(1) - factor) + sulfate[3, loc + 1, ibnd] * factor)
    end
    return τ, τ_ssa, τ_ssag
end

function compute_lookup_black_carbon_props(lkp_aero, ibnd::Int, aeromass::FT) where {FT}
    (; black_carbon) = lkp_aero
    τ = aeromass * black_carbon[1, ibnd]
    τ_ssa = τ * black_carbon[2, ibnd]
    τ_ssag = τ_ssa * black_carbon[3, ibnd]
    return τ, τ_ssa, τ_ssag
end

function compute_lookup_black_carbon_props(lkp_aero, ibnd::Int, aeromass::FT, rh::FT) where {FT}
    (; rh_levels, black_carbon_rh) = lkp_aero
    @inbounds begin
        loc, factor = interp1d_loc_factor(rh, rh_levels)
        τ = aeromass * (black_carbon_rh[1, loc, ibnd] * (FT(1) - factor) + black_carbon_rh[1, loc + 1, ibnd] * factor)
        τ_ssa = τ * (black_carbon_rh[2, loc, ibnd] * (FT(1) - factor) + black_carbon_rh[2, loc + 1, ibnd] * factor)
        τ_ssag = τ_ssa * (black_carbon_rh[3, loc, ibnd] * (FT(1) - factor) + black_carbon_rh[3, loc + 1, ibnd] * factor)
    end
    return τ, τ_ssa, τ_ssag
end

function compute_lookup_organic_carbon_props(lkp_aero, ibnd::Int, aeromass::FT) where {FT}
    (; organic_carbon) = lkp_aero
    τ = aeromass * organic_carbon[1, ibnd]
    τ_ssa = τ * organic_carbon[2, ibnd]
    τ_ssag = τ_ssa * organic_carbon[3, ibnd]
    return τ, τ_ssa, τ_ssag
end

function compute_lookup_organic_carbon_props(lkp_aero, ibnd::Int, aeromass::FT, rh::FT) where {FT}
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

function locate_merra_size_bin(size_bin_limits, aerosize)
    nbins = size(size_bin_limits, 2)
    bin = 1
    for ibin in 1:nbins
        if size_bin_limits[1, ibin] ≤ aerosize ≤ size_bin_limits[2, ibin]
            bin = ibin
            break
        end
    end
    return bin
end
