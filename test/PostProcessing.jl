#### PostProcessing
using CLIMAParameters
using CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
using CLIMAParameters.Planet: cp_d
const param_set = EarthParameterSet()

"""
    compute_heating_rate(flux_up::Array{FT},
                         flux_dn::Array{FT},
                         as::AtmosphericState{FT}) where {FT<:AbstractFloat}

Compute heating rate from fluxes

Heating rate H [K/sec] = 1/(ρ c_p) d f_net/d z
Here use hydrostatic equation for density and heat capacity of dry air
"""
function compute_heating_rate(
    flux_up::Array{FT},
    flux_dn::Array{FT},
    as::AtmosphericState{FT},
) where {FT<:AbstractFloat}

    ncol = as.ncol
    nlay = as.nlay
    heating_rate = Array{FT}(undef, ncol, nlay)
    z = convert(Array, as.p_lay[1, :])
    for ilay = 1:nlay
        heating_rate[:, ilay] .=
            (
                flux_up[:, ilay+1] .- flux_up[:, ilay] .- flux_dn[:, ilay+1] .+
                flux_dn[:, ilay]
            ) .* grav(FT) ./
            (FT(cp_d(param_set)) .* (as.p_lev[:, ilay+1] .- as.p_lev[:, ilay]))
    end
    heating_rate = convert(Array, sum(heating_rate, dims = 1)' / ncol)
    return heating_rate, z
end

function rel_err(flux_up, flux_dn, ref_up, ref_dn)
    @assert size(flux_up) == size(ref_up)
    @assert size(flux_dn) == size(ref_dn)
    @assert size(ref_up) == size(ref_dn)

    FT = eltype(flux_up)
    tol_FT, rel_err_up, rel_err_dn = sqrt(eps(FT)), FT(0), FT(0)

    for k = 1:size(ref_up, 3), j = 1:size(ref_up, 2), i = 1:size(ref_up, 1)
        diff_up = abs(flux_up[i, j, k] - ref_up[i, j, k])
        diff_dn = abs(flux_dn[i, j, k] - ref_dn[i, j, k])

        if abs(ref_up[i, j, k]) > tol_FT
            rel_err_up = max(rel_err_up, diff_up / ref_up[i, j, k])
        else
            rel_err_up = max(rel_err_up, diff_up)
        end

        if abs(ref_dn[i, j, k]) > tol_FT
            rel_err_dn = max(rel_err_dn, diff_dn / ref_dn[i, j, k])
        else
            rel_err_dn = max(rel_err_dn, diff_dn)
        end
    end
    return rel_err_up, rel_err_dn
end
