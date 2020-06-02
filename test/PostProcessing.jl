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
