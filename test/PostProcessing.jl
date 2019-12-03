#### PostProcessing

using RRTMGP.PhysicalConstants

"""
    compute_heating_rate(flux_up::Array{FT},
                         flux_dn::Array{FT},
                         as::AtmosphericState{FT}) where {FT<:AbstractFloat}

Compute heating rate from fluxes

Heating rate H [K/sec] = 1/(ρ c_p) d f_net/d z
Here use hydrostatic equation for density and heat capacity of dry air
"""
function compute_heating_rate(flux_up::Array{FT},
                              flux_dn::Array{FT},
                              as::AtmosphericState{FT}) where {FT<:AbstractFloat}

  ncol = as.ncol
  nlay = as.nlay
  heating_rate = Array{FT}(undef, ncol, nlay)
  z = convert(Array,as.p_lay[1,:])
  for ilay in 1:nlay
    heating_rate[:,ilay] .= (flux_up[:,ilay+1] .- flux_up[:,ilay] .-
                             flux_dn[:,ilay+1] .+ flux_dn[:,ilay]) .*
                             grav(FT) ./ (cp_dry(FT) .* (as.p_lev[:,ilay+1] .- as.p_lev[:,ilay]))
  end
  heating_rate = convert(Array,sum(heating_rate,dims=1)'/ncol)
  return heating_rate, z
end