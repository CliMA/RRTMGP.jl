module Parameters

abstract type AbstractRRTMGPParameters end
const ARP = AbstractRRTMGPParameters

Base.@kwdef struct RRTMGPParameters{FT} <: ARP
    grav::FT
    molmass_dryair::FT
    molmass_water::FT
    gas_constant::FT
    kappa_d::FT
    Stefan::FT
    avogad::FT
    # valid temperature range of the optics lookup tables [K]; `clip!` clamps
    # the spectral solvers' temperature inputs to it
    optics_lookup_temperature_min::FT = 160
    optics_lookup_temperature_max::FT = 355
end

# Method wrappers
for var in fieldnames(RRTMGPParameters)
    @eval $var(ps::ARP) = ps.$var
end

# Derived parameters
R_d(ps::ARP) = gas_constant(ps) / molmass_dryair(ps)
cp_d(ps::ARP) = R_d(ps) / kappa_d(ps)

end # module
