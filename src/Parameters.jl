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
end

# Method wrappers
for var in fieldnames(RRTMGPParameters)
    @eval $var(ps::ARP) = ps.$var
end

# Derived parameters
R_d(ps::ARP) = gas_constant(ps) / molmass_dryair(ps)
cp_d(ps::ARP) = R_d(ps) / kappa_d(ps)

Base.eltype(::RRTMGPParameters{FT}) where {FT} = FT

end # module
