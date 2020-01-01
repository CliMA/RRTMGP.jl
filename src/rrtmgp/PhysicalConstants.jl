"""
    PhysicalConstants

Physical and mathematical constants used in gas optics calculation
"""
module PhysicalConstants

export k_boltz,
       m_h2o,
       avogad,
       R_univ_gconst,
       m_dry,
       grav,
       cp_dry

"""
    k_boltz(::Type{FT})

Physical constants, 2018 SI definition of metric system

Boltzmann constant [J/K] = [(kg m^2)/(K s^2)]

Refs:

 - `doi:10.1088/1681-7575/aa950a`
 - `https://www.nist.gov/si-redefinition/meet-constants`

"""
k_boltz(::Type{FT}) where FT = FT(1.380649e-23)

"""
    m_h2o(::Type{FT})

molecular weight of water [kg/mol]
"""
m_h2o(::Type{FT}) where FT = FT(0.018016)

"""
    avogad(::Type{FT})

Avogadro's number [molec/mol]
"""
avogad(::Type{FT}) where FT = FT(6.02214076e23)

"""
    R_univ_gconst(::Type{FT})

Universal gas constant [J/(mol K)]
"""
R_univ_gconst(::Type{FT}) where FT = avogad(FT) * k_boltz(FT)

#### Constants specific to the earth's
# atmosphere -- changeable because they
# might be different on, e.g., other planets

"""
    m_dry(::Type{FT})

molecular weight of dry air [kg/mol]
"""
m_dry(::Type{FT}) where FT = FT(0.028964)

"""
    grav(::Type{FT})

Gravity at Earth's surface [m/s2]
"""
grav(::Type{FT}) where FT = FT(9.80665)

"""
    cp_dry(::Type{FT})

Specific heat at constant pressure for dry air [J/(K kg)]
"""
cp_dry(::Type{FT}) where FT = FT(1004.64)

end
