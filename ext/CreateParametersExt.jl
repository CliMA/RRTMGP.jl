module CreateParametersExt

import RRTMGP.Parameters.RRTMGPParameters
import CLIMAParameters as CP

RRTMGPParameters(::Type{FT}, overrides = NamedTuple()) where {FT <: AbstractFloat} =
    RRTMGPParameters(CP.create_toml_dict(FT), overrides)

function RRTMGPParameters(toml_dict::CP.AbstractTOMLDict, overrides = NamedTuple())
    name_map = (;
        :gravitational_acceleration => :grav,
        :molar_mass_dry_air => :molmass_dryair,
        :molar_mass_water => :molmass_water,
        :gas_constant => :gas_constant,
        :adiabatic_exponent_dry_air => :kappa_d,
        :stefan_boltzmann_constant => :Stefan,
        :avogadro_constant => :avogad,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "RRTMGP")
    parameters = merge(parameters, overrides)
    FT = CP.float_type(toml_dict)
    return RRTMGPParameters{FT}(; parameters...)
end

end
