import CLIMAParameters as CP
import RRTMGP.Parameters as RP

#=
include(joinpath(pkgdir(RRTMGP), "parameters", "create_parameters.jl"))
param_set = create_insolation_parameters(FT)
=#
function create_insolation_parameters(FT, overrides::NamedTuple = NamedTuple())
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    aliases = string.(fieldnames(RP.RRTMGPParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "RRTMGP")
    params = merge((; pairs...), overrides) # overrides
    param_set = RP.RRTMGPParameters{FT}(; params...)
    return param_set
end
