using RRTMGP.Gases
using RRTMGP.GasConcentrations
using RRTMGP.GasOptics
using RRTMGP.ReferenceStates
using NCDatasets

function read_char_vec(ds, var_name)
    tmp = ds[var_name][:]
    return String[strip(join(tmp[:, i])) for i in 1:size(tmp, 2)]
end

read_gases(ds, var_name) = convert.(AbstractGas, read_char_vec(ds, var_name))

"""
    load_and_init(ds, ::Type{FT}, gases_prescribed::Vector{AbstractGas}) where {FT<:AbstractFloat}

Initialize the gas optics class with data. The calls look slightly different depending
  on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
"""
function load_and_init(ds, ::Type{FT}, gases_prescribed::Vector{AbstractGas}) where {FT <: AbstractFloat}

    # Reading the properties from the NetCDF file

    I = Int
    kmajor = Array{FT}(ds["kmajor"][:])
    key_species = Array{I}(ds["key_species"][:])
    band2gpt = Array{I}(ds["bnd_limits_gpt"][:])
    band_lims_wavenum = Array{FT}(ds["bnd_limits_wavenumber"][:])
    gases_in_database = read_gases(ds, "gas_names")
    gas_minor = read_gases(ds, "gas_minor")
    identifier_minor = read_gases(ds, "identifier_minor")

    rayl_lower = haskey(ds, "rayl_lower") ? Array{FT}(ds["rayl_lower"][:]) : nothing
    rayl_upper = haskey(ds, "rayl_upper") ? Array{FT}(ds["rayl_upper"][:]) : nothing
    @assert haskey(ds, "rayl_lower") == haskey(ds, "rayl_upper")

    lower = GasOpticsVars{FT, I}(
        Array{I}(ds["minor_limits_gpt_lower"][:]),
        Array{Bool}(ds["minor_scales_with_density_lower"][:]),
        Array{Bool}(ds["scale_by_complement_lower"][:]),
        Array{I}(ds["kminor_start_lower"][:]),
        Array{FT}(ds["kminor_lower"][:]),
        read_gases(ds, "scaling_gas_lower"),
        read_gases(ds, "minor_gases_lower"),
    )
    upper = GasOpticsVars{FT, I}(
        Array{I}(ds["minor_limits_gpt_upper"][:]),
        Array{Bool}(ds["minor_scales_with_density_upper"][:]),
        Array{Bool}(ds["scale_by_complement_upper"][:]),
        Array{I}(ds["kminor_start_upper"][:]),
        Array{FT}(ds["kminor_upper"][:]),
        read_gases(ds, "scaling_gas_upper"),
        read_gases(ds, "minor_gases_upper"),
    )

    ref = ReferenceState(
        Array{FT}(ds["press_ref"][:]),
        Array{FT}(ds["temp_ref"][:]),
        FT(ds["press_ref_trop"][:]),
        Array{FT}(ds["vmr_ref"][:]),
        gases_prescribed,
        gases_in_database,
    )

    optical_props = OpticalPropsBase("GasOptics optical props", band_lims_wavenum, band2gpt)

    args = (
        gases_prescribed,
        gases_in_database,
        key_species,
        optical_props,
        kmajor,
        gas_minor,
        identifier_minor,
        lower,
        upper,
    )

    if haskey(ds, "totplnk")
        totplnk = Array{FT}(ds["totplnk"][:])
        planck_frac = Array{FT}(ds["plank_fraction"][:])
        kdist = get_k_dist_lw(totplnk, planck_frac, rayl_lower, rayl_upper, ref, args...)
    else
        solar_src = ds["solar_source"][:]
        kdist = get_k_dist_sw(solar_src, rayl_lower, rayl_upper, ref, args...)
    end
    return kdist

end
