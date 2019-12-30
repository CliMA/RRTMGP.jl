"""
    GasOptics

Computes spectrally-resolved gas optical properties and source functions
 given atmospheric physical properties (profiles of temperature, pressure, and gas concentrations)
 The `KDistributionLongwave` and `KDistributionShortwave` must be initialized with data (provided
 as a netCDF file) before being used.

Two variants apply to internal Planck sources (longwave radiation in the Earth's atmosphere) and to
 external stellar radiation (shortwave radiation in the Earth's atmosphere).
 The variant is chosen based on what information is supplied during initialization.
 (It might make more sense to define two sub-classes)
"""
module GasOptics

using TimerOutputs
using DocStringExtensions
using StaticArrays

const to_gor = TimerOutput()

using ..Gases
using ..ReferenceStates
using ..FortranIntrinsics
using ..PhysicalConstants
using ..Utilities
using ..OpticalProps
using ..SourceFunctions
using ..GasConcentrations
using ..OpticalProps
using ..AtmosphericStates

export gas_optics!, get_k_dist_lw, get_k_dist_sw
export source_is_internal, source_is_external
export GasOpticsVars

"""
    GasOpticsAux{I}

Auxiliary indexes for variables in the upper
and lower levels of the atmosphere.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GasOpticsAux{I<:Int}
  "indexes for determining `col_gas`"
  idx_minor::Vector{I}
  "indexes that have special treatment in density scaling"
  idx_minor_scaling::Vector{I}
  function GasOpticsAux(gas_names_present, gas_minor, identifier_minor, reduced_atmos, ::Type{I}) where I
    return new{I}(create_idx_minor(gas_names_present,
                                   gas_minor,
                                   identifier_minor,
                                   reduced_atmos.minor_gases),
                  create_idx_minor_scaling(gas_names_present,
                                           reduced_atmos.scaling_gas))
  end
end

"""
    GasOpticsVars{FT,I}

Variables defined at

  - upper  (log(p) < ref.press_trop_log)
  - lower  (log(p) > ref.press_trop_log)

levels of the atmosphere for both full and reduced sets of gases.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GasOpticsVars{FT<:AbstractFloat,I<:Int}
  "minor g-point limits"
  minor_limits_gpt::Array{I,2}
  "minor scales with density"
  minor_scales_with_density::Vector{Bool}
  "scale by complement"
  scale_by_complement::Vector{Bool}
  "kminor start"
  kminor_start::Vector{I}
  "absorption coefficient of minor species (n_minor,η,temperature)"
  kminor::Array{FT,3}
  "scaling gas"
  scaling_gas::Array{AbstractGas}
  "minor gases"
  minor_gases::Array{AbstractGas}
end

abstract type AbstractGasOptics{FT,I} <: AbstractOpticalProps{FT,I} end

"""
    KDistributionLongwave{FT,I} <: AbstractGasOptics{FT,I}

Gas optics with internal sources (for longwave radiation)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct KDistributionLongwave{FT,I} <: AbstractGasOptics{FT,I}
  "Reference state data"
  ref::ReferenceState{FT}
  "Base optical properties"
  optical_props::OpticalPropsBase{FT,I}
  "GasOpticsVars in the lower atmosphere"
  lower::GasOpticsVars{FT,I}
  "GasOpticsVars in the upper atmosphere"
  upper::GasOpticsVars{FT,I}
  "Auxiliary variables (index maps) in the lower atmosphere"
  lower_aux::GasOpticsAux{I}
  "Auxiliary variables (index maps) in the upper atmosphere"
  upper_aux::GasOpticsAux{I}
  "Present gas names"
  gas_names::Vector{AbstractGas}
  "Indexes of present gas names"
  i_gases::Vector{I}
  "Absorption coefficient for major species (g-point,η,pressure,temperature)"
  kmajor::Array{FT,4}
  "major species pair; [2, nflav]"
  flavor::Array{I,2}
  "Flavor per g-point: `lower.flavor = gpoint_flavor[1, 1:ngpt]`, `upper.flavor = gpoint_flavor[2, 1:ngpt]`"
  gpoint_flavor::Array{I,2}
  "Indicates whether a key species is in any band"
  is_key::Vector{Bool}
  "Stored fraction of Planck irradiance in band for given g-point"
  planck_frac::Array{FT,4}
  "integrated Planck irradiance by band; [Planck temperatures,band]"
  totplnk::Array{FT,2}
  "temperature steps in totplnk"
  totplnk_delta::FT
  "Absorption coefficient for Rayleigh scattering [g-point,η,temperature,upper/lower atmosphere]"
  krayl::Union{Array{FT,4},Nothing}
end

"""
    KDistributionShortwave{FT,I} <: AbstractGasOptics{FT,I}

Gas optics with external sources (for shortwave radiation)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct KDistributionShortwave{FT,I} <: AbstractGasOptics{FT,I}
  "Reference state data"
  ref::ReferenceState{FT}
  "Base optical properties"
  optical_props::OpticalPropsBase{FT,I}
  "GasOpticsVars in the lower atmosphere"
  lower::GasOpticsVars{FT,I}
  "GasOpticsVars in the upper atmosphere"
  upper::GasOpticsVars{FT,I}
  "Auxiliary variables (index maps) in the lower atmosphere"
  lower_aux::GasOpticsAux{I}
  "Auxiliary variables (index maps) in the upper atmosphere"
  upper_aux::GasOpticsAux{I}
  "Present gas names"
  gas_names::Vector{AbstractGas}
  "Indexes of present gas names"
  i_gases::Vector{I}
  "Absorption coefficient for major species (g-point,η,pressure,temperature)"
  kmajor::Array{FT,4}
  "major species pair; [2, nflav]"
  flavor::Array{I,2}
  "Flavor per g-point: lower.flavor = gpoint_flavor[1, g-point], upper.flavor = gpoint_flavor[2, g-point]"
  gpoint_flavor::Array{I,2}
  "Indicates whether a key species is in any band"
  is_key::Vector{Bool}
  "incoming solar irradiance (g-point)"
  solar_src::Vector{FT} # incoming solar irradiance(g-point)
  "absorption coefficient for Rayleigh scattering (g-point,η,temperature,upper/lower atmosphere)"
  krayl::Union{Array{FT,4},Nothing}
end

"""
    InterpolationCoefficients{FT,I}

Interpolation coefficients

# Fields
$(DocStringExtensions.FIELDS)
"""
struct InterpolationCoefficients{FT<:AbstractFloat,I<:Int}
  "index for temperature"
  jtemp::Array{I,2}
  "index for pressure"
  jpress::Array{I,2}
  "index for binary species parameter (η)"
  j_η::Array{I,4}
  "fractions for major species"
  fmajor::Array{FT,6}
  "fractions for minor species. [reference η level (temperature dependent), reference temperature level, flavor, layer]"
  fminor::Array{FT,5}
  "combination of major specie's column amounts"
  col_mix::Array{FT,4}
end
function InterpolationCoefficients(::Type{FT}, ::Type{I}, ncol, nlay, nflav) where {I<:Int, FT<:AbstractFloat}
  jtemp = Array{I}(undef,                ncol, nlay)
  jpress = Array{I}(undef,               ncol, nlay)
  j_η = Array{I}(undef, 2, nflav,        ncol, nlay)
  fmajor = zeros(FT, 2,2,2, nflav,       ncol, nlay)
  fminor  = Array{FT}(undef, 2,2, nflav, ncol, nlay)
  col_mix  = Array{FT}(undef, 2, nflav,  ncol, nlay)
  return InterpolationCoefficients{FT,I}(jtemp,jpress,j_η,fmajor,fminor,col_mix)
end

"""
    InterpolationCoefficientsPGP{FT,I}

Interpolation coefficients per grid point

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct InterpolationCoefficientsPGP{FT<:AbstractFloat,I<:Int}
  "index for temperature"
  jtemp::I
  "index for pressure"
  jpress::I
  "index for binary species parameter (η)"
  j_η::Array{I,2}
  "fractions for major species"
  fmajor::Array{FT,4}
  "fractions for minor species. [reference η level (temperature dependent), reference temperature level, flavor, layer]"
  fminor::Array{FT,3}
  "combination of major specie's column amounts"
  col_mix::Array{FT,2}
end

function Base.convert(::Type{InterpolationCoefficients}, data::Array{InterpolationCoefficientsPGP{FT,I}}) where {FT,I}
  s = size(data)
  s_j_η = size(first(data).j_η)
  s_fmajor = size(first(data).fmajor)
  s_fminor = size(first(data).fminor)
  s_col_mix = size(first(data).col_mix)
  jtemp  = Array{I}([data[i,j].jtemp  for i in 1:s[1], j in 1:s[2]])
  jpress = Array{I}([data[i,j].jpress for i in 1:s[1], j in 1:s[2]])

  j_η    = Array{I}([data[i,j].j_η[k,p] for k in 1:s_j_η[1],
                                             p in 1:s_j_η[2],
                                             i in 1:s[1],
                                             j in 1:s[2]])
  fmajor = Array{FT}([data[i,j].fmajor[k,p,q,r] for k in 1:s_fmajor[1],
                                                    p in 1:s_fmajor[2],
                                                    q in 1:s_fmajor[3],
                                                    r in 1:s_fmajor[4],
                                                    i in 1:s[1],
                                                    j in 1:s[2]])
  fminor = Array{FT}([data[i,j].fminor[k,p,q] for k in 1:s_fminor[1],
                                                  p in 1:s_fminor[2],
                                                  q in 1:s_fminor[3],
                                                  i in 1:s[1],
                                                  j in 1:s[2]])
  col_mix = Array{FT}([data[i,j].col_mix[k,p] for k in 1:s_col_mix[1],
                                                  p in 1:s_col_mix[2],
                                                  i in 1:s[1],
                                                  j in 1:s[2]])
  return InterpolationCoefficients{FT,I}(jtemp,jpress,j_η,fmajor,fminor,col_mix)
end

Base.convert(::Type{Array{InterpolationCoefficientsPGP}}, data::InterpolationCoefficients{FT,I}) where {FT,I} =
  [InterpolationCoefficientsPGP{FT,I}(data.jtemp[i,j],
                                      data.jpress[i,j],
                                      data.j_η[:,:,i,j],
                                      data.fmajor[:,:,:,:,i,j],
                                      data.fminor[:,:,:,i,j],
                                      data.col_mix[:,:,i,j]) for i in 1:size(data.jtemp,1), j in 1:size(data.jtemp,2)]

"""
    gas_optics!(this::KDistributionLongwave{FT,I},
                as::AtmosphericState{FT,I},
                optical_props::AbstractOpticalPropsArry{FT,I},
                sources::SourceFuncLongWave{FT,I}) where {FT<:AbstractFloat,I<:Int}

Compute gas optical depth and Planck source functions given:

 - `this` gas optics, see [`KDistributionLongwave`](@ref)
 - `as` atmospheric state, see [`AtmosphericState`](@ref)
 - `optical_props` optical properties, see [`AbstractOpticalPropsArry`](@ref)
 - `sources` longwave sources, see [`SourceFuncLongWave`](@ref)
"""
function gas_optics!(this::KDistributionLongwave{FT,I},
                     as::AtmosphericState{FT,I},
                     optical_props::AbstractOpticalPropsArry{FT,I},
                     sources::SourceFuncLongWave{FT,I}) where {FT<:AbstractFloat,I<:Int}
  ics = InterpolationCoefficients(FT, I, as.ncol, as.nlay, get_nflav(this))
  ics = convert(Array{InterpolationCoefficientsPGP}, ics)
  ics = convert(InterpolationCoefficients, ics)
  compute_gas_τs!(ics, this, as, optical_props)
  @assert get_ncol(sources) == as.ncol
  @assert get_nlay(sources) == as.nlay
  @assert get_ngpt(sources) == get_ngpt(this.optical_props)
  source!(sources, as, ics, this)
  return nothing
end

"""
    gas_optics!(this::KDistributionShortwave{FT,I},
                as::AtmosphericState{FT,I},
                optical_props::AbstractOpticalPropsArry{FT,I},
                last_call=false) where {FT<:AbstractFloat,I<:Int}

Compute gas optical depth given:

 - `this` gas optics, see [`KDistributionShortwave`](@ref)
 - `as` atmospheric state, see [`AtmosphericState`](@ref)
 - `optical_props` optical properties, see [`AbstractOpticalPropsArry`](@ref)
"""
function gas_optics!(this::KDistributionShortwave{FT,I},
                     as::AtmosphericState{FT,I},
                     optical_props::AbstractOpticalPropsArry{FT,I},
                     last_call=false) where {FT<:AbstractFloat,I<:Int}
  ics = InterpolationCoefficients(FT, I, as.ncol,as.nlay, get_nflav(this))
  ics = convert(Array{InterpolationCoefficientsPGP}, ics)
  ics = convert(InterpolationCoefficients, ics)
  @timeit to_gor "compute_gas_τs!" compute_gas_τs!(ics, this, as, optical_props, last_call)
  last_call && @show to_gor
  return nothing
end


"""
    compute_gas_τs!(ics::InterpolationCoefficients{FT,I},
                    this::AbstractGasOptics{FT},
                    as::AtmosphericState{FT,I},
                    optical_props::AbstractOpticalPropsArry{FT,I},
                    last_call=false) where {FT<:AbstractFloat,I<:Int}

Returns optical properties and interpolation coefficients

 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)
 - `this` gas optics, see [`KDistributionLongwave`](@ref) or [`KDistributionShortwave`](@ref)
 - `as` atmospheric state, see [`AtmosphericState`](@ref)
 - `optical_props` optical properties, see [`AbstractOpticalPropsArry`](@ref)

Local variables
 - `τ`, `τ_Rayleigh` [ngpt,nlay,ncol] absorption, Rayleigh scattering optical depths
 - `col_mix` combination of major specie's column amounts [reference temperature level, flavor, column]
"""
function compute_gas_τs!(ics::InterpolationCoefficients{FT,I},
                         this::AbstractGasOptics{FT},
                         as::AtmosphericState{FT,I},
                         optical_props::AbstractOpticalPropsArry{FT,I},
                         last_call=false) where {FT<:AbstractFloat,I<:Int}
  ncol  = get_ncol(optical_props)
  nlay  = get_nlay(optical_props)
  ngpt  = get_ngpt(optical_props)
  τ          = Array{FT}(undef, ngpt,nlay,ncol) # absorption, Rayleigh scattering optical depths
  τ_Rayleigh = Array{FT}(undef, ngpt,nlay,ncol) # absorption, Rayleigh scattering optical depths
  # Check for presence of key species in GasConcs; return error if any key species are not present
  check_key_species_present(this, as.gas_conc.gas_names)
  @timeit to_gor "interpolation!" interpolation!(ics, this, as)
  @timeit to_gor "compute_τ_absorption!" compute_τ_absorption!(τ, this, ics, as, last_call)
  if this.krayl ≠ nothing
    @timeit to_gor "compute_τ_Rayleigh!" compute_τ_Rayleigh!(τ_Rayleigh, this, ics, as)
  end
  # Combine optical depths and reorder for radiative transfer solver.
  @timeit to_gor "combine_and_reorder!" combine_and_reorder!(optical_props, τ, τ_Rayleigh, this.krayl ≠ nothing)
  return nothing
end

"""
    source!(args...)

Compute source functions, currently set to Planck source
"""
source!(args...) = compute_Planck_source!(args...)

#####
##### Initialization
#####

"""
    get_k_dist_lw(totplnk,
                  planck_frac,
                  rayl_lower,
                  rayl_upper,
                  ref::ReferenceState{FT},
                  args...) where {FT<:AbstractFloat}

Initialize object based on data read from netCDF file however the user desires.
Rayleigh scattering tables may or may not be present; this is indicated with allocation status
This interface is for the internal-sources object -- includes Planck functions and fractions
"""
function get_k_dist_lw(totplnk,
                       planck_frac,
                       rayl_lower,
                       rayl_upper,
                       ref::ReferenceState{FT},
                       args...) where {FT<:AbstractFloat}
  abs_coeffs = init_abs_coeffs(args...)

  if rayl_lower ≠ nothing && rayl_upper ≠ nothing
    krayl = zeros(FT, size(rayl_lower)...,2)
    krayl[:,:,:,1] = rayl_lower
    krayl[:,:,:,2] = rayl_upper
  else
    krayl = nothing
  end
  totplnk_delta =  (ref.temp_max-ref.temp_min) / (size(totplnk, 1)-1)

  return KDistributionLongwave{FT,Int}(ref,
                                       abs_coeffs...,
                                       planck_frac,
                                       totplnk,
                                       totplnk_delta,
                                       krayl)

end

"""
    get_k_dist_sw(solar_src,
                  rayl_lower,
                  rayl_upper,
                  ref::ReferenceState{FT},
                  args...) where {FT<:AbstractFloat}

Initialize object based on data read from netCDF file however the user desires.
Rayleigh scattering tables may or may not be present; this is indicated with allocation status
This interface is for the external-sources object -- includes TOA source function table
"""
function get_k_dist_sw(solar_src,
                       rayl_lower,
                       rayl_upper,
                       ref::ReferenceState{FT},
                       args...) where {FT<:AbstractFloat}
  abs_coeffs = init_abs_coeffs(args...)

  if rayl_lower ≠ nothing && rayl_upper ≠ nothing
    krayl = zeros(FT, size(rayl_lower)...,2)
    krayl[:,:,:,1] = rayl_lower
    krayl[:,:,:,2] = rayl_upper
  else
    krayl = nothing
  end

  return KDistributionShortwave{FT,Int}(ref,
                                        abs_coeffs...,
                                        solar_src,
                                        krayl)

end

"""
    init_abs_coeffs(gases_prescribed::Vector{AbstractGas},
                    gases_in_database::Vector{AbstractGas},
                    key_species::Array{I,3},
                    optical_props::OpticalPropsBase{FT,I},
                    kmajor::Array{FT,4},
                    gas_minor::Vector{AbstractGas},
                    identifier_minor::Vector{AbstractGas},
                    lower::GasOpticsVars{FT},
                    upper::GasOpticsVars{FT}) where {FT<:AbstractFloat,I<:Int}

Initialize absorption coefficient arrays

 - `lower` lower atmospheric gas optics variables, see [`GasOpticsVars`](@ref)
 - `upper` upper atmospheric gas optics variables, see [`GasOpticsVars`](@ref)
 - `lower_aux` lower atmospheric gas optics auxiliary variables, see [`GasOpticsAux`](@ref)
 - `upper_aux` upper atmospheric gas optics auxiliary variables, see [`GasOpticsAux`](@ref)
 - `gases_in_database` gases available in database
 - `optical_props` optical properties, see [`OpticalProps`](@ref)
"""
function init_abs_coeffs(gases_prescribed::Vector{AbstractGas},
                         gases_in_database::Vector{AbstractGas},
                         key_species::Array{I,3},
                         optical_props::OpticalPropsBase{FT,I},
                         kmajor::Array{FT,4},
                         gas_minor::Vector{AbstractGas},
                         identifier_minor::Vector{AbstractGas},
                         lower::GasOpticsVars{FT},
                         upper::GasOpticsVars{FT}) where {FT<:AbstractFloat,I<:Int}

  # Which gases known to the gas optics are present in the host model (gases_prescribed)?
  gas_names_present = intersect(gases_in_database, gases_prescribed)
  i_gases = I[loc_in_array(gas, gas_names_present) for gas in gas_names_present]

  reduced_lower = reduce_minor_arrays(gases_prescribed, gas_minor, identifier_minor, lower)
  reduced_upper = reduce_minor_arrays(gases_prescribed, gas_minor, identifier_minor, upper)

  lower_aux = GasOpticsAux(gas_names_present, gas_minor, identifier_minor, reduced_lower, I)
  upper_aux = GasOpticsAux(gas_names_present, gas_minor, identifier_minor, reduced_upper, I)

  # create flavor list
  # Reduce (remap) key_species list; checks that all key gases are present in incoming
  key_species_red = create_key_species_reduce(gases_in_database, gas_names_present, key_species)

  flavor = create_flavor(key_species_red)

  # create gpoint_flavor list
  gpoint_flavor = create_gpoint_flavor(key_species_red, get_gpoint_bands(optical_props), flavor)

  # Which species are key in one or more bands?
  #   flavor is an index into gas_names_present
  is_key = Bool[false for i in 1:length(gas_names_present)]
  for j in 1:size(flavor, 2)
    for i in 1:size(flavor, 1) # should be 2
      if flavor[i,j] ≠ 0
        is_key[flavor[i,j]] = true
      end
    end
  end

  return (optical_props,
          reduced_lower,
          reduced_upper,
          lower_aux,
          upper_aux,
          gas_names_present,
          i_gases,
          kmajor,
          flavor,
          gpoint_flavor,
          is_key)

end

"""
    check_key_species_present(this::AbstractGasOptics, gas_names::Vector{AbstractGas})

Ensure that every key gas required by the k-distribution is
present in the gas concentration object
"""
function check_key_species_present(this::AbstractGasOptics, gas_names::Vector{AbstractGas})
  key_gas_names = pack(this.gas_names, this.is_key)
  for igas = 1:length(key_gas_names)
    @assert key_gas_names[igas] in gas_names
  end
end

"""
    get_minor_list(this::AbstractGasOptics, gas_names::Vector{AbstractGas}, names_spec::Vector{AbstractGas})

List of minor gases to be used in gas_optics()
Function to define names of key and minor gases to be used by gas_optics().
The final list gases includes those that are defined in gas_optics_specification
and are provided in GasConcs.

!!! Not yet tested
"""
get_minor_list(this::AbstractGasOptics, gas_names::Vector{AbstractGas}, names_spec::Vector{AbstractGas}) =
  pack(this.gas_names, map(x->x in gas_names, names_spec))

#####
##### Inquiry functions
#####

"""
    source_is_internal(this::AbstractGasOptics)

Bool indicating if initialized for internal sources
"""
source_is_internal(::KDistributionShortwave) = false
source_is_internal(::KDistributionLongwave) = true

"""
    source_is_external(this::AbstractGasOptics)

Bool indicating if initialized for external sources
"""
source_is_external(::KDistributionShortwave) = true
source_is_external(::KDistributionLongwave) = false

"""
    rewrite_key_species_pair(key_species_pair)

(0,0) becomes (2,2) -- because absorption coefficients for these g-points will be 0.
"""
rewrite_key_species_pair(key_species_pair) =
  all(key_species_pair .== [0,0]) ? [2,2] : key_species_pair

"""
    key_species_pair_exists(key_species_list, key_species_pair)

True is key_species_pair exists in key_species_list

integer, dimension(:,:), intent(in) :: key_species_list
integer, dimension(2),   intent(in) :: key_species_pair
"""
key_species_pair_exists(key_species_list, key_species_pair) =
  any([all(key_species_list[:,i] .== key_species_pair) for i=1:size(key_species_list,2)])


"""
    create_flavor(key_species)

Flavor list -- an unordered array of extent (2,:)
containing all possible pairs of key species
used in either upper or lower atmosphere
"""
function create_flavor(key_species::Array{I,3}) where {I<:Int}

  key_species_list = Array{Int}(undef, 2,size(key_species,3)*2)

  # prepare list of key_species
  i = 1
  for ibnd=1:size(key_species,3)
    for iatm=1:size(key_species,1)
      key_species_list[:,i] .= key_species[:,iatm,ibnd]
      i += 1
    end
  end
  # rewrite single key_species pairs
  for i=1:size(key_species_list,2)
    key_species_list[:,i] .= rewrite_key_species_pair(key_species_list[:,i])
  end
  # count unique key species pairs
  x = key_species_list
  key_species_pair_absent = [!key_species_pair_exists(x[:,1:i-1],x[:,i]) for i=1:size(x,2)]
  iflavor = count(key_species_pair_absent)
  # fill flavors
  flavor = Array{Int}(undef, 2, iflavor)
  iflavor = 0
  for i=1:size(key_species_list,2)
    if key_species_pair_absent[i]
      iflavor += 1
      flavor[:,iflavor] = key_species_list[:,i]
    end
  end
  return flavor
end

"""
    create_idx_minor(gas_names::VG,
                     gas_minor::VG,
                     identifier_minor::VG,
                     minor_gases_atm::VG) where {VG<:Vector{AbstractGas}}

create index list for extracting col_gas needed for minor gas optical depth calculations

 - `gas_names`
 - `gas_minor`
 - `identifier_minor`
 - `minor_gases_atm`
"""
function create_idx_minor(gas_names::VG,
                          gas_minor::VG,
                          identifier_minor::VG,
                          minor_gases_atm::VG) where {VG<:Vector{AbstractGas}}
  # Find identifying gas for minor species in list of possible identifiers (e.g. h2o_slf)
  idx_mnr = map(x->loc_in_array(x, identifier_minor), minor_gases_atm)
  # Find name of gas associated with minor species identifier (e.g. h2o)
  return map(x->loc_in_array(gas_minor[x], gas_names), idx_mnr)
end

"""
    create_idx_minor_scaling(gas_names::VG, scaling_gas_atm::VG) where {VG<:Vector{AbstractGas}}

Create index for special treatment in density scaling of minor gases

 - `gas_names` gas names
 - `scaling_gas_atm`
"""
create_idx_minor_scaling(gas_names::VG, scaling_gas_atm::VG) where {VG<:Vector{AbstractGas}} =
  map(x->loc_in_array(x, gas_names), scaling_gas_atm)

"""
    create_key_species_reduce(gas_names::VG,
                              gas_names_red::VG,
                              key_species::Array{I,3}) where {VS<:Vector{AbstractGas},I<:Int}

create flavor list
Reduce (remap) key_species list; checks that all key gases are present in incoming
"""
function create_key_species_reduce(gas_names::VG,
                                   gas_names_red::VG,
                                   key_species::Array{I,3}) where {VG<:Vector{AbstractGas},I<:Int}
  key_species_red = map(x-> x≠0 ? loc_in_array(gas_names[x],gas_names_red) : x, key_species)
  @assert !any(key_species_red .== -1)
  return key_species_red
end

"""
    reduce_minor_arrays(gases_prescribed::Vector{AbstractGas},
                        gas_minor::Vector{AbstractGas},
                        identifier_minor::Vector{AbstractGas},
                        atmos::GasOpticsVars{FT}) where FT

A reduced `GasOpticsVars` for minor species arrays.
Variables only contain minor gases that are available.

 - `gases_prescribed` array of available gases
 - `gas_minor` array of minor gases
 - `identifier_minor`
 - `atmos` original gas optics for minor species, see [`GasOpticsVars`](@ref)
"""
function reduce_minor_arrays(gases_prescribed::Vector{AbstractGas},
                             gas_minor::Vector{AbstractGas},
                             identifier_minor::Vector{AbstractGas},
                             atmos::GasOpticsVars{FT}) where FT

  I = Int
  nm = length(atmos.minor_gases)
  tot_g=0

  mask = map(x->loc_in_array(x, identifier_minor), atmos.minor_gases)
  gas_is_present = map(x->x in gases_prescribed, gas_minor[mask])

  for i = 1:length(atmos.minor_gases)
    if gas_is_present[i]
      tot_g = tot_g + (atmos.minor_limits_gpt[2,i]-atmos.minor_limits_gpt[1,i]+1)
    end
  end
  red_nm = count(gas_is_present)

  if red_nm == nm
    atmos_red_minor_gases               = atmos.minor_gases
    atmos_red_scaling_gas               = atmos.scaling_gas

    atmos_red_minor_scales_with_density = atmos.minor_scales_with_density
    atmos_red_scale_by_complement       = atmos.scale_by_complement
    atmos_red_kminor_start              = atmos.kminor_start

    atmos_red_kminor                    = atmos.kminor
    atmos_red_minor_limits_gpt          = atmos.minor_limits_gpt
  else
    atmos_red_minor_gases               = pack(atmos.minor_gases              , gas_is_present)
    atmos_red_scaling_gas               = pack(atmos.scaling_gas              , gas_is_present)

    atmos_red_minor_scales_with_density = pack(atmos.minor_scales_with_density, gas_is_present)
    atmos_red_scale_by_complement       = pack(atmos.scale_by_complement      , gas_is_present)
    atmos_red_kminor_start              = pack(atmos.kminor_start       , gas_is_present)

    atmos_red_kminor = zeros(FT, tot_g, size(atmos.kminor,2), size(atmos.kminor,3))
    atmos_red_minor_limits_gpt = Array{I}(undef, 2, red_nm)

    icnt = 0
    n_elim = 0
    for i = 1:nm
      ng = atmos.minor_limits_gpt[2,i]-atmos.minor_limits_gpt[1,i]+1
      if gas_is_present[i]
        icnt = icnt + 1
        atmos_red_minor_limits_gpt[1:2,icnt] = atmos.minor_limits_gpt[1:2,i]
        atmos_red_kminor_start[icnt] = atmos.kminor_start[i]-n_elim
        for j = 1:ng
          atmos_red_kminor[atmos_red_kminor_start[icnt]+j-1,:,:] = atmos.kminor[atmos.kminor_start[i]+j-1,:,:]
        end
      else
        n_elim = n_elim + ng
      end
    end
  end
  return GasOpticsVars{FT,I}(atmos_red_minor_limits_gpt,
                             atmos_red_minor_scales_with_density,
                             atmos_red_scale_by_complement,
                             atmos_red_kminor_start,
                             atmos_red_kminor,
                             atmos_red_scaling_gas,
                             atmos_red_minor_gases)

end


"""
    key_species_pair2flavor(flavor, key_species_pair)

Flavor index; -1 if not found
"""
function key_species_pair2flavor(flavor::Array{I,2},
                                 key_species_pair::Vector{I}) where {I<:Int}
  iflav = findfirst(mapslices(x->all(key_species_pair .== x), flavor; dims=1)[:])
  return iflav==nothing ? -1 : iflav
end

"""
    create_gpoint_flavor(key_species, gpt2band, flavor)

create gpoint_flavor list a map pointing from each
g-point to the corresponding entry in the "flavor list"

integer, dimension(:,:,:), intent(in) :: key_species
integer, dimension(:), intent(in) :: gpt2band
integer, dimension(:,:), intent(in) :: flavor
integer, dimension(:,:), intent(out) :: gpoint_flavor
"""
function create_gpoint_flavor(key_species, gpt2band, flavor)
  ngpt = length(gpt2band)
  gpoint_flavor = Array{Int}(undef, 2,ngpt)
  for igpt=1:ngpt
    for iatm=1:2
      gpoint_flavor[iatm,igpt] = key_species_pair2flavor(flavor,
                                                         rewrite_key_species_pair(key_species[:,iatm,gpt2band[igpt]])
                                                         )
    end
  end
  return gpoint_flavor
end

#####
##### Sizes of tables: pressure, temperate, η (mixing fraction), nflav
#####

"""
    get_nflav(this::AbstractGasOptics)

Number of distinct major gas pairs in the spectral bands (referred to as
"flavors" - all bands have a flavor even if there is one or no major gas)
"""
get_nflav(this::AbstractGasOptics) = size(this.flavor, 2)

"""
    get_neta(this::AbstractGasOptics)

size of η dimension
"""
get_neta(this::AbstractGasOptics) = size(this.kmajor,2)

"""
    get_npres(this::AbstractGasOptics)

Number of pressures in reference profile absorption
coefficient table is one bigger since a pressure is
repeated in upper/lower atmosphere.
"""
get_npres(this::AbstractGasOptics) = size(this.kmajor,3)-1

"""
    get_ntemp(this::AbstractGasOptics)

Number of temperatures
"""
get_ntemp(this::AbstractGasOptics) = size(this.kmajor,4)

include("GasOptics_kernels.jl")

end #module
