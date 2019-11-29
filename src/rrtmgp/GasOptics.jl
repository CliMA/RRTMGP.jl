"""
    GasOptics

Computes spectrally-resolved gas optical properties and source functions
 given atmospheric physical properties (profiles of temperature, pressure, and gas concentrations)
 The struct must be initialized with data (provided as a netCDF file) before being used.

Two variants apply to internal Planck sources (longwave radiation in the Earth's atmosphere) and to
 external stellar radiation (shortwave radiation in the Earth's atmosphere).
 The variant is chosen based on what information is supplied during initialization.
 (It might make more sense to define two sub-classes)
"""
module GasOptics

using OffsetArrays
using TimerOutputs
using DocStringExtensions

const to_gor = TimerOutput()

using ..FortranIntrinsics
using ..PhysicalConstants
using ..ArrayUtilities
using ..OpticalProps
using ..SourceFunctions
using ..Utilities
using ..GasConcentrations
using ..OpticalProps

export gas_optics_int!, gas_optics_ext!, load_totplnk, load_solar_source
export source_is_internal, source_is_external, get_press_min
export get_col_dry
export GasOpticsVars

include("Reference.jl")

"""
    GasOpticsAux{I}

Auxiliary indexes for variables in the upper
and lower levels of the atmosphere.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GasOpticsAux{I}
  "indexes for determining `col_gas`"
  idx_minor::Vector{I}
  "indexes that have special treatment in density scaling"
  idx_minor_scaling::Vector{I}
  function GasOpticsAux(gas_names_present, gas_minor, identifier_minor, reduced_atmos, ::Type{I}) where I
    return new{I}(create_idx_minor(I,
                                   gas_names_present,
                                   gas_minor,
                                   identifier_minor,
                                   reduced_atmos.minor_gases),
                  create_idx_minor_scaling(I,
                                           gas_names_present,
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
struct GasOpticsVars{FT,I}
  "minor g-point limits"
  minor_limits_gpt::Array{I,2}
  "minor scales with density"
  minor_scales_with_density::Vector{Bool}
  "scale by complement"
  scale_by_complement::Vector{Bool}
  "kminor start"
  kminor_start::Vector{I}
  "absorption coefficient of minor species"
  kminor::Array{FT,3} #, (n_minor,η,temperature)
  "scaling gas"
  scaling_gas::Array{String}
  "minor gases"
  minor_gases::Array{String}
end

abstract type AbstractGasOptics{T,I} <: AbstractOpticalProps{T,I} end

"""
    InternalSourceGasOptics{FT,I} <: AbstractGasOptics{FT,I}

Gas optics with internal sources (for longwave radiation)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct InternalSourceGasOptics{FT,I} <: AbstractGasOptics{FT,I}
  "Reference data"
  ref::Reference{FT}
  "Base optical properties"
  optical_props::OpticalPropsBase{FT,I}
  "GasOpticsVars in the lower atmosphere"
  lower::GasOpticsVars
  "GasOpticsVars in the upper atmosphere"
  upper::GasOpticsVars
  "Auxiliary variables (index maps) in the lower atmosphere"
  lower_aux::GasOpticsAux{I}
  "Auxiliary variables (index maps) in the upper atmosphere"
  upper_aux::GasOpticsAux{I}
  "Present gas names"
  gas_names::Vector{String}
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
    ExternalSourceGasOptics{FT,I} <: AbstractGasOptics{FT,I}

Gas optics with external sources (for shortwave radiation)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ExternalSourceGasOptics{FT,I} <: AbstractGasOptics{FT,I}
  "Reference data"
  ref::Reference{FT}
  "Base optical properties"
  optical_props::OpticalPropsBase{FT,I}
  "GasOpticsVars in the lower atmosphere"
  lower::GasOpticsVars
  "GasOpticsVars in the upper atmosphere"
  upper::GasOpticsVars
  "Auxiliary variables (index maps) in the lower atmosphere"
  lower_aux::GasOpticsAux{I}
  "Auxiliary variables (index maps) in the upper atmosphere"
  upper_aux::GasOpticsAux{I}
  "Present gas names"
  gas_names::Vector{String}
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
struct InterpolationCoefficients{FT,I}
  "index for temperature"
  jtemp::Array{I,2}
  "index for pressure"
  jpress::Array{I,2}
  "index for binary species parameter (η)"
  j_η::Array{I,4}
  "troposphere mask: itropo = merge(1,2,tropo[icol,ilay]); itropo = 1 lower atmosphere; itropo = 2 upper atmosphere"
  tropo::Array{Bool,2}
  "fractions for major species"
  fmajor::Array{FT,6}
  "fractions for minor species"
  fminor::Array{FT,5}
end
function InterpolationCoefficients(::Type{FT}, ::Type{I}, ncol, nlay, nflav) where {I<:Int, FT<:AbstractFloat}
  jtemp = Array{I}(undef, ncol, nlay)
  jpress = Array{I}(undef, ncol, nlay)
  j_η = Array{I}(undef, 2, nflav, ncol, nlay)
  tropo = Array{Bool}(undef, ncol, nlay)
  fmajor = zeros(FT, 2,2,2, nflav,ncol, nlay)
  fminor  = Array{FT}(undef, 2,2, nflav,ncol,nlay)
  return InterpolationCoefficients{FT,I}(jtemp,jpress,j_η,tropo,fmajor,fminor)
end

"""
    get_ngas(this::AbstractGasOptics)

Number of gases registered in the spectral configuration
"""
get_ngas(this::AbstractGasOptics) = length(this.gas_names)


"""
    get_nflav(this::AbstractGasOptics)

Number of distinct major gas pairs in the spectral bands (referred to as
"flavors" - all bands have a flavor even if there is one or no major gas)
"""
get_nflav(this::AbstractGasOptics) = size(this.flavor, 2)


"""
    gas_optics_int!(this::InternalSourceGasOptics,
                     play,
                     plev,
                     tlay,
                     tsfc,
                     gas_desc::GasConcs,
                     optical_props::AbstractOpticalPropsArry,
                     sources::SourceFuncLW;
                     col_dry=nothing,
                     tlev=nothing)

Compute gas optical depth and Planck source functions,
#  given temperature, pressure, and composition

 - `col_dry` Number of molecules per cm-2 of dry air
 - `this` gas optics, see [`InternalSourceGasOptics`](@ref)

inputs
real(FT), dimension(:,:), intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
                                          plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
                                          tlay      # layer temperatures [K]; (ncol,nlay)
real(FT), dimension(:),   intent(in   ) :: tsfc      # surface skin temperatures [K]; (ncol)
type(GasConcs),       intent(in   ) :: gas_desc  # Gas volume mixing ratios
# output
class(AbstractOpticalPropsArry), intent(inout) :: optical_props # Optical properties
class(SourceFuncLW    ), intent(inout) :: sources       # Planck sources
# Optional inputs
real(FT), dimension(:,:),   intent(in   ), optional, target :: col_dry,   # Column dry amount; dim(ncol,nlay)
                                                                tlev        # level temperatures [K]; (ncol,nlay+1)
"""
function gas_optics_int!(this::InternalSourceGasOptics,
                     play,
                     plev,
                     tlay,
                     tsfc,
                     gas_desc::GasConcs,
                     optical_props::AbstractOpticalPropsArry,
                     sources::SourceFuncLW;
                     col_dry=nothing,
                     tlev=nothing)
  FT = eltype(play)

  ncol,nlay  = size(play)
  ics = InterpolationCoefficients(FT, Int, ncol, nlay, get_nflav(this))

  ngpt  = get_ngpt(this.optical_props)
  nband = get_nband(this.optical_props)

  # Gas optics
  compute_gas_τs!(ics, this,
                   ncol, nlay, ngpt, nband,
                   play, plev, tlay, gas_desc,
                   optical_props,
                   col_dry)

  # External source -- check arrays sizes and values
  # input data sizes and values
  check_extent(tsfc, ncol, "tsfc")
  check_range(tsfc, this.ref.temp_min,  this.ref.temp_max,  "tsfc")
  if present(tlev)
    check_extent(tlev, (ncol, nlay+1), "tlev")
    check_range(tlev, this.ref.temp_min, this.ref.temp_max, "tlev")
  end

  #   output extents
  @assert get_ncol(sources) == ncol
  @assert get_nlay(sources) == nlay
  @assert get_ngpt(sources) == ngpt

  # Interpolate source function
  source!(this,
         ncol, nlay, nband, ngpt,
         play, plev, tlay, tsfc,
         ics,
         sources,
         tlev)
  nothing
end

"""
    gas_optics_ext!(this::ExternalSourceGasOptics,
                     play,
                     plev,
                     tlay,
                     gas_desc::GasConcs,    # mandatory inputs
                     optical_props::AbstractOpticalPropsArry,
                     toa_src,        # mandatory outputs
                     col_dry=nothing)

Compute gas optical depth, `toa_src`, given temperature, pressure, and composition

# class(ExternalSourceGasOptics), intent(in) :: this
# real(FT), dimension(:,:), intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
#                                            plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
#                                            tlay      # layer temperatures [K]; (ncol,nlay)
# type(GasConcs),       intent(in   ) :: gas_desc  # Gas volume mixing ratios
# # output
# class(AbstractOpticalPropsArry), intent(inout) :: optical_props
# real(FT), dimension(:,:), intent(  out) :: toa_src     # Incoming solar irradiance(ncol,ngpt)

# Optional inputs
# real(FT), dimension(:,:), intent(in   ), optional, target :: col_dry # Column dry amount; dim(ncol,nlay)
"""
function gas_optics_ext!(this::ExternalSourceGasOptics,
                     play,
                     plev,
                     tlay,
                     gas_desc::GasConcs,    # mandatory inputs
                     optical_props::AbstractOpticalPropsArry,
                     toa_src,        # mandatory outputs
                     col_dry=nothing, last_call=false)

  FT = eltype(play)

  ncol,nlay  = size(play)
  ics = InterpolationCoefficients(FT, Int, ncol,nlay, get_nflav(this))
  ngpt  = get_ngpt(this.optical_props)
  nband = get_nband(this.optical_props)
  ngas  = get_ngas(this)
  nflav = get_nflav(this)

  # Gas optics
  @timeit to_gor "compute_gas_τs!" compute_gas_τs!(ics, this,
                   ncol, nlay, ngpt, nband,
                   play, plev, tlay, gas_desc,
                   optical_props,
                   col_dry, last_call)

  # External source function is constant
  check_extent(toa_src,     (ncol,         ngpt), "toa_src")
  for igpt in 1:ngpt
     for icol in 1:ncol
        toa_src[icol,igpt] = this.solar_src[igpt]
     end
  end
  last_call && @show to_gor
end


# Returns optical properties and interpolation coefficients
"""
    compute_gas_τs!(this::AbstractGasOptics,
                          ncol, nlay, ngpt, nband,
                          play, plev, tlay, gas_desc::GasConcs,
                          optical_props::AbstractOpticalPropsArry,
                          col_dry=nothing)

 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)

class(AbstractGasOptics), intent(in   ) :: this
integer,                          intent(in   ) :: ncol, nlay, ngpt, nband
real(FT), dimension(:,:),         intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
                                                   plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
                                                   tlay      # layer temperatures [K]; (ncol,nlay)
type(GasConcs),               intent(in   ) :: gas_desc  # Gas volume mixing ratios
class(AbstractOpticalPropsArry),     intent(inout) :: optical_props #inout because components are allocated
Optional inputs
real(FT), dimension(:,:), intent(in   ), optional, target :: col_dry # Column dry amount; dim(ncol,nlay)
----------------------------------------------------------
Local variables
real(FT), dimension(ngpt,nlay,ncol) :: τ, τ_Rayleigh  # absorption, Rayleigh scattering optical depths
# integer :: igas, idx_h2o # index of some gases
# Number of molecules per cm^2
real(FT), dimension(ncol,nlay), target  :: col_dry_arr
real(FT), dimension(:,:),       pointer :: col_dry_wk
#
# Interpolation variables used in major gas but not elsewhere, so don't need exporting
#
real(FT), dimension(ncol,nlay,  this%get_ngas()) :: vmr     # volume mixing ratios
real(FT), dimension(ncol,nlay,0:this%get_ngas()) :: col_gas # column amounts for each gas, plus col_dry
real(FT), dimension(2,    get_nflav(this),ncol,nlay) :: col_mix # combination of major species's column amounts
                                                    # index(1) : reference temperature level
                                                    # index(2) : flavor
                                                    # index(3) : layer
real(FT), dimension(2,2,  get_nflav(this),ncol,nlay) :: fminor # interpolation fractions for minor species
                                                     # index(1) : reference η level (temperature dependent)
                                                     # index(2) : reference temperature level
                                                     # index(3) : flavor
                                                     # index(4) : layer
integer :: ngas, nflav, neta, npres, ntemp
----------------------------------------------------------
"""
function compute_gas_τs!(ics::InterpolationCoefficients,
                           this::AbstractGasOptics{FT},
                           ncol, nlay, ngpt, nband,
                           play, plev, tlay, gas_desc::GasConcs,
                           optical_props::AbstractOpticalPropsArry,
                           col_dry=nothing, last_call=false) where FT

  τ          = Array{FT}(undef, ngpt,nlay,ncol)          # absorption, Rayleigh scattering optical depths
  τ_Rayleigh = Array{FT}(undef, ngpt,nlay,ncol) # absorption, Rayleigh scattering optical depths
  col_dry_arr  = Array{FT}(undef, ncol, nlay)
  col_dry_wk   = Array{FT}(undef, ncol, nlay)

  vmr     = Array{FT}(undef, ncol,nlay,  get_ngas(this)) # volume mixing ratios
  col_gas = OffsetArray{FT}(undef, 1:ncol,1:nlay,0:get_ngas(this)) # column amounts for each gas, plus col_dry
  col_mix = Array{FT}(undef, 2,    get_nflav(this),ncol,nlay) # combination of major species's column amounts

  # Check for presence of key species in GasConcs; return error if any key species are not present
  check_key_species_present(this, gas_desc.gas_names)

  # Check input data sizes and values
  check_extent(play, (ncol, nlay  ), "play")
  check_extent(plev, (ncol, nlay+1), "plev")
  check_extent(tlay, (ncol, nlay  ), "tlay")
  check_range(play, this.ref.press_min, this.ref.press_max, "play")
  check_range(plev, this.ref.press_min, this.ref.press_max, "plev")
  check_range(tlay, this.ref.temp_min,  this.ref.temp_max,  "tlay")
  if present(col_dry)
    check_extent(col_dry, (ncol, nlay), "col_dry")
    check_range(col_dry, FT(0), floatmax(FT), "col_dry")
  end

  ngas  = get_ngas(this)
  nflav = get_nflav(this)
  neta  = get_neta(this)
  npres = get_npres(this)
  ntemp = get_ntemp(this)

  # Fill out the array of volume mixing ratios
  for igas in 1:ngas
    # Get vmr if  gas is provided in GasConcs
    if this.gas_names[igas] in gas_desc.gas_names
      get_vmr!(@view(vmr[:,:,igas]), gas_desc, this.gas_names[igas])
    end
  end

  # Compute dry air column amounts (number of molecule per cm^2) if user hasn't provided them
  idx_h2o = loc_in_array("h2o", this.gas_names)
  if present(col_dry)
    col_dry_wk = col_dry
  else
    col_dry_arr = get_col_dry(vmr[:,:,idx_h2o], plev, tlay) # dry air column amounts computation
    col_dry_wk = col_dry_arr
  end

  # compute column gas amounts [molec/cm^2]
  col_gas[1:ncol,1:nlay,0] .= col_dry_wk[1:ncol,1:nlay]
  for igas = 1:ngas
    col_gas[1:ncol,1:nlay,igas] .= vmr[1:ncol,1:nlay,igas] .* col_dry_wk[1:ncol,1:nlay]
  end

  # ---- calculate gas optical depths ----
  τ .= 0
  @timeit to_gor "interpolation!" interpolation!(ics,
          col_mix,
          ncol,nlay,                        # problem dimensions
          nflav, neta, npres, ntemp,  # interpolation dimensions
          this.flavor,
          this.ref,
          play,
          tlay,
          col_gas)
  @timeit to_gor "compute_τ_absorption!" compute_τ_absorption!(τ,
          ncol,nlay,nband,ngpt,                      # dimensions
          idx_h2o,
          this.gpoint_flavor,
          get_band_lims_gpoint(this.optical_props),
          this.kmajor,
          this.lower,
          this.upper,
          this.lower_aux,
          this.upper_aux,
          ics,
          col_mix,
          play,
          tlay,
          col_gas)
  if allocated(this.krayl)

    @timeit to_gor "compute_τ_Rayleigh!" compute_τ_Rayleigh!(          #Rayleigh scattering optical depths
          ncol,nlay,nband,ngpt,# dimensions
          this.gpoint_flavor,
          get_band_lims_gpoint(this.optical_props),
          this.krayl,                   # inputs from object
          idx_h2o,
          col_dry_wk,
          col_gas,
          ics,
          τ_Rayleigh)

  end

  # Combine optical depths and reorder for radiative transfer solver.
  @timeit to_gor "combine_and_reorder!" combine_and_reorder!(τ, τ_Rayleigh, allocated(this.krayl), optical_props)

end

"""
    source!(this::InternalSourceGasOptics,
                ncol, nlay, nbnd, ngpt,
                play, plev, tlay, tsfc,
                ics,
                sources::SourceFuncLW,          # Planck sources
                tlev)

Compute Planck source functions at layer centers and levels

 - `ics` interpolation coefficients, see [`InterpolationCoefficients`](@ref)

# inputs
class(InternalSourceGasOptics),    intent(in ) :: this
integer,                               intent(in   ) :: ncol, nlay, nbnd, ngpt
real(FT), dimension(ncol,nlay),        intent(in   ) :: play   # layer pressures [Pa, mb]
real(FT), dimension(ncol,nlay+1),      intent(in   ) :: plev   # level pressures [Pa, mb]
real(FT), dimension(ncol,nlay),        intent(in   ) :: tlay   # layer temperatures [K]
real(FT), dimension(ncol),             intent(in   ) :: tsfc   # surface skin temperatures [K]
class(SourceFuncLW    ),          intent(inout) :: sources
real(FT), dimension(ncol,nlay+1),      intent(in   ), optional, target :: tlev          # level temperatures [K]
----------------------------------------------------------
integer                                      :: icol, ilay, igpt
real(FT), dimension(ngpt,nlay,ncol)          :: lay_source_t, lev_source_inc_t, lev_source_dec_t
real(FT), dimension(ngpt,     ncol)          :: sfc_source_t
# Variables for temperature at layer edges [K] (ncol, nlay+1)
real(FT), dimension(   ncol,nlay+1), target  :: tlev_arr
real(FT), dimension(:,:),            pointer :: tlev_wk
"""
function source!(this::InternalSourceGasOptics{FT},
                 ncol, nlay, nbnd, ngpt,
                 play, plev, tlay, tsfc,
                 ics::InterpolationCoefficients{FT,I},
                 sources::SourceFuncLW,
                 tlev=nothing) where {FT<:AbstractFloat,I<:Int}

  lay_source_t     = zeros(FT, ngpt,nlay,ncol)
  lev_source_inc_t = zeros(FT, ngpt,nlay,ncol)
  lev_source_dec_t = zeros(FT, ngpt,nlay,ncol)
  sfc_source_t     = zeros(FT, ngpt,ncol)
  tlev_arr         = zeros(FT, ncol,nlay+1)

  # Source function needs temperature at interfaces/levels and at layer centers
  if present(tlev)
    #   Users might have provided these
    tlev_wk = tlev
  else
    tlev_wk = tlev_arr
    #
    # Interpolate temperature to levels if not provided
    #   Interpolation and extrapolation at boundaries is weighted by pressure
    #
    for icol = 1:ncol
       tlev_arr[icol,1] = tlay[icol,1] +
                         (plev[icol,1]-play[icol,1])*
                         (tlay[icol,2]-tlay[icol,1])/
                         (play[icol,2]-play[icol,1])
    end
    for ilay in 2:nlay
      for icol in 1:ncol
         tlev_arr[icol,ilay] = (play[icol,ilay-1]*tlay[icol,ilay-1]*
                               (plev[icol,ilay]-play[icol,ilay]) +
                                play[icol,ilay]*tlay[icol,ilay]*
                               (play[icol,ilay-1]-plev[icol,ilay]))/
                               (plev[icol,ilay]*(play[icol,ilay-1] - play[icol,ilay]))
      end
    end
    for icol = 1:ncol
       tlev_arr[icol,nlay+1] = tlay[icol,nlay] +
                              (plev[icol,nlay+1]-play[icol,nlay])*
                              (tlay[icol,nlay]-tlay[icol,nlay-1])/
                              (play[icol,nlay]-play[icol,nlay-1])
    end
  end

  # Compute internal (Planck) source functions at layers and levels,
  #  which depend on mapping from spectral space that creates k-distribution.
  compute_Planck_source!(ncol,
              nlay,
              nbnd,
              ngpt,
              tlay,
              tlev_wk,
              tsfc,
              fmerge(1,nlay,play[1,1] > play[1,nlay]),
              ics,
              get_gpoint_bands(this.optical_props),
              get_band_lims_gpoint(this.optical_props),
              this.planck_frac,
              this.ref.temp_min,
              this.totplnk_delta,
              this.totplnk,
              this.gpoint_flavor,
              sfc_source_t,
              lay_source_t,
              lev_source_inc_t,
              lev_source_dec_t)

  sources.sfc_source .= convert(Array,transpose(sfc_source_t))

  permutedims!(sources.lay_source, lay_source_t, [3,2,1])
  permutedims!(sources.lev_source_inc, lev_source_inc_t, [3,2,1])
  permutedims!(sources.lev_source_dec, lev_source_dec_t, [3,2,1])
  return nothing
end

#####
##### Initialization
#####

"""
    load_totplnk(totplnk, planck_frac, rayl_lower, rayl_upper, ref, args...)

Initialize object based on data read from netCDF file however the user desires.
Rayleigh scattering tables may or may not be present; this is indicated with allocation status
This interface is for the internal-sources object -- includes Plank functions and fractions
"""
function load_totplnk(totplnk, planck_frac, rayl_lower, rayl_upper, ref, args...)
  abs_coeffs = init_abs_coeffs(args...)

  FT = Float64
  if rayl_lower ≠ nothing && rayl_upper ≠ nothing
    krayl = zeros(FT, size(rayl_lower)...,2)
    krayl[:,:,:,1] = rayl_lower
    krayl[:,:,:,2] = rayl_upper
  else
    krayl = nothing
  end
  totplnk_delta =  (ref.temp_max-ref.temp_min) / (size(totplnk, 1)-1)

  return InternalSourceGasOptics{FT,Int}(ref,
                                         abs_coeffs...,
                                         planck_frac,
                                         totplnk,
                                         totplnk_delta,
                                         krayl)

end

"""
    load_solar_source(solar_src, rayl_lower, rayl_upper, ref, args...)

Initialize object based on data read from netCDF file however the user desires.
Rayleigh scattering tables may or may not be present; this is indicated with allocation status
This interface is for the external-sources object -- includes TOA source function table
"""
function load_solar_source(solar_src, rayl_lower, rayl_upper, ref, args...)
  abs_coeffs = init_abs_coeffs(args...)
  FT = Float64

  if rayl_lower ≠ nothing && rayl_upper ≠ nothing
    krayl = zeros(FT, size(rayl_lower)...,2)
    krayl[:,:,:,1] = rayl_lower
    krayl[:,:,:,2] = rayl_upper
  else
    krayl = nothing
  end

  return ExternalSourceGasOptics{FT,Int}(ref,
                                         abs_coeffs...,
                                         solar_src,
                                         krayl)

end

"""
    init_abs_coeffs(...)

Initialize absorption coefficient arrays

 - `lower` lower atmospheric gas optics variables, see [`GasOpticsVars`](@ref)
 - `upper` upper atmospheric gas optics variables, see [`GasOpticsVars`](@ref)
 - `lower_aux` lower atmospheric gas optics auxiliary variables, see [`GasOpticsAux`](@ref)
 - `upper_aux` upper atmospheric gas optics auxiliary variables, see [`GasOpticsAux`](@ref)

character(len=*), dimension(:),       intent(in) :: gas_names
integer,  dimension(:,:,:),   intent(in) :: key_species
integer,  dimension(:,:),     intent(in) :: band2gpt
real(FT), dimension(:,:),     intent(in) :: band_lims_wavenum
real(FT), dimension(:,:,:,:), intent(in) :: kmajor
"""
function init_abs_coeffs(available_gases::Vector{String},
                         gas_names,
                         key_species,
                         band2gpt,
                         band_lims_wavenum,
                         kmajor,
                         gas_minor,
                         identifier_minor,
                         lower,
                         upper)
  FT = eltype(kmajor)

  # Which gases known to the gas optics are present in the host model (available_gases)?
  gas_names_present = intersect(gas_names, available_gases)

  reduced_lower = reduce_minor_arrays(available_gases, gas_minor, identifier_minor, lower)
  reduced_upper = reduce_minor_arrays(available_gases, gas_minor, identifier_minor, upper)

  lower_aux = GasOpticsAux(gas_names_present, gas_minor, identifier_minor, reduced_lower, Int)
  upper_aux = GasOpticsAux(gas_names_present, gas_minor, identifier_minor, reduced_upper, Int)

  # create flavor list
  # Reduce (remap) key_species list; checks that all key gases are present in incoming
  key_species_red = create_key_species_reduce(gas_names, gas_names_present, key_species)

  flavor = create_flavor(key_species_red)
  optical_props = OpticalPropsBase("AbstractGasOptics optical props", band_lims_wavenum, band2gpt)

  # create gpoint_flavor list
  gpoint_flavor = create_gpoint_flavor(key_species_red, get_gpoint_bands(optical_props), flavor)

  # Which species are key in one or more bands?
  #   flavor is an index into gas_names_present
  is_key = [false for i in 1:length(gas_names_present)]
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
          kmajor,
          flavor,
          gpoint_flavor,
          is_key)

end

"""
    check_key_species_present(this::AbstractGasOptics, gas_names::Vector{String})

Ensure that every key gas required by the k-distribution is
present in the gas concentration object
"""
function check_key_species_present(this::AbstractGasOptics, gas_names::Vector{String})
  key_gas_names = pack(this.gas_names, this.is_key)
  for igas = 1:length(key_gas_names)
    @assert key_gas_names[igas] in gas_names
  end
end

"""
    get_minor_list(this::AbstractGasOptics, gas_desc::GasConcs, ngas, names_spec)

List of minor gases to be used in gas_optics()
Function to define names of key and minor gases to be used by gas_optics().
The final list gases includes those that are defined in gas_optics_specification
and are provided in GasConcs.
"""
function get_minor_list(this::AbstractGasOptics, gas_names::Vector{String}, ngas, names_spec::Vector{String})
  gas_is_present = map(x->x in gas_names, names_spec)
  for igas = 1:get_ngas(this)
    gas_is_present[igas] = names_spec[igas] in gas_names
  end
  return pack(this.gas_names, gas_is_present)
end

#####
##### Inquiry functions
#####

"""
    source_is_internal(this::AbstractGasOptics)

Bool indicating if initialized for internal sources
"""
source_is_internal(::ExternalSourceGasOptics) = false
source_is_internal(::InternalSourceGasOptics) = true

"""
    source_is_external(this::AbstractGasOptics)

Bool indicating if initialized for external sources
"""
source_is_external(::ExternalSourceGasOptics) = true
source_is_external(::InternalSourceGasOptics) = false

"""
    get_gases(this::AbstractGasOptics)

Gas names
"""
get_gases(this::AbstractGasOptics) = this.gas_names

"""
    get_press_min(this::AbstractGasOptics)

Minimum pressure on the interpolation grids
"""
get_press_min(this::AbstractGasOptics) = this.ref.press_min

"""
    get_press_max(this::AbstractGasOptics)

Maximum pressure on the interpolation grids
"""
get_press_max(this::AbstractGasOptics) = this.ref.press_max

"""
    get_temp_min(this::AbstractGasOptics)

Minimum temparature on the interpolation grids
"""
get_temp_min(this::AbstractGasOptics) = this.ref.temp_min

"""
    get_temp_max(this::AbstractGasOptics)

Maximum temperature on the interpolation grids
"""
get_temp_max(this::AbstractGasOptics) = this.ref.temp_max

"""
    get_col_dry(vmr_h2o, plev, tlay, latitude=nothing)

Utility function, provided for user convenience
computes column amounts of dry air using hydrostatic equation

# input
real(FT), dimension(:,:), intent(in) :: vmr_h2o  # volume mixing ratio of water vapor to dry air; (ncol,nlay)
real(FT), dimension(:,:), intent(in) :: plev     # Layer boundary pressures [Pa] (ncol,nlay+1)
real(FT), dimension(:,:), intent(in) :: tlay     # Layer temperatures [K] (ncol,nlay)
real(FT), dimension(:),   optional, intent(in) :: latitude # Latitude [degrees] (ncol)
"""
function get_col_dry(vmr_h2o, plev, tlay, latitude=nothing)
  FT = eltype(plev)

  # first and second term of Helmert formula
  helmert1 = FT(9.80665)
  helmert2 = FT(0.02586)
  # local variables
  g0         = zeros(FT, size(tlay,1)             ) # (ncol)
  delta_plev = zeros(FT, size(tlay,1),size(tlay,2)) # (ncol,nlay)
  m_air      = zeros(FT, size(tlay,1),size(tlay,2)) # average mass of air; (ncol,nlay)

  nlay = size(tlay, 2)
  nlev = size(plev, 2)

  if present(latitude)
    g0 .= helmert1 - helmert2 * cos(FT(2) * π * latitude / FT(180)) # acceleration due to gravity [m/s^2]
  else
    g0 .= grav(FT)
  end
  delta_plev .= abs.(plev[:,1:nlev-1] .- plev[:,2:nlev])

  # Get average mass of moist air per mole of moist air
  m_air .= (m_dry(FT) .+ m_h2o(FT) .* vmr_h2o) ./ (1 .+ vmr_h2o)

  # Hydrostatic equation
  col_dry = zeros(FT, size(tlay,1),size(tlay,2))
  col_dry .= FT(10) .* delta_plev .* avogad(FT) ./ (FT(1000)*m_air .* FT(100) .* spread(g0, 2, nlay))
  col_dry .= col_dry ./ (FT(1) .+ vmr_h2o)
  return col_dry
end

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
    create_idx_minor(::Type{I}, gas_names, gas_minor, identifier_minor, minor_gases_atm) where {I<:Integer}

create index list for extracting col_gas needed for minor gas optical depth calculations

 - `gas_names`
 - `gas_minor`
 - `identifier_minor`
 - `minor_gases_atm`
"""
function create_idx_minor(::Type{I},
                          gas_names::VS,
                          gas_minor::VS,
                          identifier_minor::VS,
                          minor_gases_atm::VS) where {I<:Integer,VS<:Vector{String}}
  # Find identifying string for minor species in list of possible identifiers (e.g. h2o_slf)
  idx_mnr = map(x->loc_in_array(x, identifier_minor), minor_gases_atm)
  # Find name of gas associated with minor species identifier (e.g. h2o)
  return map(x->loc_in_array(gas_minor[x], gas_names), idx_mnr)
end

"""
    create_idx_minor_scaling(::Type{I}, gas_names::Vector{S}, scaling_gas_atm::Vector{S}) where {I<:Integer, S<:String}

Create index for special treatment in density scaling of minor gases

 - `gas_names` gas names
 - `scaling_gas_atm`
"""
create_idx_minor_scaling(::Type{I}, gas_names::VS, scaling_gas_atm::VS) where {I<:Integer, VS<:Vector{String}} =
  map(x->loc_in_array(x, gas_names), scaling_gas_atm)

"""
    create_key_species_reduce(gas_names::VS,
                              gas_names_red::VS,
                              key_species::Array{I,3}) where {VS<:Vector{String},I<:Int}

create flavor list
Reduce (remap) key_species list; checks that all key gases are present in incoming
"""
function create_key_species_reduce(gas_names::VS,
                                   gas_names_red::VS,
                                   key_species::Array{I,3}) where {VS<:Vector{String},I<:Int}
  key_species_red = map(x-> x≠0 ? loc_in_array(gas_names[x],gas_names_red) : x, key_species)
  @assert !any(key_species_red .== -1)
  return key_species_red
end

"""
    reduce_minor_arrays(available_gases::Vector{String},
                        gas_minor,
                        identifier_minor,
                        atmos::GasOpticsVars{FT},
                        atmos_red::GasOpticsVars{FT})

Reduce minor arrays so variables only contain minor gases that are available

 - `available_gases` array of available gases
 - `gas_minor` array of minor gases
 - `identifier_minor`
 - `atmos` original minor `GasOpticsVars` (in)
"""
function reduce_minor_arrays(available_gases::Vector{String},
                             gas_minor,
                             identifier_minor,
                             atmos::GasOpticsVars{FT}) where FT

  I = Int
  nm = length(atmos.minor_gases)
  tot_g=0

  mask = map(x->loc_in_array(x, identifier_minor), atmos.minor_gases)
  gas_is_present = map(x->x in available_gases, gas_minor[mask])

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
integer, dimension(:,:), intent(out), allocatable :: gpoint_flavor
integer :: ngpt, igpt, iatm
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

"""
    combine_and_reorder!(τ::Array{FT,3},
                              τ_Rayleigh::Array{FT,3},
                              has_Rayleigh::Bool,
                              optical_props::AbstractOpticalPropsArry{FT}) where FT

Utility function to combine optical depths from gas absorption and Rayleigh scattering
(and reorder them for convenience, while we're at it)
"""
function combine_and_reorder!(τ::Array{FT,3},
                              τ_Rayleigh::Array{FT,3},
                              has_Rayleigh::Bool,
                              optical_props::AbstractOpticalPropsArry{FT}) where FT
  if !has_Rayleigh
    # index reorder (ngpt, nlay, ncol) -> (ncol,nlay,gpt)
    permutedims!(optical_props.τ, τ, [3,2,1])
    if optical_props isa TwoStream
      optical_props.ssa .= FT(0)
      optical_props.g   .= FT(0)
    end
  else
    # combine optical depth and Rayleigh scattering
    if optical_props isa OneScalar
      # User is asking for absorption optical depth
      permutedims!(optical_props.τ, τ, [3,2,1])

    elseif optical_props isa TwoStream
      combine_and_reorder_2str!(optical_props, τ, τ_Rayleigh)
    end
  end
end

#####
##### Sizes of tables: pressure, temperate, η (mixing fraction)
#####   Equivalent routines for the number of gases and flavors
#####   (get_ngas(), get_nflav()) are defined above because they're
#####   used in function definitions
#####

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

"""
    get_nPlanckTemp(this::InternalSourceGasOptics)

Number of temperatures for Planck function
"""
get_nPlanckTemp(this::InternalSourceGasOptics) = size(this.totplnk,1) # dimensions are Planck-temperature, band

"""
    check_extent(array, s, label)

Assert array has correct size
"""
check_extent(array, s, label) = @assert all(size(array).==s)

"""
    check_range(val, minV, maxV, label)

Assert range of values is valid
"""
check_range(val, minV, maxV, label) = any(val .< minV) || any(val .> maxV) ? error(strip(label) * " values out of range.") : nothing

include("GasOptics_kernels.jl")

end #module
