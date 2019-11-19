"""
    mo_gas_optics_rrtmgp

Computes spectrally-resolved gas optical properties and source functions
 given atmopsheric physical properties (profiles of temperature, pressure, and gas concentrations)
 The struct must be initialized with data (provided as a netCDF file) before being used.

Two variants apply to internal Planck sources (longwave radiation in the Earth's atmosphere) and to
 external stellar radiation (shortwave radiation in the Earth's atmosphere).
 The variant is chosen based on what information is supplied during initialization.
 (It might make more sense to define two sub-classes)
"""
module mo_gas_optics_rrtmgp
using OffsetArrays
using TimerOutputs
const to_gor = TimerOutput()
using ..fortran_intrinsics
using ..mo_rrtmgp_constants
using ..mo_util_array
using ..mo_optical_props
using ..mo_source_functions

using ..Utilities
using ..mo_gas_concentrations
using ..mo_optical_props
export gas_optics_int!, gas_optics_ext!, ty_gas_optics_rrtmgp, load_totplnk, load_solar_source
export source_is_internal, source_is_external, get_press_min

export Reference
export AtmosVars

abstract type ty_gas_optics{T,I} <: ty_optical_props{T,I} end

struct Reference{FT}
  press#::Vector{T}
  press_log#::Vector{T}
  temp#::Vector{T}
  press_min#::T
  press_max#::T
  temp_min#::T
  temp_max#::T
  press_log_delta#::T
  temp_delta#::T
  press_trop_log#::T
  vmr#::Array{T,3}   # vmr(lower or upper atmosphere, gas, temp)
  function Reference(press::Array{FT},
                     temp::Array{FT},
                     press_ref_trop::FT,
                     vmr_ref::Array{FT},
                     available_gases::Array{S},
                     gas_names::Array{S}) where {FT<:AbstractFloat,S<:AbstractString}

    gas_is_present = map(x->x in available_gases, gas_names)
    ngas = count(gas_is_present)

    press_log = log.(press)
    # TODO: remove assumption of ordering
    temp_min = temp[1]
    temp_max = temp[length(temp)]
    press_min = press[length(press)]
    press_max = press[1]

    press_trop_log = log(press_ref_trop)

    press_log_delta = (log(press_min)-log(press_max))/(length(press)-1)
    temp_delta      = (temp_max-temp_min)/(length(temp)-1)

    vmr_ref_red = OffsetArray{FT}(undef, 1:size(vmr_ref, 1),0:ngas, 1:size(vmr_ref, 3))

    # Gas 0 is used in single-key species method, set to 1.0 (col_dry)
    vmr_ref_red[:,0,:] = vmr_ref[:,1,:]
    for i = 1:ngas
      idx = loc_in_array(available_gases[i], gas_names)
      vmr_ref_red[:,i,:] = vmr_ref[:,idx+1,:]
    end
    return new{FT}(press,
                   press_log,
                   temp,
                   press_min,
                   press_max,
                   temp_min,
                   temp_max,
                   press_log_delta,
                   temp_delta,
                   press_trop_log,
                   vmr_ref_red)
  end
end

mutable struct AuxVars{FT}
  idx_minor#::Vector{I}
  idx_minor_scaling#::Vector{I}
end
AuxVars(FT) = AuxVars{FT}(ntuple(i->nothing,2)...)

mutable struct AtmosVars{FT}
  minor_limits_gpt#::Array{I,2}
  minor_scales_with_density#::Vector{Bool}
  scale_by_complement#::Vector{Bool}
  kminor_start#::Vector{I}
  kminor#::Array{T,3}
  scaling_gas
  minor_gases
  # AtmosVars(args::FT...) = AtmosVars{FT}(args..., nothing)
end
AtmosVars(FT) = AtmosVars{FT}(ntuple(i->nothing,7)...)

mutable struct ty_gas_optics_rrtmgp{T,I} <: ty_gas_optics{T,I}
  optical_props#::ty_optical_props
  ref#::Reference
  lower
  upper
  lower_aux
  upper_aux
  gas_names#::Vector{String}     # gas names
  flavor#::Array{I, 2}        # major species pair; (2,nflav)
  gpoint_flavor#::Array{I, 2} # flavor = gpoint_flavor(2, g-point)
  kmajor#::Array{T,4}        #  kmajor(g-point,eta,pressure,temperature)
  krayl#::Array{T, 4} # krayl(g-point,eta,temperature,upper/lower atmosphere)
  planck_frac#::Array{T, 4}   # stored fraction of Planck irradiance in band for given g-point
  totplnk#::Array{T,2}       # integrated Planck irradiance by band; (Planck temperatures,band)
  totplnk_delta#::T # temperature steps in totplnk
  solar_src#::Vector{T} # incoming solar irradiance(g-point)
  is_key#::Vector{Bool}
end

# col_dry is the number of molecules per cm-2 of dry air
export get_col_dry

"""
    get_ngas(this::ty_gas_optics_rrtmgp)

Number of gases registered in the spectral configuration
"""
get_ngas(this::ty_gas_optics_rrtmgp) = length(this.gas_names)


"""
    get_nflav(this::ty_gas_optics_rrtmgp)

Number of distinct major gas pairs in the spectral bands (referred to as
"flavors" - all bands have a flavor even if there is one or no major gas)
"""
get_nflav(this::ty_gas_optics_rrtmgp) = size(this.flavor, 2)


"""
    gas_optics_int!(this::ty_gas_optics_rrtmgp,
                     play,
                     plev,
                     tlay,
                     tsfc,
                     gas_desc::ty_gas_concs,
                     optical_props::ty_optical_props_arry,
                     sources::ty_source_func_lw;
                     col_dry=nothing,
                     tlev=nothing)

Compute gas optical depth and Planck source functions,
#  given temperature, pressure, and composition
inputs
class(ty_gas_optics_rrtmgp), intent(in) :: this
real(FT), dimension(:,:), intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
                                          plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
                                          tlay      # layer temperatures [K]; (ncol,nlay)
real(FT), dimension(:),   intent(in   ) :: tsfc      # surface skin temperatures [K]; (ncol)
type(ty_gas_concs),       intent(in   ) :: gas_desc  # Gas volume mixing ratios
# output
class(ty_optical_props_arry),
                         intent(inout) :: optical_props # Optical properties
class(ty_source_func_lw    ),
                         intent(inout) :: sources       # Planck sources
character(len=128)                      :: error_msg
# Optional inputs
real(FT), dimension(:,:),   intent(in   ),
                      optional, target :: col_dry,   # Column dry amount; dim(ncol,nlay)
                                          tlev        # level temperatures [K]; (ncol,nlay+1)
"""
function gas_optics_int!(this::ty_gas_optics_rrtmgp,
                     play,
                     plev,
                     tlay,
                     tsfc,
                     gas_desc::ty_gas_concs,
                     optical_props::ty_optical_props_arry,
                     sources::ty_source_func_lw;
                     col_dry=nothing,
                     tlev=nothing)
  FT = eltype(play)

  jpress = Array{Int}(undef, size(play))
  jtemp = Array{Int}(undef, size(play))
  tropo = Array{Bool}(undef, size(play))
  fmajor = zeros(FT, 2,2,2,get_nflav(this),size(play)...)
  jeta = Array{Int}(undef, 2,    get_nflav(this), size(play)...)

  ncol  = size(play, 1)
  nlay  = size(play, 2)
  ngpt  = get_ngpt(this.optical_props)
  nband = get_nband(this.optical_props)

  # Gas optics
  compute_gas_taus!(jtemp, jpress, jeta, tropo, fmajor, this,
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
         jtemp, jpress, jeta, tropo, fmajor,
         sources,
         tlev)
  nothing
end

"""
    gas_optics_ext!(this::ty_gas_optics_rrtmgp,
                     play,
                     plev,
                     tlay,
                     gas_desc::ty_gas_concs,    # mandatory inputs
                     optical_props::ty_optical_props_arry,
                     toa_src,        # mandatory outputs
                     col_dry=nothing)

Compute gas optical depth, `toa_src`, given temperature, pressure, and composition

# class(ty_gas_optics_rrtmgp), intent(in) :: this
# real(FT), dimension(:,:), intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
#                                            plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
#                                            tlay      # layer temperatures [K]; (ncol,nlay)
# type(ty_gas_concs),       intent(in   ) :: gas_desc  # Gas volume mixing ratios
# # output
# class(ty_optical_props_arry),
#                           intent(inout) :: optical_props
# real(FT), dimension(:,:), intent(  out) :: toa_src     # Incoming solar irradiance(ncol,ngpt)
# character(len=128)                      :: error_msg

# Optional inputs
# real(FT), dimension(:,:), intent(in   ),
#                        optional, target :: col_dry # Column dry amount; dim(ncol,nlay)
"""
function gas_optics_ext!(this::ty_gas_optics_rrtmgp,
                     play,
                     plev,
                     tlay,
                     gas_desc::ty_gas_concs,    # mandatory inputs
                     optical_props::ty_optical_props_arry,
                     toa_src,        # mandatory outputs
                     col_dry=nothing, last_call=false)

  FT = eltype(play)

  jpress = Array{Int}( undef, size(play))
  jtemp  = Array{Int}( undef, size(play))
  tropo  = Array{Bool}(undef, size(play))
  fmajor = zeros(FT, 2,2,2, get_nflav(this), size(play)...)
  jeta   = Array{Int}( undef, 2,     get_nflav(this), size(play)...)

  ncol  = size(play, 1)
  nlay  = size(play, 2)
  ngpt  = get_ngpt(this.optical_props)
  nband = get_nband(this.optical_props)
  ngas  = get_ngas(this)
  nflav = get_nflav(this)

  # Gas optics
  @timeit to_gor "compute_gas_taus!" compute_gas_taus!(jtemp, jpress, jeta, tropo, fmajor, this,
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
  last_call && println("--------------- last_call")
  last_call && @show to_gor
end


# Returns optical properties and interpolation coefficients
"""
    compute_gas_taus!(this::ty_gas_optics_rrtmgp,
                          ncol, nlay, ngpt, nband,
                          play, plev, tlay, gas_desc::ty_gas_concs,
                          optical_props::ty_optical_props_arry,
                          col_dry=nothing)

class(ty_gas_optics_rrtmgp),
                                 intent(in   ) :: this
integer,                          intent(in   ) :: ncol, nlay, ngpt, nband
real(FT), dimension(:,:),         intent(in   ) :: play,    # layer pressures [Pa, mb]; (ncol,nlay)
                                                  plev,    # level pressures [Pa, mb]; (ncol,nlay+1)
                                                  tlay      # layer temperatures [K]; (ncol,nlay)
type(ty_gas_concs),               intent(in   ) :: gas_desc  # Gas volume mixing ratios
class(ty_optical_props_arry),     intent(inout) :: optical_props #inout because components are allocated
# Interpolation coefficients for use in internal source function
integer,     dimension(                       ncol, nlay), intent(  out) :: jtemp, jpress
integer,     dimension(2,    get_nflav(this),ncol, nlay), intent(  out) :: jeta
logical(wl), dimension(                       ncol, nlay), intent(  out) :: tropo
real(FT),    dimension(2,2,2,get_nflav(this),ncol, nlay), intent(  out) :: fmajor
character(len=128)                                         :: error_msg
Optional inputs
real(FT), dimension(:,:), intent(in   ),
                      optional, target :: col_dry # Column dry amount; dim(ncol,nlay)
----------------------------------------------------------
Local variables
real(FT), dimension(ngpt,nlay,ncol) :: tau, tau_rayleigh  # absorption, Rayleigh scattering optical depths
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
                                                     # index(1) : reference eta level (temperature dependent)
                                                     # index(2) : reference temperature level
                                                     # index(3) : flavor
                                                     # index(4) : layer
integer :: ngas, nflav, neta, npres, ntemp
integer :: nminorlower, nminorklower,nminorupper, nminorkupper
logical :: use_rayl
----------------------------------------------------------
"""
function compute_gas_taus!(jtemp, jpress, jeta, tropo, fmajor, this::ty_gas_optics_rrtmgp{FT},
                          ncol, nlay, ngpt, nband,
                          play, plev, tlay, gas_desc::ty_gas_concs,
                          optical_props::ty_optical_props_arry,
                          col_dry=nothing, last_call=false) where FT

  tau          = Array{FT}(undef, ngpt,nlay,ncol)          # absorption, Rayleigh scattering optical depths
  tau_rayleigh = Array{FT}(undef, ngpt,nlay,ncol) # absorption, Rayleigh scattering optical depths
  col_dry_arr  = Array{FT}(undef, ncol, nlay)
  col_dry_wk   = Array{FT}(undef, ncol, nlay)

  vmr     = Array{FT}(undef, ncol,nlay,  get_ngas(this)) # volume mixing ratios
  col_gas = OffsetArray{FT}(undef, 1:ncol,1:nlay,0:get_ngas(this)) # column amounts for each gas, plus col_dry
  col_mix = Array{FT}(undef, 2,    get_nflav(this),ncol,nlay) # combination of major species's column amounts
  fminor  = Array{FT}(undef, 2,2,  get_nflav(this),ncol,nlay) # interpolation fractions for minor species

  # Error checking
  use_rayl = allocated(this.krayl)

  # Check for presence of key species in ty_gas_concs; return error if any key species are not present
  check_key_species_present(this, gas_desc)

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
  # number of minor contributors, total num absorption coeffs
  nminorlower  = size(this.lower.minor_scales_with_density)
  nminorklower = size(this.lower.kminor, 1)
  nminorupper  = size(this.upper.minor_scales_with_density)
  nminorkupper = size(this.upper.kminor, 1)

  # Fill out the array of volume mixing ratios
  for igas in 1:ngas
    # Get vmr if  gas is provided in ty_gas_concs
    if this.gas_names[igas] in gas_desc.gas_name
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
  tau .= 0
  @timeit to_gor "interpolation!" interpolation!(jtemp,fmajor,fminor,col_mix,tropo,jeta,jpress,
          ncol,nlay,                        # problem dimensions
          ngas, nflav, neta, npres, ntemp,  # interpolation dimensions
          this.flavor,
          this.ref,
          play,
          tlay,
          col_gas)
  @timeit to_gor "compute_tau_absorption!" compute_tau_absorption!(tau,
          ncol,nlay,nband,ngpt,                      # dimensions
          ngas,nflav,neta,npres,ntemp,
          nminorlower, nminorklower,                # number of minor contributors, total num absorption coeffs
          nminorupper, nminorkupper,
          idx_h2o,
          this.gpoint_flavor,
          get_band_lims_gpoint(this.optical_props),
          this.kmajor,
          this.lower,
          this.upper,
          this.lower_aux,
          this.upper_aux,
          tropo,
          col_mix,fmajor,fminor,
          play,tlay,col_gas,
          jeta,jtemp,jpress)
  if allocated(this.krayl)

    @timeit to_gor "compute_tau_rayleigh!" compute_tau_rayleigh!(          #Rayleigh scattering optical depths
          ncol,nlay,nband,ngpt,
          ngas,nflav,neta,npres,ntemp,  # dimensions
          this.gpoint_flavor,
          get_band_lims_gpoint(this.optical_props),
          this.krayl,                   # inputs from object
          idx_h2o, col_dry_wk,col_gas,
          fminor,jeta,tropo,jtemp,      # local input
          tau_rayleigh)

  end

  # Combine optical depths and reorder for radiative transfer solver.
  @timeit to_gor "combine_and_reorder!" combine_and_reorder!(tau, tau_rayleigh, allocated(this.krayl), optical_props)

end

"""
    source!(this::ty_gas_optics_rrtmgp,
                ncol, nlay, nbnd, ngpt,
                play, plev, tlay, tsfc,
                jtemp, jpress, jeta, tropo, fmajor,
                sources::ty_source_func_lw,          # Planck sources
                tlev)

Compute Planck source functions at layer centers and levels

# inputs
class(ty_gas_optics_rrtmgp),    intent(in ) :: this
integer,                               intent(in   ) :: ncol, nlay, nbnd, ngpt
real(FT), dimension(ncol,nlay),        intent(in   ) :: play   # layer pressures [Pa, mb]
real(FT), dimension(ncol,nlay+1),      intent(in   ) :: plev   # level pressures [Pa, mb]
real(FT), dimension(ncol,nlay),        intent(in   ) :: tlay   # layer temperatures [K]
real(FT), dimension(ncol),             intent(in   ) :: tsfc   # surface skin temperatures [K]
# Interplation coefficients
integer,     dimension(ncol,nlay),     intent(in   ) :: jtemp, jpress
logical(wl), dimension(ncol,nlay),     intent(in   ) :: tropo
real(FT),    dimension(2,2,2,get_nflav(this),ncol,nlay),
                                      intent(in   ) :: fmajor
integer,     dimension(2,    get_nflav(this),ncol,nlay),
                                      intent(in   ) :: jeta
class(ty_source_func_lw    ),          intent(inout) :: sources
real(FT), dimension(ncol,nlay+1),      intent(in   ),
                                 optional, target :: tlev          # level temperatures [K]
character(len=128)                                 :: error_msg
----------------------------------------------------------
integer                                      :: icol, ilay, igpt
real(FT), dimension(ngpt,nlay,ncol)          :: lay_source_t, lev_source_inc_t, lev_source_dec_t
real(FT), dimension(ngpt,     ncol)          :: sfc_source_t
# Variables for temperature at layer edges [K] (ncol, nlay+1)
real(FT), dimension(   ncol,nlay+1), target  :: tlev_arr
real(FT), dimension(:,:),            pointer :: tlev_wk
"""
function source!(this::ty_gas_optics_rrtmgp,
                ncol, nlay, nbnd, ngpt,
                play, plev, tlay, tsfc,
                jtemp, jpress, jeta, tropo, fmajor,
                sources::ty_source_func_lw,          # Planck sources
                tlev)                                # optional input
  FT = eltype(this.ref.vmr) # Float64

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
                           (plev[icol,1]-play[icol,1])*(tlay[icol,2]-tlay[icol,1]) /
                                                       (play[icol,2]-play[icol,1])
    end
    for ilay in 2:nlay
      for icol in 1:ncol
         tlev_arr[icol,ilay] = (play[icol,ilay-1]*tlay[icol,ilay-1]*(plev[icol,ilay]-play[icol,ilay]) +
                                play[icol,ilay]*tlay[icol,ilay]*(play[icol,ilay-1]-plev[icol,ilay])) /
                                (plev[icol,ilay]*(play[icol,ilay-1] - play[icol,ilay]))
      end
    end
    for icol = 1:ncol
       tlev_arr[icol,nlay+1] = tlay[icol,nlay] +
                               (plev[icol,nlay+1]-play[icol,nlay])*(tlay[icol,nlay]-tlay(icol,nlay-1)) /
                                                                     (play(icol,nlay)-play(icol,nlay-1))
    end
  end

  # Compute internal (Planck) source functions at layers and levels,
  #  which depend on mapping from spectral space that creates k-distribution.
  compute_Planck_source!(ncol, nlay, nbnd, ngpt,
              get_nflav(this), get_neta(this), get_npres(this), get_ntemp(this), get_nPlanckTemp(this),
              tlay, tlev_wk, tsfc, fmerge(1,nlay,play[1,1] > play[1,nlay]),
              fmajor, jeta, tropo, jtemp, jpress,
              get_gpoint_bands(this.optical_props), get_band_lims_gpoint(this.optical_props), this.planck_frac, this.ref.temp_min,
              this.totplnk_delta, this.totplnk, this.gpoint_flavor,
              sfc_source_t, lay_source_t, lev_source_inc_t, lev_source_dec_t)

  sources.sfc_source .= convert(Array,transpose(sfc_source_t))

  permutedims!(sources.lay_source, lay_source_t, [3,2,1])
  permutedims!(sources.lev_source_inc, lev_source_inc_t, [3,2,1])
  permutedims!(sources.lev_source_dec, lev_source_dec_t, [3,2,1])
end

#--------------------------------------------------------------------------------------------------------------------
# Initialization
#--------------------------------------------------------------------------------------------------------------------
# Initialize object based on data read from netCDF file however the user desires.
#  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
# This interface is for the internal-sources object -- includes Plank functions and fractions
#
function load_totplnk(totplnk, planck_frac, rayl_lower, rayl_upper, args...)
  this = init_abs_coeffs(rayl_lower, rayl_upper, args...)
  # Planck function tables
  #
  this.totplnk = totplnk
  this.planck_frac = planck_frac
  # Temperature steps for Planck function interpolation
  #   Assumes that temperature minimum and max are the same for the absorption coefficient grid and the
  #   Planck grid and the Planck grid is equally spaced
  this.totplnk_delta =  (this.ref.temp_max-this.ref.temp_min) / (size(this.totplnk, 1)-1)
  return this
end

#--------------------------------------------------------------------------------------------------------------------
#
# Initialize object based on data read from netCDF file however the user desires.
#  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
# This interface is for the external-sources object -- includes TOA source function table
#
function load_solar_source(solar_src, rayl_lower, rayl_upper, args...)
  this = init_abs_coeffs(rayl_lower, rayl_upper, args...)
  #
  # Solar source table init
  #
  this.solar_src = solar_src
  return this

end
#--------------------------------------------------------------------------------------------------------------------
#
# Initialize absorption coefficient arrays,
#   including Rayleigh scattering tables if provided (allocated)
#
"""
    init_abs_coeffs(...)

class(ty_gas_optics_rrtmgp), intent(inout) :: this
class(ty_gas_concs),                intent(in   ) :: available_gases
character(len=*),
         dimension(:),       intent(in) :: gas_names
integer,  dimension(:,:,:),   intent(in) :: key_species
integer,  dimension(:,:),     intent(in) :: band2gpt
real(FT), dimension(:,:),     intent(in) :: band_lims_wavenum
real(FT), dimension(:),       intent(in) :: press_ref, temp_ref
real(FT),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
real(FT), dimension(:,:,:),   intent(in) :: vmr_ref
real(FT), dimension(:,:,:,:), intent(in) :: kmajor
real(FT), dimension(:,:,:),   intent(in) :: kminor_lower, kminor_upper
character(len=*),   dimension(:),
                             intent(in) :: gas_minor,
                                           identifier_minor
character(len=*),   dimension(:),
                             intent(in) :: minor_gases_lower,
                                           minor_gases_upper
integer,  dimension(:,:),     intent(in) :: minor_limits_gpt_lower,
                                           minor_limits_gpt_upper
logical(wl), dimension(:),    intent(in) :: lower.minor_scales_with_density,
                                           upper.minor_scales_with_density
character(len=*),   dimension(:),
                             intent(in) :: scaling_gas_lower,
                                           scaling_gas_upper
logical(wl), dimension(:),    intent(in) :: lower.scale_by_complement,
                                           upper.scale_by_complement
integer,  dimension(:),       intent(in) :: lower.kminor_start,
                                           upper.kminor_start
real(FT), dimension(:,:,:),   intent(in),
                            allocatable :: rayl_lower, rayl_upper
character(len=128)                       :: err_message
--------------------------------------------------------------------------
logical,  dimension(:),     allocatable :: gas_is_present
logical,  dimension(:),     allocatable :: key_species_present_init
integer,  dimension(:,:,:), allocatable :: key_species_red
real(FT), dimension(:,:,:), allocatable :: vmr_ref_red
character(len=256),
         dimension(:),     allocatable :: minor_gases_lower_red,
                                          minor_gases_upper_red
character(len=256),
         dimension(:),     allocatable :: scaling_gas_lower_red,
                                          scaling_gas_upper_red
integer :: i, j, idx
integer :: ngas
# --------------------------------------
"""
function init_abs_coeffs(rayl_lower, rayl_upper,
available_gases::ty_gas_concs,
gas_names,
key_species,
band2gpt,
band_lims_wavenum,
ref,
kmajor,
gas_minor,
identifier_minor,
lower,
upper)
  FT = eltype(kmajor)

  this = ty_gas_optics_rrtmgp{FT,Int}(
    ty_optical_props_base("ty_gas_optics_rrtmgp optical props", band_lims_wavenum, band2gpt),
    ref,
    AtmosVars(FT),
    AtmosVars(FT),
    AuxVars(FT),
    AuxVars(FT),
    ntuple(i->nothing, 10)...)
  #
  # Which gases known to the gas optics are present in the host model (available_gases)?
  #
  ngas = length(gas_names)
  gas_is_present = Vector{Bool}(undef, ngas...)

  for i in 1:ngas
    gas_is_present[i] = gas_names[i] in available_gases.gas_name
  end
  #
  # Now the number of gases is the union of those known to the k-distribution and provided
  #   by the host model
  #
  ngas = count(gas_is_present)
  #
  # Initialize the gas optics object, keeping only those gases known to the
  #   gas optics and also present in the host model
  #
  this.gas_names = pack(gas_names, gas_is_present)

  #
  # Reduce minor arrays so variables only contain minor gases that are available
  # Reduce size of minor Arrays
  #

  this.lower.kminor,
  minor_gases_lower_red,
  this.lower.minor_limits_gpt,
  this.lower.minor_scales_with_density,
  scaling_gas_lower_red,
  this.lower.scale_by_complement,
  this.lower.kminor_start              = reduce_minor_arrays(available_gases,
gas_names,
gas_minor,
identifier_minor,
lower.kminor,
lower.minor_gases,
lower.minor_limits_gpt,
lower.minor_scales_with_density,
lower.scaling_gas,
lower.scale_by_complement,
lower.kminor_start
)

  this.upper.kminor,
  minor_gases_upper_red,
  this.upper.minor_limits_gpt,
  this.upper.minor_scales_with_density,
  scaling_gas_upper_red,
  this.upper.scale_by_complement,
  this.upper.kminor_start = reduce_minor_arrays(available_gases,
                           gas_names,
                           gas_minor,
                           identifier_minor,
                           upper.kminor,
                           upper.minor_gases,
                           upper.minor_limits_gpt,
                           upper.minor_scales_with_density,
                           upper.scaling_gas,
                           upper.scale_by_complement,
                           upper.kminor_start
                           )

  # Arrays not reduced by the presence, or lack thereof, of a gas
  this.kmajor    = kmajor
  FT = eltype(kmajor)

  @assert allocated(rayl_lower) == allocated(rayl_upper)

  if allocated(rayl_lower)
    this.krayl = zeros(FT, size(rayl_lower,1),size(rayl_lower,2),size(rayl_lower,3),2)
    this.krayl[:,:,:,1] = rayl_lower
    this.krayl[:,:,:,2] = rayl_upper
  end

  # ---- post processing ----

  # Get index of gas (if present) for determining col_gas
  this.lower_aux.idx_minor = create_idx_minor(this.gas_names, gas_minor, identifier_minor, minor_gases_lower_red)
  this.upper_aux.idx_minor = create_idx_minor(this.gas_names, gas_minor, identifier_minor, minor_gases_upper_red)
  # Get index of gas (if present) that has special treatment in density scaling
  this.lower_aux.idx_minor_scaling = create_idx_minor_scaling(this.gas_names, scaling_gas_lower_red)
  this.upper_aux.idx_minor_scaling = create_idx_minor_scaling(this.gas_names, scaling_gas_upper_red)

  # create flavor list
  # Reduce (remap) key_species list; checks that all key gases are present in incoming
  key_species_red,key_species_present_init = create_key_species_reduce(gas_names,this.gas_names, key_species)
  check_key_species_present_init(gas_names,key_species_present_init)

  # create flavor list
  this.flavor = create_flavor(key_species_red)
  # create gpoint_flavor list
  this.gpoint_flavor = create_gpoint_flavor(key_species_red, get_gpoint_bands(this.optical_props), this.flavor)

  # Which species are key in one or more bands?
  #   this%flavor is an index into this%gas_names
  #
  this.is_key = [false for i in 1:get_ngas(this)]
  for j in 1:size(this.flavor, 2)
    for i in 1:size(this.flavor, 1) # should be 2
      if this.flavor[i,j] ≠ 0
        this.is_key[this.flavor[i,j]] = true
      end
    end
  end
  return this

end

check_key_species_present_init(gas_names, key_species_present_init) = @assert all(key_species_present_init)

# Ensure that every key gas required by the k-distribution is
#    present in the gas concentration object
function check_key_species_present(this::ty_gas_optics_rrtmgp, gas_desc::ty_gas_concs)
  # class(ty_gas_optics_rrtmgp), intent(in) :: this
  # class(ty_gas_concs),                intent(in) :: gas_desc
  key_gas_names = pack(this.gas_names, this.is_key)
  for igas = 1:length(key_gas_names)
    @assert key_gas_names[igas] in gas_desc.gas_name
  end
end

# Function to define names of key and minor gases to be used by gas_optics().
# The final list gases includes those that are defined in gas_optics_specification
# and are provided in ty_gas_concs.
#
function get_minor_list(this::ty_gas_optics_rrtmgp, gas_desc::ty_gas_concs, ngas, names_spec)
  # class(ty_gas_optics_rrtmgp), intent(in)       :: this
  # class(ty_gas_concs), intent(in)                      :: gas_desc
  # integer, intent(in)                                  :: ngas
  # character(32), dimension(ngas), intent(in)           :: names_spec

  # # List of minor gases to be used in gas_optics()
  # character(len=32), dimension(:), allocatable         :: get_minor_list
  # # Logical flag for minor species in specification (T = minor; F = not minor)
  # logical, dimension(size(names_spec))                 :: gas_is_present
  # integer                                              :: igas, icnt

  for igas = 1:get_ngas(this)
    gas_is_present[igas] = names_spec[igas] in gas_desc.gas_name
  end
  icnt = count(gas_is_present)
  minor_list = Vector{String}(undef, icnt)
  minor_list .= pack(this.gas_names, gas_is_present)
  return minor_list
end

#### Inquiry functions

"""
    source_is_internal(this::ty_gas_optics_rrtmgp)

Bool indicating if initialized for internal sources
"""
source_is_internal(this::ty_gas_optics_rrtmgp) = allocated(this.totplnk) && allocated(this.planck_frac)

"""
    source_is_external(this::ty_gas_optics_rrtmgp)

Bool indicating if initialized for external sources
"""
source_is_external(this::ty_gas_optics_rrtmgp) = allocated(this.solar_src)

"""
    get_gases(this::ty_gas_optics_rrtmgp)

Gas names
"""
get_gases(this::ty_gas_optics_rrtmgp) = this.gas_names

"""
    get_press_min(this::ty_gas_optics_rrtmgp)

Minimum pressure on the interpolation grids
"""
get_press_min(this::ty_gas_optics_rrtmgp) = this.ref.press_min

"""
    get_press_max(this::ty_gas_optics_rrtmgp)

Maximum pressure on the interpolation grids
"""
get_press_max(this::ty_gas_optics_rrtmgp) = this.ref.press_max

"""
    get_temp_min(this::ty_gas_optics_rrtmgp)

Minimum temparature on the interpolation grids
"""
get_temp_min(this::ty_gas_optics_rrtmgp) = this.ref.temp_min

"""
    get_temp_max(this::ty_gas_optics_rrtmgp)

Maximum temparature on the interpolation grids
"""
get_temp_max(this::ty_gas_optics_rrtmgp) = this.ref.temp_max

"""
    get_col_dry(vmr_h2o, plev, tlay, latitude=nothing)

Utility function, provided for user convenience
computes column amounts of dry air using hydrostatic equation

# input
real(FT), dimension(:,:), intent(in) :: vmr_h2o  # volume mixing ratio of water vapor to dry air; (ncol,nlay)
real(FT), dimension(:,:), intent(in) :: plev     # Layer boundary pressures [Pa] (ncol,nlay+1)
real(FT), dimension(:,:), intent(in) :: tlay     # Layer temperatures [K] (ncol,nlay)
real(FT), dimension(:),   optional,
                         intent(in) :: latitude # Latitude [degrees] (ncol)
# output
real(FT), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: col_dry # Column dry amount (ncol,nlay)
------------------------------------------------
# first and second term of Helmert formula
real(FT), parameter :: helmert1 = 9.80665_wp
real(FT), parameter :: helmert2 = 0.02586_wp
# local variables
real(FT), dimension(size(tlay,dim=1)                 ) :: g0 # (ncol)
real(FT), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: delta_plev # (ncol,nlay)
real(FT), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: m_air # average mass of air; (ncol,nlay)
integer :: nlev, nlay
# ------------------------------------------------
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

#
# Internal procedures
#
function rewrite_key_species_pair(key_species_pair)
  # (0,0) becomes (2,2) -- because absorption coefficients for these g-points will be 0.
  # integer, dimension(2) :: rewrite_key_species_pair
  # integer, dimension(2), intent(in) :: key_species_pair
  result = key_species_pair
  if all(key_species_pair .== [0,0])
    result .= [2,2]
  end
  return result
end

# ---------------------------------------------------------------------------------------
# true is key_species_pair exists in key_species_list
function key_species_pair_exists(key_species_list, key_species_pair)
  # logical                             :: key_species_pair_exists
  # integer, dimension(:,:), intent(in) :: key_species_list
  # integer, dimension(2),   intent(in) :: key_species_pair
  # integer :: i
  for i=1:size(key_species_list,2)
    if all(key_species_list[:,i] .== key_species_pair)
      result = true
      return result
    end
  end
  result = false
  return result
end
# ---------------------------------------------------------------------------------------
# create flavor list --
#   an unordered array of extent (2,:) containing all possible pairs of key species
#   used in either upper or lower atmos
#
function create_flavor(key_species)
  # integer, dimension(:,:,:), intent(in) :: key_species
  # integer, dimension(:,:), allocatable, intent(out) :: flavor
  # integer, dimension(2,size(key_species,3)*2) :: key_species_list
  # integer :: ibnd, iatm, i, iflavor

  key_species_list = Array{Int}(undef, 2,size(key_species,3)*2)

  # prepare list of key_species
  i = 1
  for ibnd=1:size(key_species,3)
    for iatm=1:size(key_species,1)
      key_species_list[:,i] .= key_species[:,iatm,ibnd]
      i = i + 1
    end
  end
  # rewrite single key_species pairs
  for i=1:size(key_species_list,2)
      key_species_list[:,i] = rewrite_key_species_pair(key_species_list[:,i])
  end
  # count unique key species pairs
  iflavor = 0
  for i=1:size(key_species_list,2)
    if !key_species_pair_exists(key_species_list[:,1:i-1],key_species_list[:,i])
      iflavor = iflavor + 1
    end
  end
  # fill flavors
  flavor = Array{Int}(undef, 2,iflavor)
  iflavor = 0
  for i=1:size(key_species_list,2)
    if !key_species_pair_exists(key_species_list[:,1:i-1],key_species_list[:,i])
      iflavor = iflavor + 1
      flavor[:,iflavor] = key_species_list[:,i]
    end
  end
  return flavor
end
# ---------------------------------------------------------------------------------------
#
# create index list for extracting col_gas needed for minor gas optical depth calculations
#
function create_idx_minor(gas_names, gas_minor, identifier_minor, minor_gases_atm)
  # character(len=*), dimension(:), intent(in) :: gas_names
  # character(len=*), dimension(:), intent(in) ::
  #                                               gas_minor,
  #                                               identifier_minor
  # character(len=*), dimension(:), intent(in) :: minor_gases_atm
  # integer, dimension(:), allocatable,
  #                                intent(out) :: idx_minor_atm

  # # local
  # integer :: imnr
  # integer :: idx_mnr
  idx_minor_atm = Vector{Int}(undef, size(minor_gases_atm,1))
  for imnr = 1:size(minor_gases_atm,1) # loop over minor absorbers in each band
        # Find identifying string for minor species in list of possible identifiers (e.g. h2o_slf)
        idx_mnr     = loc_in_array(minor_gases_atm[imnr], identifier_minor)
        # Find name of gas associated with minor species identifier (e.g. h2o)
        idx_minor_atm[imnr] = loc_in_array(gas_minor[idx_mnr],    gas_names)
  end
  return idx_minor_atm

end

# ---------------------------------------------------------------------------------------
#
# create index for special treatment in density scaling of minor gases
#
function create_idx_minor_scaling(gas_names, scaling_gas_atm)
  # character(len=*), dimension(:), intent(in) :: gas_names
  # character(len=*), dimension(:), intent(in) :: scaling_gas_atm
  # integer, dimension(:), allocatable,
  #                                intent(out) :: idx_minor_scaling_atm

  # # local
  # integer :: imnr
  idx_minor_scaling_atm = Vector{Int}(undef, size(scaling_gas_atm,1))
  for imnr = 1:size(scaling_gas_atm,1) # loop over minor absorbers in each band
        # This will be -1 if there's no interacting gas
        idx_minor_scaling_atm[imnr] = loc_in_array(scaling_gas_atm[imnr], gas_names)
  end
  return idx_minor_scaling_atm

end
# ---------------------------------------------------------------------------------------
function create_key_species_reduce(gas_names, gas_names_red, key_species)
  # character(len=*),
  #           dimension(:),       intent(in) :: gas_names
  # character(len=*),
  #           dimension(:),       intent(in) :: gas_names_red
  # integer,  dimension(:,:,:),   intent(in) :: key_species
  # integer,  dimension(:,:,:), allocatable, intent(out) :: key_species_red

  # logical, dimension(:), allocatable, intent(out) :: key_species_present_init
  # integer :: ip, ia, it, np, na, nt

  np = size(key_species,1)
  na = size(key_species,2)
  nt = size(key_species,3)
  key_species_red = Array{Int}(undef, size(key_species,1),
                                      size(key_species,2),
                                      size(key_species,3))
  key_species_present_init = Vector{Bool}(undef, size(gas_names))
  key_species_present_init = true

  for ip = 1:np
    for ia = 1:na
      for it = 1:nt
        if key_species[ip,ia,it] ≠ 0
          key_species_red[ip,ia,it] = loc_in_array(gas_names[key_species[ip,ia,it]],gas_names_red)
          if key_species_red[ip,ia,it] == -1
            key_species_present_init[key_species[ip,ia,it]] = false
          end
        else
          key_species_red[ip,ia,it] = key_species[ip,ia,it]
        end
      end
    end
  end
  return key_species_red,key_species_present_init

end

# ---------------------------------------------------------------------------------------
function reduce_minor_arrays(available_gases::ty_gas_concs,
gas_names,
gas_minor,
identifier_minor,
kminor_atm,
minor_gases_atm,
minor_limits_gpt_atm,
minor_scales_with_density_atm,
scaling_gas_atm,
scale_by_complement_atm,
kminor_start_atm
# kminor_atm_red,
# minor_gases_atm_red,
# minor_limits_gpt_atm_red,
# minor_scales_with_density_atm_red,
# scaling_gas_atm_red,
# scale_by_complement_atm_red,
# kminor_start_atm_red
)

  # class(ty_gas_concs),                intent(in   ) :: available_gases
  # character(len=*), dimension(:),     intent(in) :: gas_names
  # real(FT),         dimension(:,:,:), intent(in) :: kminor_atm
  # character(len=*), dimension(:),     intent(in) :: gas_minor,
  #                                                   identifier_minor
  # character(len=*), dimension(:),     intent(in) :: minor_gases_atm
  # integer,          dimension(:,:),   intent(in) :: minor_limits_gpt_atm
  # logical(wl),      dimension(:),     intent(in) :: minor_scales_with_density_atm
  # character(len=*), dimension(:),     intent(in) :: scaling_gas_atm
  # logical(wl),      dimension(:),     intent(in) :: scale_by_complement_atm
  # integer,          dimension(:),     intent(in) :: kminor_start_atm
  # real(FT),         dimension(:,:,:), allocatable, intent(out) :: kminor_atm_red
  # character(len=*), dimension(:),     allocatable, intent(out) :: minor_gases_atm_red
  # integer,          dimension(:,:),   allocatable, intent(out) :: minor_limits_gpt_atm_red
  # logical(wl),      dimension(:),     allocatable, intent(out) :: minor_scales_with_density_atm_red
  # character(len=*), dimension(:),     allocatable, intent(out) :: scaling_gas_atm_red
  # logical(wl),      dimension(:),     allocatable, intent(out) :: scale_by_complement_atm_red
  # integer,          dimension(:),     allocatable, intent(out) :: kminor_start_atm_red

  # # Local variables
  # integer :: i, j
  # integer :: idx_mnr, nm, tot_g, red_nm
  # integer :: icnt, n_elim, ng
  # logical, dimension(:), allocatable :: gas_is_present
  FT = eltype(kminor_atm)

  nm = length(minor_gases_atm)
  tot_g=0
  gas_is_present = Vector{Bool}(undef, nm)
  for i = 1:length(minor_gases_atm)
    idx_mnr = loc_in_array(minor_gases_atm[i], identifier_minor)
    gas_is_present[i] = gas_minor[idx_mnr] in available_gases.gas_name
    if gas_is_present[i]
      tot_g = tot_g + (minor_limits_gpt_atm[2,i]-minor_limits_gpt_atm[1,i]+1)
    end
  end
  red_nm = count(gas_is_present)

  if red_nm == nm
    kminor_atm_red = kminor_atm
    minor_gases_atm_red = minor_gases_atm
    minor_limits_gpt_atm_red = minor_limits_gpt_atm
    minor_scales_with_density_atm_red = minor_scales_with_density_atm
    scaling_gas_atm_red = scaling_gas_atm
    scale_by_complement_atm_red = scale_by_complement_atm
    kminor_start_atm_red = kminor_start_atm
  else
    minor_gases_atm_red= pack(minor_gases_atm, gas_is_present)
    minor_scales_with_density_atm_red = pack(minor_scales_with_density_atm,
      gas_is_present)
    scaling_gas_atm_red = pack(scaling_gas_atm,
      gas_is_present)
    scale_by_complement_atm_red = pack(scale_by_complement_atm,
      gas_is_present)
    kminor_start_atm_red = pack(kminor_start_atm,
      gas_is_present)

    minor_limits_gpt_atm_red = Array{Int}(undef, 2, red_nm)
    kminor_atm_red = zeros(FT, tot_g, size(kminor_atm,2), size(kminor_atm,3))

    icnt = 0
    n_elim = 0
    for i = 1:nm
      ng = minor_limits_gpt_atm[2,i]-minor_limits_gpt_atm[1,i]+1
      if gas_is_present[i]
        icnt = icnt + 1
        minor_limits_gpt_atm_red[1:2,icnt] = minor_limits_gpt_atm[1:2,i]
        kminor_start_atm_red[icnt] = kminor_start_atm[i]-n_elim
        for j = 1:ng
          kminor_atm_red[kminor_start_atm_red[icnt]+j-1,:,:] = kminor_atm[kminor_start_atm[i]+j-1,:,:]
        end
      else
        n_elim = n_elim + ng
      end
    end
  end
  return kminor_atm_red,
         minor_gases_atm_red,
         minor_limits_gpt_atm_red,
         minor_scales_with_density_atm_red,
         scaling_gas_atm_red,
         scale_by_complement_atm_red,
         kminor_start_atm_red

end

# ---------------------------------------------------------------------------------------
# returns flavor index; -1 if not found
function key_species_pair2flavor(flavor, key_species_pair)
  # integer :: key_species_pair2flavor
  # integer, dimension(:,:), intent(in) :: flavor
  # integer, dimension(2), intent(in) :: key_species_pair
  # integer :: iflav
  for iflav=1:size(flavor,2)
    if all(key_species_pair .== flavor[:,iflav])
      return iflav
    end
  end
  return -1
end

# ---------------------------------------------------------------------------------------
#
# create gpoint_flavor list
#   a map pointing from each g-point to the corresponding entry in the "flavor list"
#
function create_gpoint_flavor(key_species, gpt2band, flavor)
  # integer, dimension(:,:,:), intent(in) :: key_species
  # integer, dimension(:), intent(in) :: gpt2band
  # integer, dimension(:,:), intent(in) :: flavor
  # integer, dimension(:,:), intent(out), allocatable :: gpoint_flavor
  # integer :: ngpt, igpt, iatm
  ngpt = length(gpt2band)
  gpoint_flavor = Array{Int}(undef, 2,ngpt)
  for igpt=1:ngpt
    for iatm=1:2
      gpoint_flavor[iatm,igpt] = key_species_pair2flavor(
        flavor,
        rewrite_key_species_pair(key_species[:,iatm,gpt2band[igpt]])
      )
    end
  end
  return gpoint_flavor
end

#--------------------------------------------------------------------------------------------------------------------
#
# Utility function to combine optical depths from gas absorption and Rayleigh scattering
#   (and reorder them for convenience, while we're at it)
#
function combine_and_reorder!(tau, tau_rayleigh, has_rayleigh, optical_props::ty_optical_props_arry{FT}) where FT
  # real(FT), dimension(:,:,:),   intent(in) :: tau
  # real(FT), dimension(:,:,:),   intent(in) :: tau_rayleigh
  # logical,                      intent(in) :: has_rayleigh
  # class(ty_optical_props_arry), intent(inout) :: optical_props

  # integer :: ncol, nlay, ngpt, nmom

  ncol = size(tau, 3)
  nlay = size(tau, 2)
  ngpt = size(tau, 1)
  if !has_rayleigh
    # index reorder (ngpt, nlay, ncol) -> (ncol,nlay,gpt)
    permutedims!(optical_props.tau, tau, [3,2,1])
      if optical_props isa ty_optical_props_2str
        optical_props.ssa .= FT(0)
        optical_props.g   .= FT(0)
      end
  else
    # combine optical depth and rayleigh scattering
      if optical_props isa ty_optical_props_1scl
        # User is asking for absorption optical depth
        permutedims!(optical_props.tau, tau, [3,2,1])

      elseif optical_props isa ty_optical_props_2str
        combine_and_reorder_2str!(optical_props, ncol, nlay, ngpt,       tau, tau_rayleigh)
      end
  end
end

#--------------------------------------------------------------------------------------------------------------------
# Sizes of tables: pressure, temperate, eta (mixing fraction)
#   Equivalent routines for the number of gases and flavors (get_ngas(), get_nflav()) are defined above because they're
#   used in function defintions
# Table kmajor has dimensions (ngpt, neta, npres, ntemp)
#--------------------------------------------------------------------------------------------------------------------

# return extent of eta dimension
get_neta(this::ty_gas_optics_rrtmgp) = size(this.kmajor,2)

# return the number of pressures in reference profile
#   absorption coefficient table is one bigger since a pressure is repeated in upper/lower atmos
get_npres(this::ty_gas_optics_rrtmgp) = size(this.kmajor,3)-1

# return the number of temperatures
get_ntemp(this::ty_gas_optics_rrtmgp) = size(this.kmajor,4)

# return the number of temperatures for Planck function
get_nPlanckTemp(this::ty_gas_optics_rrtmgp) = size(this.totplnk,1) # dimensions are Planck-temperature, band

#### Generic procedures for checking sizes, limits

# Extents
check_extent(array, s, label) = @assert all(size(array).==s)

# Values
check_range(val, minV, maxV, label) = any(val .< minV) || any(val .> maxV) ? trim(label) * " values out of range." : ""

include("mo_gas_optics_kernels.jl")

end #module
