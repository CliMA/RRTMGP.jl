"""
    RTESolver

# Shortwave

Contains a single routine to compute direct and diffuse fluxes of solar radiation given
 - atmospheric optical properties on a spectral grid
 - information about vertical ordering
 - boundary conditions
 - solar zenith angle
 - spectrally-resolved incident collimated flux
 - surface albedos for direct and diffuse radiation
 - (optionally) a boundary condition for incident diffuse radiation

It is the user's responsibility to ensure that
boundary conditions (incident fluxes, surface
albedos) are on the same spectral grid as the
optical properties.

Final output is via user-extensible `AbstractFluxes`
which must reduce the detailed spectral fluxes to
whatever summary the user needs.

The routine does error checking and choses which
lower-level kernel to invoke based on what kinds
of optical properties are supplied

# Longwave

Contains a single routine to compute direct and diffuse fluxes of solar radiation given
 - atmospheric optical properties, spectrally-resolved
 - information about vertical ordering
 - internal Planck source functions, defined per g-point on the same spectral grid at the atmosphere
 - boundary conditions: surface emissivity defined per band
 - optionally, a boundary condition for incident diffuse radiation
 - optionally, an integer number of angles at which to do Gaussian quadrature if scattering is neglected

If optical properties are supplied via a `OneScalar` (absorption optical thickenss only)
  then an emission/absorption solver is called
  If optical properties are supplied via a `TwoStream` fluxes are computed via
  two-stream calculations and adding.

It is the user's responsibility to ensure that emissivity is on the same
 spectral grid as the optical properties.

Final output is via user-extensible AbstractFluxes which must reduce the detailed spectral fluxes to
 whatever summary the user needs.

The routine does error checking and choses which lower-level kernel to invoke based on
 what kinds of optical properties are supplied

"""
module RTESolver

using ..RadiativeBoundaryConditions
using ..AngularDiscretizations
using ..Utilities
using ..OpticalProps
using ..Fluxes
using ..SourceFunctions
using ..FortranIntrinsics

export rte_lw!, rte_sw!, expand_and_transpose

"""
    rte_sw!(atmos::AbstractOpticalPropsArry{FT},
            top_at_1::Bool,
            μ_0::Array{FT},
            bcs::ShortwaveBCs{FT},
            fluxes::FluxesBroadBand{FT}) where {FT<:AbstractFloat}

Computes

 - `fluxes` broadband fluxes, see [`AbstractFluxes`](@ref)

given

 - `atmos` Optical properties provided as arrays (`AbstractOpticalPropsArry`)
 - `top_at_1` indicates that the top of the domain is at index 1
 - `μ_0` cosine of solar zenith angle (ncol)
 - `bcs` boundary conditions, see [`ShortwaveBCs`](@ref)
"""
function rte_sw!(atmos::AbstractOpticalPropsArry{FT},
                 top_at_1::Bool,
                 μ_0::Array{FT},
                 bcs::ShortwaveBCs{FT},
                 fluxes::FluxesBroadBand{FT}) where {FT<:AbstractFloat}

  ncol  = get_ncol(atmos)
  nlay  = get_nlay(atmos)
  ngpt  = get_ngpt(atmos)

  # Error checking -- consistency / sizes / validity of values
  @assert size(μ_0, 1) ==  ncol
  @assert !any_vals_outside(μ_0, FT(0), FT(1))

  gpt_flux_up  = zeros(FT, ncol, nlay+1, ngpt)
  gpt_flux_dn  = zeros(FT, ncol, nlay+1, ngpt)
  gpt_flux_dir = zeros(FT, ncol, nlay+1, ngpt)

  # Lower boundary condition -- expand surface albedos by band to gpoints
  #   and switch dimension ordering
  sfc_alb_dir_gpt = expand_and_transpose(atmos, bcs.sfc_alb_direct)
  sfc_alb_dif_gpt = expand_and_transpose(atmos, bcs.sfc_alb_diffuse)

  # Compute the radiative transfer...
  #
  # Apply boundary conditions
  #   On input flux_dn is the diffuse component; the last action in each solver is to add
  #   direct and diffuse to represent the total, consistent with the LW
  #
  apply_BC!(gpt_flux_dir, ncol, nlay, ngpt, top_at_1,   bcs.toa_flux, μ_0)

  if present(bcs.inc_flux_diffuse)
    apply_BC!(gpt_flux_dn, ncol, nlay, ngpt, top_at_1, bcs.inc_flux_diffuse)
  else
    apply_BC!(gpt_flux_dn, ncol, nlay, ngpt, top_at_1)
  end

  validate!(atmos)
  if atmos isa OneScalar
    # Direct beam only
    sw_solver_noscat!(ncol, nlay, ngpt, top_at_1,
                      atmos.τ, μ_0,
                      gpt_flux_dir)
    # No diffuse flux
    # gpt_flux_up .= FT(0)
    # gpt_flux_dn .= FT(0)
  elseif atmos isa TwoStream
    # two-stream calculation with scattering
    sw_solver!(ncol, nlay, ngpt, top_at_1,
               atmos,
               μ_0,
               sfc_alb_dir_gpt,
               sfc_alb_dif_gpt,
               gpt_flux_up,
               gpt_flux_dn,
               gpt_flux_dir)
  end

  # ...and reduce spectral fluxes to desired output quantities
  reduce!(fluxes,
          gpt_flux_up,
          gpt_flux_dn,
          atmos,
          top_at_1,
          gpt_flux_dir)
end

"""
    rte_lw!(optical_props::AbstractOpticalPropsArry{FT},
            top_at_1::Bool,
            sources::SourceFuncLongWave{FT, I},
            bcs::LongwaveBCs{FT},
            fluxes::FluxesBroadBand{FT},
            n_gauss_angles=nothing) where {FT<:AbstractFloat,I<:Int}

Interface using only optical properties and source functions as inputs; fluxes as outputs.

 - `optical_props` optical properties
 - `top_at_1` indicates that the top of the domain is at index 1
 - `sources` radiation sources
 - `bcs` boundary conditions, see [`LongwaveBCs`](@ref)
 - `fluxes` broadband fluxes, see [`AbstractFluxes`](@ref)
 - `n_gauss_angles` Number of angles used in Gaussian quadrature (no-scattering solution)
"""
function rte_lw!(optical_props::AbstractOpticalPropsArry{FT},
                 top_at_1::Bool,
                 sources::SourceFuncLongWave{FT, I},
                 bcs::LongwaveBCs{FT},
                 fluxes::FluxesBroadBand{FT},
                 angle_disc::Union{GaussQuadrature{FT,I},Nothing}=nothing) where {FT<:AbstractFloat,I<:Int}

  # Error checking
  ncol  = get_ncol(optical_props)
  nlay  = get_nlay(optical_props)
  ngpt  = get_ngpt(optical_props)

  # Error checking -- consistency of sizes and validity of values
  @assert get_ncol(sources) == ncol
  @assert get_nlay(sources) == nlay
  @assert get_ngpt(sources) == ngpt

  if angle_disc == nothing
    angle_disc = GaussQuadrature(FT, 1)
  end

  # Ensure values of τ, ssa, and g are reasonable
  validate!(optical_props)

  # Lower boundary condition -- expand surface emissivity by band to gpoints
  gpt_flux_up  = Array{FT}(undef,ncol,nlay+1,ngpt)
  gpt_flux_dn  = Array{FT}(undef,ncol,nlay+1,ngpt)

  sfc_emis_gpt = expand_and_transpose(optical_props, bcs.sfc_emis)

  #   Upper boundary condition
  if bcs.inc_flux ≠ nothing
    apply_BC!(gpt_flux_dn, ncol, nlay, ngpt, top_at_1, bcs.inc_flux)
  else
    # Default is zero incident diffuse flux
    apply_BC!(gpt_flux_dn, ncol, nlay, ngpt, top_at_1)
  end

  # Compute the radiative transfer...
  if optical_props isa OneScalar

    # No scattering two-stream calculation
    lw_solver_noscat_GaussQuad!(ncol, nlay, ngpt, top_at_1,
                                angle_disc,
                                optical_props.τ,
                                sources,
                                sfc_emis_gpt,
                                gpt_flux_up,
                                gpt_flux_dn)

  elseif optical_props isa TwoStream
    # two-stream calculation with scattering
    lw_solver!(ncol, nlay, ngpt, top_at_1,
               optical_props,
               sources,
               sfc_emis_gpt,
               gpt_flux_up,
               gpt_flux_dn)
  end

  reduce!(fluxes, gpt_flux_up, gpt_flux_dn, optical_props, top_at_1)
end

"""
    expand_and_transpose(ops::AbstractOpticalProps,arr_in::Array{FT}) where FT

Expand from band to g-point dimension, transpose
dimensions (nband, ncol) -> (ncol,ngpt), of `arr_out`, given

 - `ops` - a `AbstractOpticalProps`
 - `arr_in` - input array
"""
function expand_and_transpose(ops::AbstractOpticalProps, arr_in::Array{FT}) where FT
  ncol  = size(arr_in,2)
  nband = get_nband(ops)
  ngpt  = get_ngpt(ops)
  arr_out = Array{FT}(undef,ncol, ngpt)
  limits = get_band_lims_gpoint(ops)
  for iband = 1:nband
    for icol = 1:ncol
      for igpt = limits[1, iband]:limits[2, iband]
        arr_out[icol, igpt] = arr_in[iband,icol]
      end
    end
  end
  return arr_out
end

include(joinpath("kernels","RTESolver_kernels.jl"))

end #module