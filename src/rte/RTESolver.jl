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

If optical properties are supplied via a `OneScalar` (absorption optical thickness only)
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

using DocStringExtensions
using ..RadiativeBoundaryConditions
using ..MeshOrientations
using ..AngularDiscretizations
using ..Utilities
using ..OpticalProps
using ..Fluxes
using ..SourceFunctions
using ..FortranIntrinsics

import ..Fluxes: reduce!
export rte_lw!, rte_sw!, expand_and_transpose

include("RTE.jl")

"""
    rte_sw!(fluxes::FluxesBroadBand{FT},
            optical_props::AbstractOpticalPropsArry{FT},
            mesh_orientation::MeshOrientation{I},
            bcs::ShortwaveBCs{FT},
            μ_0::Array{FT}) where {FT<:AbstractFloat,I<:Int}

Compute broadband radiative fluxes

 - `fluxes` broadband fluxes, see [`FluxesBroadBand`](@ref)

given

 - `optical_props` optical properties, see [`AbstractOpticalPropsArry`](@ref)
 - `mesh_orientation` mesh orientation, see [`MeshOrientation`](@ref)
 - `bcs` boundary conditions, see [`ShortwaveBCs`](@ref)
 - `μ_0` cosine of solar zenith angle (ncol)
"""
function rte_sw!(fluxes::FluxesBroadBand{FT},
                 optical_props::AbstractOpticalPropsArry{FT,I},
                 mesh_orientation::MeshOrientation{I},
                 bcs::ShortwaveBCs{FT},
                 μ_0::Array{FT}) where {FT<:AbstractFloat,I<:Int}
  base = RTEBase(fluxes, mesh_orientation, bcs, optical_props)
  rte = RTEShortWave(base, μ_0, optical_props)

  # Compute the radiative transfer...
  solve!(rte, optical_props, mesh_orientation)

  # ...and reduce spectral fluxes to desired output quantities
  reduce!(rte.base, optical_props)
  fluxes = rte.base.fluxes
  return nothing
end

"""
    rte_lw!(fluxes::FluxesBroadBand{FT},
            optical_props::AbstractOpticalPropsArry{FT},
            mesh_orientation::MeshOrientation{I},
            bcs::LongwaveBCs{FT},
            sources::SourceFuncLongWave{FT, I},
            angle_disc::Union{GaussQuadrature{FT,I},Nothing}=nothing) where {FT<:AbstractFloat,I<:Int}

Compute broadband radiative fluxes

 - `fluxes` broadband fluxes, see [`FluxesBroadBand`](@ref)

given

 - `optical_props` optical properties, see [`AbstractOpticalPropsArry`](@ref)
 - `mesh_orientation` mesh orientation, see [`MeshOrientation`](@ref)
 - `bcs` boundary conditions, see [`LongwaveBCs`](@ref)
 - `sources` radiation sources, see [`SourceFuncLongWave`](@ref)
 - `angle_disc` Gaussian quadrature for angular discretization, [`GaussQuadrature`](@ref)
"""
function rte_lw!(fluxes::FluxesBroadBand{FT},
                 optical_props::AbstractOpticalPropsArry{FT,I},
                 mesh_orientation::MeshOrientation{I},
                 bcs::LongwaveBCs{FT},
                 sources::SourceFuncLongWave{FT, I},
                 angle_disc::Union{GaussQuadrature{FT,I},Nothing}=nothing) where {FT<:AbstractFloat,I<:Int}

  base = RTEBase(fluxes, mesh_orientation, bcs, optical_props)
  rte = RTELongWave(base, sources, angle_disc, optical_props)

  # Compute the radiative transfer...
  solve!(rte, optical_props, mesh_orientation)

  reduce!(rte.base, optical_props)
  fluxes = rte.base.fluxes
  return nothing
end

"""
    reduce!(base, optical_props)

Wrapper for [`reduce!`](@ref Fluxes.reduce!)
"""
reduce!(base, optical_props) = reduce!(base.fluxes, base.gpt_flux_up, base.gpt_flux_dn, optical_props, base.gpt_flux_dir)

"""
    expand_and_transpose(ops::AbstractOpticalProps,arr_in::Array{FT}) where FT

Expand from band to g-point dimension, transpose
dimensions (`nband`, `ncol`) -> (`ncol`, `ngpt`), of `arr_out`, given

 - `ops` - optical properties, see [`AbstractOpticalProps`](@ref)
 - `arr_in` - input array
"""
function expand_and_transpose(ops::AbstractOpticalProps{FT}, arr_in::Array{FT}) where {FT<:AbstractFloat}
  ncol  = size(arr_in,2)
  nband = get_nband(ops)
  ngpt  = get_ngpt(ops)
  arr_out = Array{FT}(undef,ncol, ngpt)
  limits = get_band_lims_gpoint(ops)
  @inbounds for iband = 1:nband
    @inbounds for icol = 1:ncol
      @inbounds for igpt = limits[1, iband]:limits[2, iband]
        arr_out[icol, igpt] = arr_in[iband,icol]
      end
    end
  end
  return arr_out
end

include("RTESolver_kernels.jl")

end #module