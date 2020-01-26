# Radiative Transfer Equation (RTE)

"""
    AbstractRTE{FT<:AbstractFloat, I<:Int}

A Radiative Transfer Equation
"""
abstract type AbstractRTE{FT<:AbstractFloat, I<:Int} end

"""
    RTEBase{FT<:AbstractFloat,I<:Int}

Common fields in the RTE including
 - shortwave
 - longwave
 - with scattering
 - without scattering

# Fields
$(DocStringExtensions.FIELDS)
"""
struct RTEBase{FT<:AbstractFloat,I<:Int}
  "Broadband fluxes"
  fluxes::FluxesBroadBand{FT}
  mo::MeshOrientation{I}
  bcs::AbstractRadiativeBoundaryConditions{FT}
  "Upward flux ``[W/m^2]``"
  gpt_flux_up::Array{FT,3}
  "Downward flux ``[W/m^2]`` Top level must contain incident flux boundary condition"
  gpt_flux_dn::Array{FT,3}
 # - `sfc_emis` - surface emissivity
 # - `flux_up` upward radiance [W/m2-str]
 # - `flux_dn` downward radiance, Top level must contain incident flux boundary condition

  gpt_flux_dir::Union{Array{FT,3},Nothing}
  function RTEBase(fluxes::FluxesBroadBand{FT},
                   mo::MeshOrientation{I},
                   bcs::AbstractRadiativeBoundaryConditions{FT},
                   op::AbstractOpticalPropsArry{FT}) where {FT,I}
    validate!(op)
    ncol  = get_ncol(op)
    nlay  = get_nlay(op)
    ngpt  = get_ngpt(op)
    gpt_flux_up  = zeros(FT, ncol, nlay+1, ngpt)
    gpt_flux_dn  = zeros(FT, ncol, nlay+1, ngpt)
    gpt_flux_dir  = fluxes.flux_dn_dir == nothing ? nothing : similar(gpt_flux_dn)
    return new{FT,I}(fluxes, mo, bcs, gpt_flux_up, gpt_flux_dn, gpt_flux_dir)
  end
end

struct RTELongWaveNoScattering{FT,I} <: AbstractRTE{FT,I}
  base::RTEBase{FT,I}
  sources::SourceFuncLongWave{FT, I}
  angle_disc::Union{GaussQuadrature{FT,I},Nothing}
  sfc_emis_gpt::Array{FT,2}
end

struct RTELongWave{FT,I} <: AbstractRTE{FT,I}
  base::RTEBase{FT,I}
  sources::SourceFuncLongWave{FT, I}
  angle_disc::Union{GaussQuadrature{FT,I},Nothing}
  sfc_emis_gpt::Array{FT,2}
end

struct RTEShortWaveNoScattering{FT,I} <: AbstractRTE{FT,I}
  base::RTEBase{FT,I}
  μ_0::Array{FT}
end

struct RTEShortWave{FT,I} <: AbstractRTE{FT,I}
  base::RTEBase{FT,I}
  μ_0::Array{FT}
  sfc_alb_dir_gpt::Array{FT}
  sfc_alb_dif_gpt::Array{FT}
end

function RTELongWave(base::RTEBase{FT,I},
                     sources::SourceFuncLongWave{FT, I},
                     angle_disc::Union{GaussQuadrature{FT,I},Nothing},
                     op::OneScalar{FT,I}) where {FT,I}
  ncol  = get_ncol(op)
  nlay  = get_nlay(op)
  ngpt  = get_ngpt(op)
  @assert get_ncol(sources) == ncol
  @assert get_nlay(sources) == nlay
  @assert get_ngpt(sources) == ngpt
  angle_disc == nothing && (angle_disc = GaussQuadrature(FT, 1))
  # Lower boundary condition -- expand surface emissivity by band to gpoints
  sfc_emis_gpt = expand_and_transpose(op, base.bcs.sfc_emis)
  return RTELongWaveNoScattering{FT,I}(base, sources, angle_disc, sfc_emis_gpt)
end

function RTELongWave(base::RTEBase{FT,I},
                     sources::SourceFuncLongWave{FT, I},
                     angle_disc::Union{GaussQuadrature{FT,I},Nothing},
                     op::TwoStream{FT,I}) where {FT,I}
  ncol  = get_ncol(op)
  nlay  = get_nlay(op)
  ngpt  = get_ngpt(op)
  @assert get_ncol(sources) == ncol
  @assert get_nlay(sources) == nlay
  @assert get_ngpt(sources) == ngpt
  angle_disc == nothing && (angle_disc = GaussQuadrature(FT, 1))
  # Lower boundary condition -- expand surface emissivity by band to gpoints
  sfc_emis_gpt = expand_and_transpose(op, base.bcs.sfc_emis)
  return RTELongWave{FT,I}(base, sources, angle_disc, sfc_emis_gpt)
end

function RTEShortWave(base::RTEBase{FT,I},
                      μ_0::Array{FT},
                      op::OneScalar{FT}) where {FT,I}
  # Error checking -- consistency / sizes / validity of values
  ncol  = get_ncol(op)
  @assert size(μ_0, 1) ==  ncol
  @assert !any_vals_outside(μ_0, FT(0), FT(1))
  return RTEShortWaveNoScattering{FT,I}(base, μ_0)
end

function RTEShortWave(base::RTEBase{FT,I},
                      μ_0::Array{FT},
                      op::TwoStream{FT}) where {FT,I}
  # Lower boundary condition -- expand surface albedos by band to gpoints
  #   and switch dimension ordering
  sfc_alb_dir_gpt = expand_and_transpose(op, base.bcs.sfc_alb_direct)
  sfc_alb_dif_gpt = expand_and_transpose(op, base.bcs.sfc_alb_diffuse)
  return RTEShortWave{FT,I}(base, μ_0, sfc_alb_dir_gpt, sfc_alb_dif_gpt)
end

