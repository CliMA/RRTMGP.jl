"""
    mo_fluxes_byband

This module is for packaging output quantities from RRTMGP based on
spectral flux profiles. This implementation provides band-by-band
flux profiles.
"""
module mo_fluxes_byband

using DocStringExtensions
using ..mo_fluxes
using ..mo_optical_props
using ..mo_fluxes_byband_kernels

"""
    ty_fluxes_byband{FT} <: ty_fluxes{FT}

Contains a `fluxes_broadband` and fluxes by-band, including

 - `bnd_flux_up`     upward flux
 - `bnd_flux_dn`     downward
 - `bnd_flux_net`    net
 - `bnd_flux_dn_dir` downward direct

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct ty_fluxes_byband{FT} <: ty_fluxes{FT}
  fluxes_broadband
  bnd_flux_up#::Array{FT,3}
  bnd_flux_dn#::Array{FT,3}
  bnd_flux_net#::Array{FT,3}
  bnd_flux_dn_dir#::Array{FT,3}
end

"""
    reduce_byband(this::ty_fluxes_byband,
                       gpt_flux_up::Array{FT,3},
                       gpt_flux_dn::Array{FT,3},
                       spectral_disc::ty_optical_props,
                       top_at_1::Bool,
                       gpt_flux_dn_dir::Union{Nothing,Array{FT,3}}=nothing)

Reduces fluxes by-band to broadband in `ty_fluxes_byband` `this`, given

 - `gpt_flux_up` fluxes by gpoint [W/m2]
 - `gpt_flux_dn` fluxes by gpoint [W/m2]
 - `spectral_disc` a `ty_optical_props` struct containing spectral information
 - `top_at_1` bool indicating at top
and, optionally,
 - `gpt_flux_dn_dir` direct flux downward
"""
function reduce_byband(this::ty_fluxes_byband,
                       gpt_flux_up::Array{FT,3},
                       gpt_flux_dn::Array{FT,3},
                       spectral_disc::ty_optical_props,
                       top_at_1::Bool,
                       gpt_flux_dn_dir::Union{Nothing,Array{FT,3}}=nothing) where {FT<:AbstractFloat}
  ncol, nlev = size(gpt_flux_up)
  ngpt = get_ngpt(spectral_disc)
  nbnd = get_nband(spectral_disc)
  band_lims = deepcopy(get_band_lims_gpoint(spectral_disc))

  # Compute broadband fluxes
  #   This also checks that input arrays are consistently sized
  #
  reduce_broadband!(this.fluxes_broadband, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir)

  # Check sizes
  @assert size(gpt_flux_up, 3) == ngpt
  associated(this.bnd_flux_up)     && @assert all(size(this.bnd_flux_up) .== (ncol, nlev, nbnd))
  associated(this.bnd_flux_dn)     && @assert all(size(this.bnd_flux_dn) .== (ncol, nlev, nbnd))
  associated(this.bnd_flux_dn_dir) && @assert all(size(this.bnd_flux_dn_dir) .== (ncol, nlev, nbnd))
  associated(this.bnd_flux_net)    && @assert all(size(this.bnd_flux_net) .== (ncol, nlev, nbnd))
  # Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
  @assert !(associated(this.bnd_flux_dn_dir) && !present(gpt_flux_dn_dir))

  # Band-by-band fluxes

  # Up flux
  if associated(this.bnd_flux_up)
    this.bnd_flux_up = sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_up)
  end

  # Down flux
  if associated(this.bnd_flux_dn)
    this.bnd_flux_dn = sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn)
  end

  # Direct Down flux
  if associated(this.bnd_flux_dn_dir)
    this.bnd_flux_dn_dir = sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn_dir)
  end

  # Net flux
  if(associated(this.bnd_flux_net))
    #
    #  Reuse down and up results if possible
    #
    if(associated(this.bnd_flux_dn) && associated(this.bnd_flux_up))
      this.bnd_flux_net = net_byband(ncol, nlev,       nbnd, this.bnd_flux_dn, this.bnd_flux_up)
    else
      this.bnd_flux_net = net_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn, gpt_flux_up)
    end
  end
end

"""
    are_desired_byband(this::ty_fluxes_byband)

Boolean indicating if any fluxes desired from this
set of g-point fluxes.
"""
are_desired_byband(this::ty_fluxes_byband) =
  any([associated(this.bnd_flux_up),
       associated(this.bnd_flux_dn),
       associated(this.bnd_flux_dn_dir),
       associated(this.bnd_flux_net),
       are_desired(this.fluxes_broadband)])

end
