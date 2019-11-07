"""
    mo_fluxes

Compute output quantities from RTE based on spectrally-resolved flux profiles
  - This module contains an abstract class and a broadband implementation that sums over all spectral points
  - The abstract class defines the routines that extensions must implement: reduce() and are_desired()
  - The intent is for users to extend it as required, using mo_flxues_broadband as an example
"""
module mo_fluxes

using ..fortran_intrinsics
using ..mo_optical_props
using ..mo_fluxes_broadband_kernels

export ty_fluxes_broadband, are_desired, reduce!

"""
    ty_fluxes{FT}

Fluxes struct
  reduce! function accepts spectral flux profiles, computes desired outputs
  are_desired() returns a logical - does it makes sense to invoke reduce()?
"""
abstract type ty_fluxes{FT} end

"""
    ty_fluxes_broadband{FT} <: ty_fluxes{FT}

Class implementing broadband integration for the complete flux profile
"""
mutable struct ty_fluxes_broadband{FT} <: ty_fluxes{FT}
  flux_up#::Array{FT,2}
  flux_dn#::Array{FT,2}
  flux_net#::Array{FT,2}
  flux_dn_dir#::Array{FT,2}
end

ty_fluxes_broadband(FT) = ty_fluxes_broadband{FT}(ntuple(i->nothing, 4)...)

"""
    reduce!(this::ty_fluxes_broadband,
                 gpt_flux_up::A,
                 gpt_flux_dn::A,
                 spectral_disc::ty_optical_props,
                 top_at_1::Bool,
                 gpt_flux_dn_dir=nothing)

Sum over the spectral dimension and return the whole profile

 - `this` a `ty_fluxes_broadband` [inout]
 - `gpt_flux_up` a `Array{FT,3}(ncol, nlay+1, ngpt)` [in] Fluxes by gpoint [W/m2]
 - `gpt_flux_dn` a `Array{FT,3}(ncol, nlay+1, ngpt)` [in] Fluxes by gpoint [W/m2]
 - `spectral_disc` a  `ty_optical_props` [in] derived type with spectral information
 - `top_at_1` bool indicating at top
optional:
 - `gpt_flux_dn_dir` a `Array{FT,3}` [in] Direct flux down
"""
function reduce!(this::ty_fluxes_broadband,
                 gpt_flux_up::A,
                 gpt_flux_dn::A,
                 spectral_disc::ty_optical_props,
                 top_at_1::Bool,
                 gpt_flux_dn_dir=nothing) where A

  ncol,nlev,ngpt = size(gpt_flux_up)

  # Check array sizes
  @assert all(size(gpt_flux_dn) .== (ncol,nlev,ngpt))
  gpt_flux_dn_dir ≠ nothing &&  @assert all(size(gpt_flux_dn_dir) .== (ncol,nlev,ngpt))
  associated(this.flux_up) && @assert all(size(this.flux_up) .== (ncol,nlev))
  associated(this.flux_dn) && @assert all(size(this.flux_dn) .== (ncol,nlev))
  associated(this.flux_net) && @assert all(size(this.flux_net) .== (ncol,nlev))
  associated(this.flux_dn_dir) && @assert all(size(this.flux_dn_dir) .== (ncol,nlev))

  # Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
  @assert !(associated(this.flux_dn_dir) && !present(gpt_flux_dn_dir))

  # Broadband fluxes - call the kernels
  if associated(this.flux_up    )
    this.flux_up .= sum_broadband(ncol, nlev, ngpt, gpt_flux_up)
  end
  if associated(this.flux_dn    )
    this.flux_dn .= sum_broadband(ncol, nlev, ngpt, gpt_flux_dn)
  end
  if associated(this.flux_dn_dir)
    this.flux_dn_dir = sum_broadband(ncol, nlev, ngpt, gpt_flux_dn_dir)
  end

  if associated(this.flux_net)

    #  Reuse down and up results if possible
    if associated(this.flux_dn) && associated(this.flux_up)
      this.flux_net = net_broadband(ncol, nlev,      this.flux_dn, this.flux_up)
    else
      this.flux_net = net_broadband(ncol, nlev, ngpt, gpt_flux_dn,  gpt_flux_up)
    end
  end
end

"""
    are_desired(this::ty_fluxes_broadband)

Bool indicating if any fluxes desired from this set of
g-point fluxes. We can tell because memory will be
allocated for output
"""
are_desired(this::ty_fluxes_broadband) =
  any( [associated(this.flux_up),
        associated(this.flux_dn),
        associated(this.flux_dn_dir),
        associated(this.flux_net)] )

end
