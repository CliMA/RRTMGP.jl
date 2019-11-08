"""
    mo_rte_sw

Contains a single routine to compute direct and diffuse fluxes of solar radiation given
 - atmospheric optical properties on a spectral grid
 - information about vertical ordering
 - boundary conditions
 - solar zenith angle
 - spectrally-resolved incident colimated flux
 - surface albedos for direct and diffuse radiation
 - (optionally) a boundary condition for incident diffuse radiation

It is the user's responsibility to ensure that
boundary conditions (incident fluxes, surface
albedos) are on the same spectral grid as the
optical properties.

Final output is via user-extensible `ty_fluxes`
which must reduce the detailed spectral fluxes to
whatever summary the user needs.

The routine does error checking and choses which
lower-level kernel to invoke based on what kinds
of optical properties are supplied
"""
module mo_rte_sw

using ..mo_util_array
using ..mo_optical_props
using ..mo_fluxes
using ..mo_rte_solver_kernels
using ..fortran_intrinsics

export rte_sw!, expand_and_transpose

"""
    rte_sw!(atmos::ty_optical_props_arry,
                 top_at_1,
                 mu0::Array{FT},
                 inc_flux,
                 sfc_alb_dir,
                 sfc_alb_dif,
                 fluxes,
                 inc_flux_dif=nothing)

Compute fluxes `fluxes` (`ty_fluxes`), given

 - `atmos` Optical properties provided as arrays (`ty_optical_props_arry`)
 - `top_at_1` bool indicating if the top of the domain is at index 1
             (if not, ordering is bottom-to-top)
 - `mu0` cosine of solar zenith angle (ncol)
 - `inc_flux` incident flux at top of domain [W/m2] (ncol, ngpt)
 - `sfc_alb_dir` surface albedo for direct and
 - `sfc_alb_dif` diffuse radiation (nband, ncol)
and, optionally,
 - `inc_flux_dif` incident diffuse flux at top of domain [W/m2] (ncol, ngpt)
"""
function rte_sw!(atmos::ty_optical_props_arry,
                 top_at_1,
                 mu0::Array{FT},
                 inc_flux,
                 sfc_alb_dir,
                 sfc_alb_dif,
                 fluxes,
                 inc_flux_dif=nothing) where {FT<:AbstractFloat}

  ncol  = get_ncol(atmos)
  nlay  = get_nlay(atmos)
  ngpt  = get_ngpt(atmos)
  nband = get_nband(atmos)

  # Error checking -- consistency / sizes / validity of values
  @assert are_desired(fluxes)
  @assert size(mu0, 1) ==  ncol
  @assert !any_vals_outside(mu0, FT(0), FT(1))
  @assert all(size(inc_flux) .== (ncol, ngpt))
  @assert !any_vals_less_than(inc_flux, FT(0))
  if present(inc_flux_dif)
    @assert all(size(inc_flux_dif) .== (ncol, ngpt))
    @assert !any_vals_less_than(inc_flux_dif, FT(0))
  end
  @assert all(size(sfc_alb_dir) .== (nband, ncol))
  @assert !any_vals_outside(sfc_alb_dir,  FT(0), FT(1))
  @assert all(size(sfc_alb_dif) .== (nband, ncol))
  @assert !any_vals_outside(sfc_alb_dif, FT(0), FT(1))
  @assert !(atmos isa ty_optical_props_nstr) # not yet implemented...

  gpt_flux_up  = Array{FT}(undef,ncol, nlay+1, ngpt)
  gpt_flux_dn  = Array{FT}(undef,ncol, nlay+1, ngpt)
  gpt_flux_dir = Array{FT}(undef,ncol, nlay+1, ngpt)

  # Lower boundary condition -- expand surface albedos by band to gpoints
  #   and switch dimension ordering
  sfc_alb_dir_gpt = expand_and_transpose(atmos, sfc_alb_dir)
  sfc_alb_dif_gpt = expand_and_transpose(atmos, sfc_alb_dif)

  # Compute the radiative transfer...
  #
  # Apply boundary conditions
  #   On input flux_dn is the diffuse component; the last action in each solver is to add
  #   direct and diffuse to represent the total, consistent with the LW
  #
  gpt_flux_dir = apply_BC(ncol, nlay, ngpt, top_at_1,   inc_flux, mu0)

  if present(inc_flux_dif)
    gpt_flux_dn  = apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dif)
  else
    gpt_flux_dn  = apply_BC(ncol, nlay, ngpt, top_at_1, FT)
  end

  if atmos isa ty_optical_props_1scl
    # Direct beam only
    validate!(atmos)
    sw_solver_noscat!(ncol, nlay, ngpt, top_at_1,
                          atmos.tau, mu0,
                          gpt_flux_dir)
    # No diffuse flux
    # gpt_flux_up .= FT(0)
    # gpt_flux_dn .= FT(0)
  elseif atmos isa ty_optical_props_2str
    # two-stream calculation with scattering
    validate!(atmos)

    sw_solver_2stream!(ncol, nlay, ngpt, top_at_1,
                           atmos.tau, atmos.ssa, atmos.g, mu0,
                           sfc_alb_dir_gpt, sfc_alb_dif_gpt,
                           gpt_flux_up, gpt_flux_dn, gpt_flux_dir)
  end

  # ...and reduce spectral fluxes to desired output quantities
  reduce!(fluxes,gpt_flux_up, gpt_flux_dn, atmos, top_at_1, gpt_flux_dir)
end

"""
    expand_and_transpose(ops::ty_optical_props,arr_in::Array{FT}) where FT

Expand from band to g-point dimension, transpose
dimensions (nband, ncol) -> (ncol,ngpt), of `arr_out`, given

 - `ops` - a `ty_optical_props`
 - `arr_in` - input array
"""
function expand_and_transpose(ops::ty_optical_props,arr_in::Array{FT}) where FT
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

end #module