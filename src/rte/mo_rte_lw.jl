"""
    mo_rte_lw

Contains a single routine to compute direct and diffuse fluxes of solar radiation given
 - atmospheric optical properties, spectrally-resolved
 - information about vertical ordering
 - internal Planck source functions, defined per g-point on the same spectral grid at the atmosphere
 - boundary conditions: surface emissivity defined per band
 - optionally, a boundary condition for incident diffuse radiation
 - optionally, an integer number of angles at which to do Gaussian quadrature if scattering is neglected

If optical properties are supplied via a `ty_optical_props_1scl` (absorption optical thickenss only)
  then an emission/absorption solver is called
  If optical properties are supplied via a `ty_optical_props_2str` fluxes are computed via
  two-stream calculations and adding.

It is the user's responsibility to ensure that emissivity is on the same
 spectral grid as the optical properties.

Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
 whatever summary the user needs.

The routine does error checking and choses which lower-level kernel to invoke based on
 what kinds of optical properties are supplied

"""
module mo_rte_lw

using ..mo_util_array
using ..mo_optical_props
using ..mo_fluxes
using ..mo_rte_solver_kernels
using ..fortran_intrinsics

export rte_lw!, expand_and_transpose


"""
    rte_lw!(optical_props, top_at_1,
                sources, sfc_emis,
                fluxes,
                inc_flux=nothing, n_gauss_angles=nothing)

Interface using only optical properties and source functions as inputs; fluxes as outputs.

class(ty_optical_props_arry), intent(in   ) :: optical_props     # Array of ty_optical_props. This type is abstract
                                                                 # and needs to be made concrete, either as an array
                                                                 # (class ty_optical_props_arry) or in some user-defined way
logical,                      intent(in   ) :: top_at_1          # Is the top of the domain at index 1?
                                                                 # (if not, ordering is bottom-to-top)
type(ty_source_func_lw),      intent(in   ) :: sources
real(FT), dimension(:,:),     intent(in   ) :: sfc_emis    # emissivity at surface [] (nband, ncol)
class(ty_fluxes),             intent(inout) :: fluxes      # Array of ty_fluxes. Default computes broadband fluxes at all levels
                                                           #   if output arrays are defined. Can be extended per user desires.
real(FT), dimension(:,:),   &
          target, optional, intent(in   ) :: inc_flux    # incident flux at domain top [W/m2] (ncol, ngpts)
integer,          optional, intent(in   ) :: n_gauss_angles # Number of angles used in Gaussian quadrature
                                                            # (no-scattering solution)
"""
function rte_lw!(optical_props, top_at_1,
                sources, sfc_emis,
                fluxes,
                inc_flux=nothing, n_gauss_angles=nothing)

  #
  # Weights and angle secants for first order (k=1) Gaussian quadrature.
  #   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
  #   after Abramowitz & Stegun 1972, page 921
  #
  FT = eltype(optical_props.band_lims_wvn)

  max_gauss_pts = Integer.(4)
  # Diffusivity angle, not Gaussian angle
  gauss_Ds  = reshape([1.66000000, 0.00000000, 0.00000000, 0.00000000,
                       1.18350343, 2.81649655, 0.00000000, 0.00000000,
                       1.09719858, 1.69338507, 4.70941630, 0.00000000,
                       1.06056257, 1.38282560, 2.40148179, 7.15513024],
                       max_gauss_pts, max_gauss_pts)

  gauss_wts = reshape([0.5000000000, 0.0000000000, 0.0000000000, 0.0000000000,
                       0.3180413817, 0.1819586183, 0.0000000000, 0.0000000000,
                       0.2009319137, 0.2292411064, 0.0698269799, 0.0000000000,
                       0.1355069134, 0.2034645680, 0.1298475476, 0.0311809710],
                       max_gauss_pts, max_gauss_pts)


  # Error checking
  #   if inc_flux is present it has the right dimensions, is positive definite
  ncol  = get_ncol(optical_props)
  nlay  = get_nlay(optical_props)
  ngpt  = get_ngpt(optical_props)
  nband = get_nband(optical_props)

  # Error checking -- consistency of sizes and validity of values
  @assert are_desired(fluxes)
  @assert get_ncol(sources) == ncol
  @assert get_nlay(sources) == nlay
  @assert get_ngpt(sources) == ngpt
  @assert all(size(sfc_emis) .== (nband, ncol))
  @assert !any_vals_outside(sfc_emis, FT.(0.), FT.(1.0))
  if inc_flux ≠ nothing
    @assert all(size(inc_flux) .== (ncol, ngpt))
    @assert !any_vals_less_than(inc_flux, FT(0.0))
  end

  # Number of quadrature points for no-scattering calculation
  n_quad_angs = 1
  if n_gauss_angles ≠ nothing
    @assert !(n_gauss_angles > max_gauss_pts)
    @assert !(n_gauss_angles < 1)
    n_quad_angs = n_gauss_angles
  end

  # Ensure values of tau, ssa, and g are reasonable
  validate!(optical_props)

  #
  #    Lower boundary condition -- expand surface emissivity by band to gpoints
  #

  gpt_flux_up  = Array{FT}(undef,ncol,nlay+1,ngpt)
  gpt_flux_dn  = Array{FT}(undef,ncol,nlay+1,ngpt)
  sfc_emis_gpt = Array{FT}(undef,ncol,       ngpt)

  sfc_emis_gpt = expand_and_transpose(optical_props, sfc_emis)

  #   Upper boundary condition
  if inc_flux ≠ nothing
    gpt_flux_dn = apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux)
  else
    # Default is zero incident diffuse flux
    gpt_flux_dn = apply_BC(ncol, nlay, ngpt, top_at_1, FT)
  end

  # Compute the radiative transfer...
  @assert !(optical_props isa ty_optical_props_nstr)

  if optical_props isa ty_optical_props_1scl

    # No scattering two-stream calculation
    lw_solver_noscat_GaussQuad!(ncol, nlay, ngpt, top_at_1,
                            n_quad_angs, gauss_Ds[1:n_quad_angs,n_quad_angs], gauss_wts[1:n_quad_angs,n_quad_angs],
                            optical_props.tau,
                            sources.lay_source, sources.lev_source_inc, sources.lev_source_dec,
                            sfc_emis_gpt, sources.sfc_source,
                            gpt_flux_up, gpt_flux_dn)

  elseif optical_props isa ty_optical_props_2str
    # two-stream calculation with scattering
    lw_solver_2stream!(ncol, nlay, ngpt, top_at_1,
                               optical_props.tau, optical_props.ssa, optical_props.g,
                               sources.lay_source, sources.lev_source_inc, sources.lev_source_dec, sfc_emis_gpt, sources.sfc_source,
                               gpt_flux_up, gpt_flux_dn)
  end

  reduce!(fluxes,gpt_flux_up, gpt_flux_dn, optical_props, top_at_1)
end

# Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)

function expand_and_transpose(ops::ty_optical_props,arr_in::Array{FT}) where FT
#    class(ty_optical_props),  intent(in ) :: ops
#    real(FT), dimension(:,:), intent(in ) :: arr_in  # (nband, ncol)
#    real(FT), dimension(:,:), intent(out) :: arr_out # (ncol, igpt)

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
