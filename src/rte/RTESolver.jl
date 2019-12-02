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

using ..ArrayUtilities
using ..OpticalProps
using ..Fluxes
using ..FortranIntrinsics

export rte_lw!, rte_sw!, expand_and_transpose

"""
    rte_sw!(atmos::AbstractOpticalPropsArry,
            top_at_1::Bool,
            μ_0::Array{FT},
            inc_flux,
            sfc_alb_dir,
            sfc_alb_dif,
            fluxes,
            inc_flux_dif=nothing) where {FT<:AbstractFloat}

Compute fluxes `fluxes` (`AbstractFluxes`), given

 - `atmos` Optical properties provided as arrays (`AbstractOpticalPropsArry`)
 - `top_at_1` bool indicating if the top of the domain is at index 1
             (if not, ordering is bottom-to-top)
 - `μ_0` cosine of solar zenith angle (ncol)
 - `inc_flux` incident flux at top of domain [W/m2] (ncol, ngpt)
 - `sfc_alb_dir` surface albedo for direct and
 - `sfc_alb_dif` diffuse radiation (nband, ncol)
and, optionally,
 - `inc_flux_dif` incident diffuse flux at top of domain [W/m2] (ncol, ngpt)
"""
function rte_sw!(atmos::AbstractOpticalPropsArry,
                 top_at_1::Bool,
                 μ_0::Array{FT},
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
  @assert size(μ_0, 1) ==  ncol
  @assert !any_vals_outside(μ_0, FT(0), FT(1))
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

  gpt_flux_up  = zeros(FT, ncol, nlay+1, ngpt)
  gpt_flux_dn  = zeros(FT, ncol, nlay+1, ngpt)
  gpt_flux_dir = zeros(FT, ncol, nlay+1, ngpt)

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
  apply_BC!(gpt_flux_dir, ncol, nlay, ngpt, top_at_1,   inc_flux, μ_0)

  if present(inc_flux_dif)
    apply_BC!(gpt_flux_dn, ncol, nlay, ngpt, top_at_1, inc_flux_dif)
  else
    apply_BC!(gpt_flux_dn, ncol, nlay, ngpt, top_at_1)
  end

  if atmos isa OneScalar
    # Direct beam only
    validate!(atmos)
    sw_solver_noscat!(ncol, nlay, ngpt, top_at_1,
                          atmos.τ, μ_0,
                          gpt_flux_dir)
    # No diffuse flux
    # gpt_flux_up .= FT(0)
    # gpt_flux_dn .= FT(0)
  elseif atmos isa TwoStream
    # two-stream calculation with scattering
    validate!(atmos)

    sw_solver_2stream!(ncol, nlay, ngpt, top_at_1,
                           atmos.τ, atmos.ssa, atmos.g, μ_0,
                           sfc_alb_dir_gpt, sfc_alb_dif_gpt,
                           gpt_flux_up, gpt_flux_dn, gpt_flux_dir)
  end

  # ...and reduce spectral fluxes to desired output quantities
  reduce!(fluxes,gpt_flux_up, gpt_flux_dn, atmos, top_at_1, gpt_flux_dir)
end

"""
    rte_lw!(optical_props::AbstractOpticalPropsArry{FT},
            top_at_1,
            sources,
            sfc_emis,
            fluxes,
            inc_flux=nothing,
            n_gauss_angles=nothing) where FT

Interface using only optical properties and source functions as inputs; fluxes as outputs.

 - `optical_props` optical properties
 - `top_at_1` indicates that the top of the domain is at index 1
 - `sources` radiation sources
 - `sfc_emis` emissivity at surface (nband, ncol)
 - `fluxes` Array of AbstractFluxes. Default computes broadband fluxes at all levels
            if output arrays are defined. Can be extended per user desires.
 - [`inc_flux`] incident flux at domain top [W/m2] (ncol, ngpts)
 - [`n_gauss_angles`] Number of angles used in Gaussian quadrature (no-scattering solution)
"""
function rte_lw!(optical_props::AbstractOpticalPropsArry{FT},
                 top_at_1,
                 sources,
                 sfc_emis,
                 fluxes,
                 inc_flux=nothing,
                 n_gauss_angles=nothing) where FT

  #
  # Weights and angle secants for first order (k=1) Gaussian quadrature.
  #   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
  #   after Abramowitz & Stegun 1972, page 921
  #

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

  # Ensure values of τ, ssa, and g are reasonable
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
    apply_BC!(gpt_flux_dn, ncol, nlay, ngpt, top_at_1, inc_flux)
  else
    # Default is zero incident diffuse flux
    apply_BC!(gpt_flux_dn, ncol, nlay, ngpt, top_at_1)
  end

  # Compute the radiative transfer...
  if optical_props isa OneScalar

    # No scattering two-stream calculation
    lw_solver_noscat_GaussQuad!(ncol, nlay, ngpt, top_at_1,
                            n_quad_angs, gauss_Ds[1:n_quad_angs,n_quad_angs], gauss_wts[1:n_quad_angs,n_quad_angs],
                            optical_props.τ,
                            sources.lay_source, sources.lev_source_inc, sources.lev_source_dec,
                            sfc_emis_gpt, sources.sfc_source,
                            gpt_flux_up, gpt_flux_dn)

  elseif optical_props isa TwoStream
    # two-stream calculation with scattering
    lw_solver_2stream!(ncol, nlay, ngpt, top_at_1,
                               optical_props.τ, optical_props.ssa, optical_props.g,
                               sources.lay_source, sources.lev_source_inc, sources.lev_source_dec, sfc_emis_gpt, sources.sfc_source,
                               gpt_flux_up, gpt_flux_dn)
  end

  reduce!(fluxes,gpt_flux_up, gpt_flux_dn, optical_props, top_at_1)
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