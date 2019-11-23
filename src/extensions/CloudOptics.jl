"""
    mo_cloud_optics

Provides cloud optical properties as a function of effective radius for the RRTMGP bands
  Based on Mie calculations for liquid
    and results from doi:10.1175/JAS-D-12-039.1 for ice with variable surface roughness
  Can use either look-up tables or Pade approximates according to which data has been loaded
  Mike Iacono (AER) is the original author
The class can be used as-is but is also intended as an example of how to extend the RTE framework
"""
module mo_cloud_optics

using OffsetArrays
using ..mo_optical_props
using ..mo_util_array
using ..fortran_intrinsics

export cloud_optics!
export get_min_radius, get_max_radius

export PadeMethod, LookUpTable
export ty_cloud_optics_lut
export ty_cloud_optics_pade

abstract type AbstractInterpolationMethod end

"""
    PadeMethod <: AbstractInterpolationMethod

# Pade coefficients
 - `ext` extinction
 - `ssa` single scattering albedo
 - `asy` asymmetry parameter

Pade particle size regime boundaries
 - `sizreg_ext`
 - `sizreg_ssa`
 - `sizreg_asy`
"""
struct PadeMethod{FT} <: AbstractInterpolationMethod
  # Pade approximation coefficients
  ext::Array{FT}
  ssa::Array{FT}
  asy::Array{FT}
  # Particle size regimes for Pade formulations
  sizreg_ext::Array{FT}
  sizreg_ssa::Array{FT}
  sizreg_asy::Array{FT}
  rad_lwr::FT
  rad_upr::FT
  function PadeMethod(ext::Array{FT},
                      ssa::Array{FT},
                      asy::Array{FT},
                      sizreg_ext::Array{FT},
                      sizreg_ssa::Array{FT},
                      sizreg_asy::Array{FT}) where FT

    # # Error checking
    # nbnd         = size(ext, 1)
    # nsizereg     = size(ext, 2)
    # ncoeff_ssa_g = size(ssa, 3)
    # nbound       = length(sizreg_ext)
    # # Error checking
    # @assert nbnd == get_nband(this)
    # @assert all(size(ssa) .== (nbnd, nsizereg, ncoeff_ssa_g))
    # @assert all(size(asy) .== (nbnd, nsizereg, ncoeff_ssa_g))
    # @assert length(liq.sizreg_ssa) == nbound
    # @assert length(liq.sizreg_asy) == nbound
    # @assert nsizereg == 3

    # this.liq = liq
    # this.ice = ice
    # # Consistency among size regimes
    # @assert !any([liq.sizreg_ssa[1], liq.sizreg_asy[1]] .< this.liq.rad_lwr)
    # @assert !any([ice.sizreg_ssa[1], ice.sizreg_asy[1]] .< this.ice.rad_lwr)
    # @assert !any([liq.sizreg_ssa[nbound], liq.sizreg_asy[nbound]] .> this.liq.rad_upr)
    # @assert !any([ice.sizreg_ssa[nbound], ice.sizreg_asy[nbound]] .> this.ice.rad_upr)

    rad_lwr = sizreg_ext[1]
    rad_upr = sizreg_ext[end]
    return new{FT}(ext,ssa,asy,sizreg_ext,sizreg_ssa,sizreg_asy,rad_lwr,rad_upr)
  end
end

"""
    LookUpTable{FT,I} <: AbstractInterpolationMethod

Lookup table interpolation constants

  - `rad_lwr` particle size lower bound for interpolation
  - `rad_upr` particle size upper bound for interpolation
  - `rad_fac` constant for calculating interpolation indices

Lookup table coefficients
   `ext` extinction
   `ssa` single scattering albedo
   `asy` asymmetry parameter
"""
struct LookUpTable{FT,I} <: AbstractInterpolationMethod
  rad_lwr::FT
  rad_upr::FT
  rad_fac::FT
  nsteps::I
  step_size::FT
  ext::Array{FT}
  ssa::Array{FT}
  asy::Array{FT}
  function LookUpTable(rad_lwr::FT,
                       rad_upr::FT,
                       rad_fac::FT,
                       ext::Array{FT},
                       ssa::Array{FT},
                       asy::Array{FT}) where {FT<:AbstractFloat}

  # # Error checking
  # #
  # # LUT coefficient dimensions
  # #
  # nsize_liq = size(liq.ext,1)
  # nsize_ice = size(ice.ext,1)
  # nbnd      = size(liq.ext,2)
  # nrghice   = size(ice.ext,3)
  # #
  # # Error checking
  # #   Can we check for consistency between table bounds and _fac?
  # #
  # @assert nbnd == get_nband(this)
  # @assert size(ice.ext, 2) == nbnd
  # @assert all(size(liq.ssa) .== (nsize_liq, nbnd))
  # @assert all(size(liq.asy) .== (nsize_liq, nbnd))
  # @assert all(size(ice.ssa) .== (nsize_ice, nbnd, nrghice))
  # @assert all(size(ice.asy) .== (nsize_ice, nbnd, nrghice))

    nsteps = size(ext, 1)
    step_size = (rad_upr - rad_lwr)/FT(nsteps-1)
    return new{FT,Int}(rad_lwr,rad_upr,rad_fac,nsteps,step_size,ext,ssa,asy)
  end
end
get_min_radius(aim::AbstractInterpolationMethod) = aim.rad_lwr
get_max_radius(aim::AbstractInterpolationMethod) = aim.rad_upr

get_num_roughness_types(this::PadeMethod) = size(this.ext, 4)
get_num_roughness_types(this::LookUpTable) = size(this.ext,  3)

struct ty_cloud_optics_lut{FT, I} <: ty_optical_props{FT, I}
  base
  # Ice surface roughness category - needed for Yang (2013) ice optics parameterization
  icergh#::I
  liq#AbstractInterpolationMethod
  ice#AbstractInterpolationMethod
end

struct ty_cloud_optics_pade{FT, I} <: ty_optical_props{FT, I}
  base
  # Ice surface roughness category - needed for Yang (2013) ice optics parameterization
  icergh#::I
  liq#AbstractInterpolationMethod
  ice#AbstractInterpolationMethod
end


# Derive cloud optical properties from provided cloud physical properties

"""
    combine_optical_props!(optical_props::ty_optical_props, ltau, ltaussa, ltaussag, itau, itaussa, itaussag, nbnd, nlay, ncol)

Combine liquid and ice contributions into total cloud optical properties
   See also the `increment!` routines in `mo_optical_props_kernels`
"""
function combine_optical_props!(optical_props::ty_optical_props_1scl, ltau, ltaussa, ltaussag, itau, itaussa, itaussag, nbnd, nlay, ncol)
  # Absorption optical depth  = (1-ssa) * tau = tau - taussa
  for ibnd = 1:nbnd
    for ilay = 1:nlay
      for icol = 1:ncol
        # Absorption optical depth  = (1-ssa) * tau = tau - taussa
        optical_props.tau[icol,ilay,ibnd] = (ltau[icol,ilay,ibnd] - ltaussa[icol,ilay,ibnd]) +
                                            (itau[icol,ilay,ibnd] - itaussa[icol,ilay,ibnd])
      end
    end
  end
end
function combine_optical_props!(optical_props::ty_optical_props_2str{FT}, ltau, ltaussa, ltaussag, itau, itaussa, itaussag, nbnd, nlay, ncol) where FT
  for ibnd = 1:nbnd
    for ilay = 1:nlay
      for icol = 1:ncol
        tau    = ltau[icol,ilay,ibnd]    + itau[icol,ilay,ibnd]
        taussa = ltaussa[icol,ilay,ibnd] + itaussa[icol,ilay,ibnd]
        optical_props.g[icol,ilay,ibnd]  = (ltaussag[icol,ilay,ibnd] + itaussag[icol,ilay,ibnd]) / max(eps(FT), taussa...)
        optical_props.ssa[icol,ilay,ibnd] = taussa/max(eps(FT), tau...)
        optical_props.tau[icol,ilay,ibnd] = tau
      end
    end
  end
end

function validate_cloud_optics!(this::ty_optical_props{FT},
                      clwp, ciwp, reliq, reice,
                      optical_props) where FT
  liqmsk = BitArray(clwp .> FT(0))
  icemsk = BitArray(ciwp .> FT(0))
  # Error checking
  @assert bands_are_equal(this, optical_props)
  @assert get_nband(optical_props) == get_ngpt(optical_props)
  @assert !(this.icergh < 1 || this.icergh > get_num_roughness_types(this.ice))
  @assert !any_vals_outside(reliq, liqmsk, this.liq.rad_lwr, this.liq.rad_upr)
  @assert !any_vals_outside(reice, icemsk, this.ice.rad_lwr, this.ice.rad_upr)
  any(liqmsk) && @assert !any_vals_less_than(clwp, liqmsk, FT(0))
  any(icemsk) && @assert !any_vals_less_than(ciwp, icemsk, FT(0))
  return nothing
end

"""
    cloud_optics!(this::ty_cloud_optics_pade{FT},
                      clwp, ciwp, reliq, reice,
                      optical_props)

Compute single-scattering properties

class(ty_cloud_optics), intent(in   ) :: this
real(wp), intent(in   ) :: clwp  (:,:),    ! cloud ice water path    (units?)
                           ciwp  (:,:),    ! cloud liquid water path (units?)
                           reliq (:,:),    ! cloud ice particle effective size (microns)
                           reice (:,:)      ! cloud liquid particle effective radius (microns)
class(ty_optical_props_arry), intent(inout) :: optical_props Dimensions: (ncol,nlay,nbnd)
! ------- Local -------
logical(wl), dimension(size(clwp,1), size(clwp,2)) :: liqmsk, icemsk
real(wp),    dimension(size(clwp,1), size(clwp,2), this%get_nband()) ::
            ltau, ltaussa, ltaussag, itau, itaussa, itaussag
# ------- Local -------
type(ty_optical_props_2str) :: clouds_liq, clouds_ice
integer  :: nsizereg, ibnd, imom
"""
function cloud_optics!(this::ty_cloud_optics_pade{FT},
                      clwp, ciwp, reliq, reice,
                      optical_props) where FT
  ncol,nlay = size(clwp)
  nbnd = get_nband(this)
  liqmsk = BitArray(clwp .> FT(0))
  icemsk = BitArray(ciwp .> FT(0))
  validate_cloud_optics!(this, clwp, ciwp, reliq, reice, optical_props)

  ltau, ltaussa, ltaussag, itau, itaussa, itaussag =
    ntuple(i-> Array{FT}(undef, size(clwp,1), size(clwp,2), get_nband(this)), 6)
  #### Compute cloud optical properties.

  # Cloud optical properties from Pade coefficient method
  #   Hard coded assumptions: order of approximants, three size regimes
  nsizereg = size(this.liq.ext,2)
  compute_all_from_pade!(ncol, nlay, nbnd, nsizereg,
                            liqmsk, clwp, reliq,
                            2, 3, this.liq.sizreg_ext, this.liq.ext,
                            2, 2, this.liq.sizreg_ssa, this.liq.ssa,
                            2, 2, this.liq.sizreg_asy, this.liq.asy,
                            ltau, ltaussa, ltaussag)
  compute_all_from_pade!(ncol, nlay, nbnd, nsizereg,
                            icemsk, ciwp, reice,
                            2, 3, this.ice.sizreg_ext, this.ice.ext[:,:,:,this.icergh],
                            2, 2, this.ice.sizreg_ssa, this.ice.ssa[:,:,:,this.icergh],
                            2, 2, this.ice.sizreg_asy, this.ice.asy[:,:,:,this.icergh],
                            itau, itaussa, itaussag)

  # Copy total cloud properties onto outputs
  combine_optical_props!(optical_props, ltau, ltaussa, ltaussag, itau, itaussa, itaussag, nbnd, nlay, ncol)

end
function cloud_optics!(this::ty_cloud_optics_lut{FT},
                      clwp, ciwp, reliq, reice,
                      optical_props) where FT
  ncol = size(clwp,1)
  nlay = size(clwp,2)
  nbnd = get_nband(this)
  liqmsk = BitArray(clwp .> FT(0))
  icemsk = BitArray(ciwp .> FT(0))
  validate_cloud_optics!(this, clwp, ciwp, reliq, reice, optical_props)

  ltau, ltaussa, ltaussag, itau, itaussa, itaussag =
    ntuple(i-> Array{FT}(undef, size(clwp,1), size(clwp,2), get_nband(this)), 6)

  #### Compute cloud optical properties.

  # Liquid
  compute_all_from_table!(ncol, nlay, nbnd, liqmsk, clwp, reliq,
                              this.liq.nsteps,
                              this.liq.step_size,
                              this.liq.rad_lwr,
                              this.liq.ext,
                              this.liq.ssa,
                              this.liq.asy,
                              ltau, ltaussa, ltaussag)
  # Ice
  compute_all_from_table!(ncol, nlay, nbnd, icemsk, ciwp, reice,
                              this.ice.nsteps,
                              this.ice.step_size,
                              this.ice.rad_lwr,
                              this.ice.ext[:,:,this.icergh],
                              this.ice.ssa[:,:,this.icergh],
                              this.ice.asy[:,:,this.icergh],
                              itau, itaussa, itaussag)

  # Copy total cloud properties onto outputs
  combine_optical_props!(optical_props, ltau, ltaussa, ltaussag, itau, itaussa, itaussag, nbnd, nlay, ncol)

end

#
# Linearly interpolate values from a lookup table with "nsteps" evenly-spaced
#   elements starting at "offset." The table's second dimension is band.
# Returns 0 where the mask is false.
# We could also try gather/scatter for efficiency
#
function compute_all_from_table!(ncol, nlay, nbnd, mask, lwp, re,
                                  nsteps, step_size, offset,
                                  tau_table, ssa_table, asy_table,
                                  tau::Array{FT}, taussa, taussag) where FT
  # integer,                                intent(in) :: ncol, nlay, nbnd, nsteps
  # logical(wl), dimension(ncol,nlay),      intent(in) :: mask
  # real(wp),    dimension(ncol,nlay),      intent(in) :: lwp, re
  # real(wp),                               intent(in) :: step_size, offset
  # real(wp),    dimension(nsteps,   nbnd), intent(in) :: tau_table, ssa_table, asy_table
  # real(wp),    dimension(ncol,nlay,nbnd)             :: tau, taussa, taussag
  # # ---------------------------
  # integer  :: icol, ilay, ibnd
  # integer  :: index
  # real(wp) :: fint
  # real(wp) :: t, ts, tsg  # tau, tau*ssa, tau*ssa*g
  # # ---------------------------
  #$acc parallel loop gang vector default(present) collapse(3)
  for ibnd = 1:nbnd
    for ilay = 1:nlay
      for icol = 1:ncol
        if mask[icol,ilay]
          index = convert(Int, min(floor((re[icol,ilay] - offset)/step_size)+1, nsteps-1))
          fint = (re[icol,ilay] - offset)/step_size - (index-1)
          t   = lwp[icol,ilay] *                  (tau_table[index,  ibnd] + fint * (tau_table[index+1,ibnd] - tau_table[index,ibnd]))
          ts  = t              *                  (ssa_table[index,  ibnd] + fint * (ssa_table[index+1,ibnd] - ssa_table[index,ibnd]))
          taussag[icol,ilay,ibnd] = ts          * (asy_table[index,  ibnd] + fint * (asy_table[index+1,ibnd] - asy_table[index,ibnd]))
          taussa[icol,ilay,ibnd] = ts
          tau[icol,ilay,ibnd] = t
        else
          tau[icol,ilay,ibnd] = FT(0)
          taussa[icol,ilay,ibnd] = FT(0)
          taussag[icol,ilay,ibnd] = FT(0)
        end
      end
    end
  end
end

#
# Pade functions
#

function compute_all_from_pade!(ncol, nlay, nbnd, nsizes,
                                 mask, lwp, re,
                                 m_ext, n_ext, re_bounds_ext, coeffs_ext,
                                 m_ssa, n_ssa, re_bounds_ssa, coeffs_ssa,
                                 m_asy, n_asy, re_bounds_asy, coeffs_asy,
                                 tau::Array{FT}, taussa, taussag) where FT
  # integer,                        intent(in) :: ncol, nlay, nbnd, nsizes
  # logical(wl),
  #           dimension(ncol,nlay), intent(in) :: mask
  # real(wp), dimension(ncol,nlay), intent(in) :: lwp, re
  # real(wp), dimension(nsizes+1),  intent(in) :: re_bounds_ext, re_bounds_ssa, re_bounds_asy
  # integer,                        intent(in) :: m_ext, n_ext
  # real(wp), dimension(nbnd,nsizes,0:m_ext+n_ext),
  #                                 intent(in) :: coeffs_ext
  # integer,                        intent(in) :: m_ssa, n_ssa
  # real(wp), dimension(nbnd,nsizes,0:m_ssa+n_ssa),
  #                                 intent(in) :: coeffs_ssa
  # integer,                        intent(in) :: m_asy, n_asy
  # real(wp), dimension(nbnd,nsizes,0:m_asy+n_asy),
  #                                 intent(in) :: coeffs_asy
  # real(wp), dimension(ncol,nlay,nbnd)        :: tau, taussa, taussag
  # # ---------------------------
  # integer  :: icol, ilay, ibnd, irad, count
  # real(wp) :: t, ts

  #$acc parallel loop gang vector default(present) collapse(3)
  for ibnd = 1:nbnd
    for ilay = 1:nlay
      for icol = 1:ncol
        if mask[icol,ilay]
          #
          # Finds index into size regime table
          # This works only if there are precisely three size regimes (four bounds) and it's
          #   previously guaranteed that size_bounds(1) <= size <= size_bounds(4)
          #
          irad = convert(Int, min(floor((re[icol,ilay] - re_bounds_ext[2])/re_bounds_ext[3])+2, 3))
          t   = lwp[icol,ilay] * pade_eval(ibnd, nbnd, nsizes, m_ext, n_ext, irad, re[icol,ilay], coeffs_ext)

          irad = convert(Int,min(floor((re[icol,ilay] - re_bounds_ssa[2])/re_bounds_ssa[3])+2, 3))
          # Pade approximants for co-albedo can sometimes be negative
          ts   = t               * (FT(1) - max(FT(0), pade_eval(ibnd, nbnd, nsizes, m_ssa, n_ssa, irad, re[icol,ilay], coeffs_ssa)))
          irad = convert(Int, min(floor((re[icol,ilay] - re_bounds_asy[2])/re_bounds_asy[3])+2, 3))
          taussag[icol,ilay,ibnd] = ts             * pade_eval(ibnd, nbnd, nsizes, m_asy, n_asy, irad, re[icol,ilay], coeffs_asy)

          taussa[icol,ilay,ibnd] = ts
          tau[icol,ilay,ibnd] = t
        else
          tau[icol,ilay,ibnd] = FT(0)
          taussa[icol,ilay,ibnd] = FT(0)
          taussag[icol,ilay,ibnd] = FT(0)
        end
      end
    end
  end

end

#
# Ancillary functions
#
#

#
# Pade functions
#

#
# Evaluate Pade approximant of order [m/n]
#
function pade_eval(nbnd, nrads, m, n, irad, re, pade_coeffs)
  # integer,                intent(in) :: nbnd, nrads, m, n, irad
  # real(wp), dimension(nbnd, nrads, 0:m+n), &
  #                         intent(in) :: pade_coeffs
  # real(wp),               intent(in) :: re
  # real(wp), dimension(nbnd)          :: pade_eval_nbnd

  # integer :: iband
  # real(wp) :: numer, denom
  # integer  :: i
  FT = eltype(pade_coeffs)
  res = Vector{FT}(undef, nbnd)

  for iband = 1:nbnd
    denom = pade_coeffs[iband,irad,n+m]
    for i = n-1+m:-1:1+m
      denom = pade_coeffs[iband,irad,i]+re*denom
    end
    denom =  FT(1)                     +re*denom

    numer = pade_coeffs[iband,irad,m]
    for i = m-1:-1:1
      numer = pade_coeffs[iband,irad,i]+re*numer
    end
    numer = pade_coeffs[iband,irad,0]  +re*numer

    res[iband] = numer/denom
  end
  return res
end

#
# Evaluate Pade approximant of order [m/n]
#
function pade_eval(iband, nbnd, nrads, m, n, irad, re, pade_coeffs_array)
  # !$acc routine seq
  # !
  # integer,                intent(in) :: iband, nbnd, nrads, m, n, irad
  # real(wp), dimension(nbnd, nrads, 0:m+n), &
  #                         intent(in) :: pade_coeffs
  # real(wp),               intent(in) :: re
  # real(wp)                           :: pade_eval_1

  # real(wp) :: numer, denom
  # integer  :: i
  FT = eltype(re)
  pade_coeffs = OffsetArray{FT}(undef, 1:nbnd, 1:nrads, 0:m+n)
  pade_coeffs[:,:,:] = pade_coeffs_array[:,:,:]

  denom = pade_coeffs[iband,irad,n+m]
  for i = n-1+m:-1:1+m
    denom = pade_coeffs[iband,irad,i]+re*denom
  end
  denom =  FT(1)                     +re*denom

  numer = pade_coeffs[iband,irad,m]
  for i = m-1:-1:1
    numer = pade_coeffs[iband,irad,i]+re*numer
  end
  numer = pade_coeffs[iband,irad,0]  +re*numer

  pade_eval_1 = numer/denom
  return pade_eval_1
end

end
