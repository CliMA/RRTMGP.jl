"""
    mo_fluxes_byband_kernels

Kernels for computing byband fluxes by summing
over all elements in the spectral dimension.
"""
module mo_fluxes_byband_kernels

export sum_byband, net_byband

"""
    sum_byband(ncol, nlev, ngpt, nbnd, band_lims, spectral_flux)

Spectral reduction over all points

integer,                               intent(in ) :: ncol, nlev, ngpt, nbnd
integer,  dimension(2,          nbnd), intent(in ) :: band_lims
real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux
real(FT), dimension(ncol, nlev, nbnd), intent(out) :: byband_flux
"""
function sum_byband(ncol, nlev, ngpt, nbnd, band_lims, spectral_flux::Array{FT}) where FT

  byband_flux = Array{FT}(undef, ncol, nlev, nbnd)
  for ibnd = 1:nbnd
    for ilev = 1:nlev
      for icol = 1:ncol
        byband_flux[icol, ilev, ibnd] =  spectral_flux[icol, ilev, band_lims[1, ibnd]]
        for igpt = band_lims[1,ibnd]+1: band_lims[2,ibnd]
          byband_flux[icol, ilev, ibnd] = byband_flux[icol, ilev, ibnd] + spectral_flux[icol, ilev, igpt]
        end
      end
    end
  end
  return byband_flux
end

"""
    net_byband(ncol, nlev, ngpt, nbnd, band_lims, spectral_flux_dn, spectral_flux_up)

Net flux: Spectral reduction over all points

integer,                               intent(in ) :: ncol, nlev, ngpt, nbnd
integer,  dimension(2,          nbnd), intent(in ) :: band_lims
real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux_dn, spectral_flux_up
real(FT), dimension(ncol, nlev, nbnd), intent(out) :: byband_flux_net
"""
function net_byband(ncol, nlev, ngpt, nbnd, band_lims, spectral_flux_dn::Array{FT}, spectral_flux_up) where FT

  byband_flux_net = Array{FT}(undef, ncol, nlev, nbnd)

  for ibnd = 1:nbnd
    for ilev = 1:nlev
      for icol = 1:ncol
        igpt = band_lims[1,ibnd]
        byband_flux_net[icol, ilev, ibnd] = spectral_flux_dn[icol, ilev, igpt] -
                                            spectral_flux_up[icol, ilev, igpt]
        for igpt = band_lims[1,ibnd]+1:band_lims[2,ibnd]
          byband_flux_net[icol, ilev, ibnd] = byband_flux_net[icol, ilev, ibnd] +
                                              spectral_flux_dn[icol, ilev, igpt] -
                                              spectral_flux_up[icol, ilev, igpt]
        end
      end
    end
  end
  return byband_flux_net
end

"""
    net_byband(ncol, nlev, nbnd, byband_flux_dn, byband_flux_up)

integer,                               intent(in ) :: ncol, nlev, nbnd
real(FT), dimension(ncol, nlev, nbnd), intent(in ) :: byband_flux_dn, byband_flux_up
real(FT), dimension(ncol, nlev, nbnd), intent(out) :: byband_flux_net
"""
net_byband(ncol, nlev, nbnd, byband_flux_dn, byband_flux_up) = byband_flux_dn .- byband_flux_up

end # module
