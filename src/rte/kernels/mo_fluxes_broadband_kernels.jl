"""
    mo_fluxes_broadband_kernels

Kernels for computing broadband fluxes by summing
over all elements in the spectral dimension.
"""
module mo_fluxes_broadband_kernels

export sum_broadband, net_broadband

"""
    sum_broadband(ncol, nlev, ngpt, spectral_flux::Array{FT}) where FT
Spectral reduction over all points

integer,                               intent(in ) :: ncol, nlev, ngpt
real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux
real(FT), dimension(ncol, nlev),       intent(out) :: broadband_flux
"""
sum_broadband(ncol, nlev, ngpt, spectral_flux::Array{FT}) where FT = sum(spectral_flux, dims=3)

# net_broadband_new(ncol, nlev, ngpt, spectral_flux_dn::Array{FT}, spectral_flux_up) where FT = sum(spectral_flux_dn .- spectral_flux_up, dims=3)
"""
    net_broadband(ncol, nlev, ngpt, spectral_flux_dn::Array{FT}, spectral_flux_up) where FT

Net flux: Spectral reduction over all points

integer,                               intent(in ) :: ncol, nlev, ngpt
real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux_dn, spectral_flux_up
real(FT), dimension(ncol, nlev),       intent(out) :: broadband_flux_net
"""
function net_broadband(ncol, nlev, ngpt, spectral_flux_dn::Array{FT}, spectral_flux_up) where FT

  broadband_flux_net = Array{FT}(undef, ncol, nlev)
  for ilev = 1:nlev
    for icol = 1:ncol
      diff_ = spectral_flux_dn[icol, ilev, 1   ] - spectral_flux_up[icol, ilev,     1]
      broadband_flux_net[icol, ilev] = diff_
    end
  end
  for igpt = 2:ngpt
    for ilev = 1:nlev
      for icol = 1:ncol
        diff_ = spectral_flux_dn[icol, ilev, igpt] - spectral_flux_up[icol, ilev, igpt]
        broadband_flux_net[icol, ilev] = broadband_flux_net[icol, ilev] + diff_
      end
    end
  end
  return broadband_flux_net
end

"""
    net_broadband(ncol, nlev, flux_dn, flux_up)

Net flux when bradband flux up and down are already available
"""
net_broadband(ncol, nlev, flux_dn, flux_up) = flux_dn .- flux_up

end
