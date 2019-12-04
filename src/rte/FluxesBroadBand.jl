#####
##### Broadband fluxes
#####

export FluxesBroadBand

"""
    FluxesBroadBand{FT} <: AbstractFluxes{FT}

Contains upward, downward, net, and direct downward fluxes

# Fields
$(DocStringExtensions.FIELDS)
"""
struct FluxesBroadBand{FT} <: AbstractFluxes{FT}
  "upward flux"
  flux_up#::Array{FT,2}
  "downward flux"
  flux_dn#::Array{FT,2}
  "net flux"
  flux_net#::Array{FT,2}
  "downward direct flux"
  flux_dn_dir#::Array{FT,2}
  FluxesBroadBand(FT, s, include_direct::Bool=false) =
    include_direct ? new{FT}(ntuple(i->Array{FT}(undef, s...),4)...) :
                     new{FT}(ntuple(i->Array{FT}(undef, s...),3)...,nothing)
end

"""
    reduce!(this::FluxesBroadBand,
                 gpt_flux_up::A,
                 gpt_flux_dn::A,
                 spectral_disc::AbstractOpticalProps,
                 top_at_1::Bool,
                 gpt_flux_dn_dir=nothing)

Compute `FluxesBroadBand` `this` by summing over the
spectral dimension, given

 - `gpt_flux_up` upward fluxes by gpoint [W/m2]
 - `gpt_flux_dn` downward fluxes by gpoint [W/m2]
 - `spectral_disc` optical properties containing spectral information
 - `top_at_1` bool indicating at top
optional:
 - `gpt_flux_dn_dir` downward direct flux
"""
function reduce!(this::FluxesBroadBand,
                 gpt_flux_up::Array{FT,3},
                 gpt_flux_dn::Array{FT,3},
                 spectral_disc::AbstractOpticalProps,
                 top_at_1::Bool,
                 gpt_flux_dn_dir::Union{Nothing,Array{FT,3}}=nothing) where FT<:AbstractFloat

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
    this.flux_dn_dir .= sum_broadband(ncol, nlev, ngpt, gpt_flux_dn_dir)
  end

  if associated(this.flux_net)

    #  Reuse down and up results if possible
    if associated(this.flux_dn) && associated(this.flux_up)
      this.flux_net .= net_broadband(ncol, nlev,      this.flux_dn, this.flux_up)
    else
      this.flux_net .= net_broadband(ncol, nlev, ngpt, gpt_flux_dn,  gpt_flux_up)
    end
  end
end

#####
##### Kernels for computing broadband fluxes by summing
##### over all elements in the spectral dimension.
#####

"""
    sum_broadband(ncol, nlev, ngpt, spectral_flux::Array{FT}) where FT
Spectral reduction over all points

integer,                               intent(in ) :: ncol, nlev, ngpt
real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux
real(FT), dimension(ncol, nlev),       intent(out) :: broadband_flux
"""
sum_broadband(ncol, nlev, ngpt, spectral_flux::Array{FT}) where FT = sum(spectral_flux, dims=3)[:,:,1]

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
