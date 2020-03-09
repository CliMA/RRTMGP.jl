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
    "Upward flux"
    flux_up::Array{FT,2}
    "Downward flux"
    flux_dn::Array{FT,2}
    "Net flux"
    flux_net::Array{FT,2}
    "Downward direct flux"
    flux_dn_dir::Union{Array{FT,2},Nothing}
    FluxesBroadBand(FT, s, include_direct::Bool = false) =
        include_direct ? new{FT}(ntuple(i -> Array{FT}(undef, s...), 4)...) :
        new{FT}(ntuple(i -> Array{FT}(undef, s...), 3)..., nothing)
end

"""
    reduce!(this::FluxesBroadBand{FT},
            gpt_flux_up::Array{FT,3},
            gpt_flux_dn::Array{FT,3},
            spectral_disc::AbstractOpticalProps{FT},
            gpt_flux_dn_dir::Union{Nothing,Array{FT,3}}=nothing) where FT<:AbstractFloat

Compute `FluxesBroadBand` `this` by summing over the
spectral dimension, given

 - `gpt_flux_up` upward fluxes by gpoint [W/m2]
 - `gpt_flux_dn` downward fluxes by gpoint [W/m2]
 - `spectral_disc` optical properties containing spectral information
optional:
 - `gpt_flux_dn_dir` downward direct flux
"""
function reduce!(
    this::FluxesBroadBand{FT},
    gpt_flux_up::Array{FT,3},
    gpt_flux_dn::Array{FT,3},
    spectral_disc::AbstractOpticalProps{FT},
    gpt_flux_dn_dir::Union{Nothing,Array{FT,3}} = nothing,
) where {FT<:AbstractFloat}

    ncol, nlev, ngpt = size(gpt_flux_up)

  # Check array sizes
    @assert all(size(gpt_flux_dn) .== (ncol, nlev, ngpt))
    @assert all(size(this.flux_up) .== (ncol, nlev))
    @assert all(size(this.flux_dn) .== (ncol, nlev))
    @assert all(size(this.flux_net) .== (ncol, nlev))
    gpt_flux_dn_dir ≠ nothing && @assert all(size(gpt_flux_dn_dir) .== (
        ncol,
        nlev,
        ngpt,
    ))
    gpt_flux_dn_dir ≠ nothing && @assert all(size(this.flux_dn_dir) .== (
        ncol,
        nlev,
    ))

  # Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
    @assert !(this.flux_dn_dir ≠ nothing && gpt_flux_dn_dir == nothing)

  # Broadband fluxes - call the kernels
    this.flux_up .= sum_broadband(gpt_flux_up)
    this.flux_dn .= sum_broadband(gpt_flux_dn)
    if this.flux_dn_dir ≠ nothing
        this.flux_dn_dir .= sum_broadband(gpt_flux_dn_dir)
    end

  #  Reuse down and up results if possible
    this.flux_net .= net_broadband(this.flux_dn, this.flux_up)
end

#####
##### Kernels for computing broadband fluxes by summing
##### over all elements in the spectral dimension.
#####

"""
    sum_broadband(spectral_flux::Array{FT}) where {FT}

Broadband, computing via spectral reduction over all points
"""
sum_broadband(spectral_flux::Array{FT}) where {FT} =
    sum(spectral_flux, dims = 3)[:, :, 1]

"""
    net_broadband(flux_dn, flux_up)

Net flux when broadband flux up and down are already available
"""
net_broadband(
    flux_dn::Array{FT},
    flux_up::Array{FT},
) where {FT<:AbstractFloat} = flux_dn .- flux_up
