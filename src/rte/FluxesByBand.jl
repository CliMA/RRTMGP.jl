#####
##### Broadband by-band
#####

export FluxesByBand

"""
    FluxesByBand{FT} <: AbstractFluxes{FT}

Contains both broadband and by-band fluxes

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct FluxesByBand{FT} <: AbstractFluxes{FT}
    fluxes_broadband
    "Upward flux"
    bnd_flux_up#::Array{FT,3}
    "Downward flux"
    bnd_flux_dn#::Array{FT,3}
    "Net flux"
    bnd_flux_net#::Array{FT,3}
    "Downward direct flux"
    bnd_flux_dn_dir#::Array{FT,3}
end

"""
    reduce!(this::FluxesByBand{FT},
            gpt_flux_up::Array{FT,3},
            gpt_flux_dn::Array{FT,3},
            spectral_disc::AbstractOpticalProps{FT},
            gpt_flux_dn_dir::Union{Nothing,Array{FT,3}}=nothing)

Reduces fluxes by-band to broadband in `FluxesByBand` `this`, given

 - `gpt_flux_up` fluxes by gpoint [W/m2]
 - `gpt_flux_dn` fluxes by gpoint [W/m2]
 - `spectral_disc` spectral discretization, see [`AbstractOpticalProps`](@ref)
and, optionally,
 - `gpt_flux_dn_dir` direct flux downward
"""
function reduce!(
    this::FluxesByBand{FT},
    gpt_flux_up::Array{FT,3},
    gpt_flux_dn::Array{FT,3},
    spectral_disc::AbstractOpticalProps{FT},
    gpt_flux_dn_dir::Union{Nothing,Array{FT,3}} = nothing,
) where {FT<:AbstractFloat}
    ncol, nlev = size(gpt_flux_up)
    ngpt = get_ngpt(spectral_disc)
    nbnd = get_nband(spectral_disc)
    band_lims = deepcopy(get_band_lims_gpoint(spectral_disc))

  # Compute broadband fluxes
  #   This also checks that input arrays are consistently sized
  #
    reduce!(
        this.fluxes_broadband,
        gpt_flux_up,
        gpt_flux_dn,
        spectral_disc,
        gpt_flux_dn_dir,
    )

  # Check sizes
    @assert size(gpt_flux_up, 3) == ngpt
    this.bnd_flux_up ≠ nothing && @assert all(size(this.bnd_flux_up) .== (
        ncol,
        nlev,
        nbnd,
    ))
    this.bnd_flux_dn ≠ nothing && @assert all(size(this.bnd_flux_dn) .== (
        ncol,
        nlev,
        nbnd,
    ))
    this.bnd_flux_dn_dir ≠ nothing && @assert all(size(this.bnd_flux_dn_dir) .== (
        ncol,
        nlev,
        nbnd,
    ))
    this.bnd_flux_net ≠ nothing && @assert all(size(this.bnd_flux_net) .== (
        ncol,
        nlev,
        nbnd,
    ))
  # Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
    @assert !(this.bnd_flux_dn_dir ≠ nothing && gpt_flux_dn_dir == nothing)

  # Band-by-band fluxes

  # Up flux
    if this.bnd_flux_up ≠ nothing
        this.bnd_flux_up = sum_byband(
            ncol,
            nlev,
            ngpt,
            nbnd,
            band_lims,
            gpt_flux_up,
        )
    end

  # Down flux
    if this.bnd_flux_dn ≠ nothing
        this.bnd_flux_dn = sum_byband(
            ncol,
            nlev,
            ngpt,
            nbnd,
            band_lims,
            gpt_flux_dn,
        )
    end

  # Direct Down flux
    if this.bnd_flux_dn_dir ≠ nothing
        this.bnd_flux_dn_dir = sum_byband(
            ncol,
            nlev,
            ngpt,
            nbnd,
            band_lims,
            gpt_flux_dn_dir,
        )
    end

  # Net flux
    if (this.bnd_flux_net ≠ nothing)
    #
    #  Reuse down and up results if possible
    #
        if (this.bnd_flux_dn ≠ nothing && this.bnd_flux_up ≠ nothing)
            this.bnd_flux_net = net_byband(
                ncol,
                nlev,
                nbnd,
                this.bnd_flux_dn,
                this.bnd_flux_up,
            )
        else
            this.bnd_flux_net = net_byband(
                ncol,
                nlev,
                ngpt,
                nbnd,
                band_lims,
                gpt_flux_dn,
                gpt_flux_up,
            )
        end
    end
end

#####
##### Kernels for computing by-band fluxes by summing
##### over all elements in the spectral dimension.
#####

"""
    sum_byband(ncol, nlev, ngpt, nbnd, band_lims, spectral_flux)

Spectral reduction over all points

integer,                               intent(in ) :: ncol, nlev, ngpt, nbnd
integer,  dimension(2,          nbnd), intent(in ) :: band_lims
real(FT), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux
real(FT), dimension(ncol, nlev, nbnd), intent(out) :: byband_flux
"""
function sum_byband(
    ncol,
    nlev,
    ngpt,
    nbnd,
    band_lims,
    spectral_flux::Array{FT},
) where {FT}

    byband_flux = Array{FT}(undef, ncol, nlev, nbnd)
    for ibnd in 1:nbnd
        for ilev in 1:nlev
            for icol in 1:ncol
                byband_flux[icol, ilev, ibnd] = spectral_flux[
                    icol,
                    ilev,
                    band_lims[1, ibnd],
                ]
                for igpt in band_lims[1, ibnd]+1:band_lims[2, ibnd]
                    byband_flux[icol, ilev, ibnd] = byband_flux[
                        icol,
                        ilev,
                        ibnd,
                    ] + spectral_flux[
                        icol,
                        ilev,
                        igpt,
                    ]
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
function net_byband(
    ncol,
    nlev,
    ngpt,
    nbnd,
    band_lims,
    spectral_flux_dn::Array{FT},
    spectral_flux_up,
) where {FT}

    byband_flux_net = Array{FT}(undef, ncol, nlev, nbnd)

    for ibnd in 1:nbnd
        for ilev in 1:nlev
            for icol in 1:ncol
                igpt = band_lims[1, ibnd]
                byband_flux_net[icol, ilev, ibnd] = spectral_flux_dn[
                    icol,
                    ilev,
                    igpt,
                ] - spectral_flux_up[
                    icol,
                    ilev,
                    igpt,
                ]
                for igpt in band_lims[1, ibnd]+1:band_lims[2, ibnd]
                    byband_flux_net[icol, ilev, ibnd] = byband_flux_net[
                        icol,
                        ilev,
                        ibnd,
                    ] + spectral_flux_dn[
                        icol,
                        ilev,
                        igpt,
                    ] - spectral_flux_up[
                        icol,
                        ilev,
                        igpt,
                    ]
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
net_byband(ncol, nlev, nbnd, byband_flux_dn, byband_flux_up) =
    byband_flux_dn .- byband_flux_up
