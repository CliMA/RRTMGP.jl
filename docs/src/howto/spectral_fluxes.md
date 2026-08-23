# How to get per-band (spectral) fluxes

By default, RRTMGP retains only broadband fluxes. For applications that need the
spectral decomposition, the solver can retain per-band fluxes. Examples include
atmospheric chemistry, where photolysis rates depend on the ultraviolet and
visible fluxes; land models, where the radiation a canopy absorbs divides into
photosynthetically active (visible) and near-infrared bands with very different
vegetation albedos; band-by-band forcing diagnostics; and satellite-channel
proxies.

## Request per-band fluxes at construction

```julia
solver = RRTMGP.RRTMGPSolver(
    grid_params, method, params, bcs_lw, bcs_sw, as;
    lookups,
    spectral_fluxes = true,
)
RRTMGP.update_fluxes!(solver)
```

This is supported for the spectral (non-gray) methods — gray radiation is a
single band, and asking for per-band fluxes there is an error at construction.
The shortwave bands are always retained, because spectral shortwave optics must
be two-stream anyway. The longwave bands are retained only when the longwave is
two-stream: the no-scattering (`OneScalar`) solver does not accumulate per band.
So a mixed configuration — a no-scattering longwave with a two-stream shortwave,
as E3SM and ERF run — retains the shortwave bands, while the `spectral_lw_*`
getters keep raising their error. The per-band buffers add `(nlev, ncol, n_bnd)`
arrays per band set, so they are opt-in.

## Read them

```julia
F_lw = RRTMGP.spectral_lw_flux_up(solver)   # (nlev, ncol, nbnd_lw) view
F_sw = RRTMGP.spectral_sw_flux_dn(solver)   # (nlev, ncol, nbnd_sw) view
```

Band `b`'s slice `F_lw[:, :, b]` has the same layout as the broadband getters,
and the bands sum to the broadband flux:

```julia
sum(RRTMGP.spectral_lw_flux_up(solver); dims = 3) ≈ RRTMGP.lw_flux_up(solver)
```

The `spectral_*_flux_net` getters are views into a retained per-band net-flux
buffer, updated on every solve like the `up`/`dn` buffers.

## The shortwave direct beam, per band

The shortwave additionally retains the per-band direct beam, so the surface
downward flux can be split into direct and diffuse within each band:

```julia
dir = RRTMGP.spectral_sw_direct_flux_dn(solver)      # (nlev, ncol, nbnd_sw)
dif = RRTMGP.spectral_sw_flux_dn(solver) .- dir      # the diffuse remainder
```

It sums over bands to the broadband direct beam, like the other band buffers:

```julia
sum(RRTMGP.spectral_sw_direct_flux_dn(solver); dims = 3) ≈
    RRTMGP.sw_direct_flux_dn(solver)
```

This is what a land-surface scheme needs: Noah-MP and CLM-family canopies apply
different albedos to the direct and the diffuse beam, and to the visible and
near-infrared halves of the spectrum, so they consume four surface numbers
(direct/diffuse x visible/near-IR) rather than one broadband flux. Sum the
per-band surface values over the bands on each side of the 0.7 µm
(14286 cm⁻¹) boundary, using [`sw_band_bounds`](@ref RRTMGP.sw_band_bounds) to
classify the bands and splitting the one band that straddles it.

## Identify the bands

[`lw_band_bounds`](@ref RRTMGP.lw_band_bounds) and [`sw_band_bounds`](@ref
RRTMGP.sw_band_bounds) return the `(2, n_bnd)` wavenumber edges (cm⁻¹) of each
band:

```julia
wn = RRTMGP.lw_band_bounds(solver)
wn[:, 1]   # lower/upper wavenumber of longwave band 1
```

The RRTMGP longwave tables have 16 bands and the shortwave tables 14, following
[pincus2019](@citet).
