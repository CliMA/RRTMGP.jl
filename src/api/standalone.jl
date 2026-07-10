#####
##### Standalone (classroom) front door
#####
#
# Convenience entry points for using RRTMGP on its own, without a host model. The
# gray-atmosphere path needs neither NCDatasets nor ClimaParams, so it runs in a
# fresh session after `using RRTMGP`.

import ..GrayAtmosphere
import ..AtmosphericStates
import ..VolumeMixingRatios
import ..BCs
using ClimaComms

"""
    RadiationOutput

The result of a standalone radiation solve ([`solve`](@ref) or
[`solve_gray`](@ref)): the broadband fluxes, the heating rate, and the
underlying solver. The field names are stable — teaching material written
against them does not depend on the getter API. The flux fields are views into
the solver's buffers (they change if the solver is re-solved); `heating_rate`
is a freshly allocated array.

# Fields
- `lw_up`, `lw_dn`, `lw_net`: longwave up/down/net flux [W/m²], `(nlev, ncol)`.
- `sw_up`, `sw_dn`, `sw_direct_dn`, `sw_net`: shortwave up/down/direct-beam/net
  flux [W/m²], `(nlev, ncol)`.
- `net`: combined longwave + shortwave net flux [W/m²], `(nlev, ncol)`.
- `heating_rate`: radiative heating rate [K/s], `(nlay, ncol)`.
- `solver`: the `RRTMGPSolver`, for getter access and re-solves.
"""
struct RadiationOutput{A, H, S}
    lw_up::A
    lw_dn::A
    lw_net::A
    sw_up::A
    sw_dn::A
    sw_direct_dn::A
    sw_net::A
    net::A
    heating_rate::H
    solver::S
end

_radiation_output(solver::RRTMGPSolver) = RadiationOutput(
    lw_flux_up(solver),
    lw_flux_dn(solver),
    lw_flux_net(solver),
    sw_flux_up(solver),
    sw_flux_dn(solver),
    sw_direct_flux_dn(solver),
    sw_flux_net(solver),
    net_flux(solver),
    heating_rate(solver),
    solver,
)

"""
    default_parameters(FT)

Return Earth-like RRTMGP physical parameters as an `RRTMGPParameters{FT}`, so the
gray standalone path needs neither NCDatasets nor ClimaParams. Override individual
values by constructing `RRTMGP.Parameters.RRTMGPParameters` directly.
"""
function default_parameters(::Type{FT}) where {FT}
    return RP.RRTMGPParameters{FT}(;
        grav = FT(9.81),
        molmass_dryair = FT(0.02897),
        molmass_water = FT(0.018015),
        gas_constant = FT(8.314462618),
        kappa_d = FT(2 / 7),
        Stefan = FT(5.670374419e-8),
        avogad = FT(6.02214076e23),
    )
end

"""
    heating_rate(s::RRTMGPSolver)

Return the radiative heating rate at each layer [K/s], computed from the net
flux divergence `(g / cₚ) ∂F_net/∂p`. Call `update_fluxes!(s)` first. Allocates
and returns a fresh `(nlay, ncol)` array.
"""
function heating_rate(s::RRTMGPSolver)
    device = ClimaComms.device(s.grid_params)
    (; ncol) = s.grid_params
    FT = eltype(s.grid_params)
    # `net_flux`/`level_pressure` are domain views, so use the domain layer count.
    nlay = s.grid_params.nlay - Int(s.grid_params.isothermal_boundary_layer)
    hr = similar(net_flux(s), nlay, ncol)
    GrayAtmosphere.compute_gray_heating_rate!(
        device,
        hr,
        level_pressure(s),
        ncol,
        nlay,
        net_flux(s),
        FT(RP.cp_d(s.params)),
        FT(RP.grav(s.params)),
    )
    return hr
end

"""
    solve_gray(FT; kwargs...)

Set up and solve a gray-atmosphere radiation problem in a single call. Requires
no NetCDF lookup tables, so it runs standalone after `using RRTMGP`. Builds the
pressure/temperature profile (Schneider-2004-style), the surface boundary
conditions, and an `RRTMGPSolver`, then runs `update_fluxes!`.

# Keyword Arguments
- `context = ClimaComms.context()`: the `ClimaComms` context (CPU or GPU).
- `nlay = 60`: number of (physical) layers. The standalone path adds no isothermal
  boundary layer, so this is the whole grid; the Layer-2 `RRTMGPGridParams` constructor
  names the same physical quantity `domain_nlay`.
- `ncol = 1`: number of columns (ignored if `latitude` is given).
- `latitude = nothing`: latitudes [degrees]; defaults to the equator for a single
  column, or an evenly spaced pole-to-pole range otherwise.
- `surface_pressure = 1.0e5`: surface pressure [Pa].
- `top_pressure = 9.0e3`: top-of-atmosphere pressure [Pa].
- `optical_thickness = GrayOpticalThicknessOGorman2008(FT)`: gray optical-thickness parameters.
- `params = default_parameters(FT)`: RRTMGP physical parameters.
- `surface_emissivity = 1`: longwave surface emissivity [-].
- `cos_zenith = 0.5`: cosine of the solar zenith angle [-].
- `toa_flux = 1361`: top-of-atmosphere solar flux [W/m²].
- `surface_albedo = 0.2`: shortwave surface albedo [-].

# Returns
A [`RadiationOutput`](@ref) with the level fluxes, the layer heating rate, and
the underlying solver.

# Examples
```julia
using RRTMGP
out = RRTMGP.solve_gray(Float64; nlay = 60, ncol = 1)
out.net          # net flux at each level [W/m²]
out.heating_rate # heating rate at each layer [K/s]
```
"""
function solve_gray(
    ::Type{FT};
    context = ClimaComms.context(),
    nlay::Int = 60,
    ncol::Int = 1,
    latitude = nothing,
    surface_pressure = FT(1.0e5),
    top_pressure = FT(9.0e3),
    optical_thickness = AtmosphericStates.GrayOpticalThicknessOGorman2008(FT),
    params = default_parameters(FT),
    surface_emissivity = FT(1),
    cos_zenith = FT(0.5),
    toa_flux = FT(1361),
    surface_albedo = FT(0.2),
) where {FT <: AbstractFloat}
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    lat = if !isnothing(latitude)
        DA{FT}(collect(FT, latitude))
    elseif ncol == 1
        DA{FT}([0])
    else
        DA{FT}(range(FT(-90), FT(90); length = ncol))
    end
    ncol = length(lat)
    grid_params = RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
    as = AtmosphericStates.setup_gray_as_pr_grid(
        context,
        nlay,
        lat,
        FT(surface_pressure),
        FT(top_pressure),
        optical_thickness,
        params,
        DA,
    )
    sfc_emis = fill!(DA{FT}(undef, 1, ncol), FT(surface_emissivity))
    bcs_lw = BCs.LwBCs(sfc_emis, nothing)
    cos_zen = fill!(DA{FT}(undef, ncol), FT(cos_zenith))
    toa = fill!(DA{FT}(undef, ncol), FT(toa_flux))
    alb_dir = fill!(DA{FT}(undef, 1, ncol), FT(surface_albedo))
    alb_dif = fill!(DA{FT}(undef, 1, ncol), FT(surface_albedo))
    bcs_sw = BCs.SwBCs(cos_zen, toa, alb_dir, nothing, alb_dif)
    solver =
        RRTMGPSolver(grid_params, GrayRadiation(), params, bcs_lw, bcs_sw, as)
    update_fluxes!(solver)
    return _radiation_output(solver)
end

"""
    solve(profile::AtmosphereProfile; method = ClearSkyRadiation(false), kwargs...)

Solve the radiative-transfer problem for an [`AtmosphereProfile`](@ref) in a
single call and return a [`RadiationOutput`](@ref). Builds the atmospheric
state, boundary conditions, and an `RRTMGPSolver` internally, then runs
[`update_fluxes!`](@ref) — the one-line standalone path:

```julia
using RRTMGP, NCDatasets
profile = RRTMGP.standard_atmosphere(Float64; kind = :tropical)
out = RRTMGP.solve(profile)
out.net          # net flux at each level [W/m²]
out.heating_rate # heating rate at each layer [K/s]
```

`method` selects the radiation model: `ClearSkyRadiation(false)` (the default;
requires lookup tables, so load NCDatasets first or pass cached `lookups`) or
`GrayRadiation()` (no lookup tables — runs after a bare `using RRTMGP`; the
profile's mixing ratios are ignored). Cloud and aerosol methods are not
supported here — construct an [`RRTMGPSolver`](@ref) directly for those.

# Keyword Arguments
- `method = ClearSkyRadiation(false)`: the radiation method (see above).
- `context = ClimaComms.context()`: the `ClimaComms` context (CPU or GPU); the
  profile's host arrays are copied to the device.
- `params = default_parameters(FT)`: RRTMGP physical parameters.
- `lookups = nothing`: prebuilt [`LookupBundle`](@ref) to reuse (avoids the
  NetCDF read when solving many profiles); `nothing` builds them internally.
- `surface_emissivity = 1`: longwave surface emissivity [-].
- `cos_zenith = 0.5`: cosine of the solar zenith angle [-].
- `toa_flux = 1361`: top-of-atmosphere solar flux [W/m²].
- `surface_albedo = 0.2`: shortwave surface albedo [-].
- `optical_thickness = GrayOpticalThicknessOGorman2008(FT)`: gray
  optical-thickness parameters (`GrayRadiation` only).
"""
function solve(
    profile::AtmosphereProfile;
    method::AbstractRRTMGPMethod = ClearSkyRadiation(false),
    context = ClimaComms.context(),
    params = default_parameters(eltype(profile.p_lay)),
    lookups = nothing,
    surface_emissivity = 1,
    cos_zenith = 0.5,
    toa_flux = 1361,
    surface_albedo = 0.2,
    optical_thickness = AtmosphericStates.GrayOpticalThicknessOGorman2008(
        eltype(profile.p_lay),
    ),
)
    FT = eltype(profile.p_lay)
    (nlay, ncol) = size(profile.p_lay)
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    # Copy (never alias) the profile's host arrays: the solve mutates its state
    # in place (clipping, column amounts), and the profile should stay reusable.
    dev(x) = DA{FT}(copy(x))

    method isa Union{GrayRadiation, ClearSkyRadiation} || error(
        "`solve(profile)` supports `GrayRadiation` and `ClearSkyRadiation`; \
         for cloud or aerosol radiation, construct an `RRTMGPSolver` directly.",
    )
    method isa ClearSkyRadiation && method.aerosol_radiation && error(
        "`solve(profile)` does not support aerosol radiation (an \
         `AtmosphereProfile` carries no aerosol state); use \
         `ClearSkyRadiation(false)` or construct an `RRTMGPSolver` directly.",
    )

    grid_params = RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
    if isnothing(lookups)
        _check_lookup_support(method)
        lookups = lookup_tables(grid_params, method)
    end

    as = if method isa GrayRadiation
        AtmosphericStates.GrayAtmosphericState(
            dev(profile.lat),
            dev(profile.p_lay),
            dev(profile.p_lev),
            dev(profile.t_lay),
            dev(profile.t_lev),
            dev(profile.z_lev),
            dev(profile.t_sfc),
            optical_thickness,
        )
    else
        idx_gases = lookups.idx_gases_lw
        vmr_wm = zeros(FT, lookups.ngas_lw)
        for (gas, val) in profile.well_mixed_vmr
            gas in ("h2o", "o3") && error(
                "`well_mixed_vmr` lists \"$gas\", which is not well-mixed; it \
                 belongs in the profile's `vmr_$gas` field.",
            )
            haskey(idx_gases, gas) || error(
                "unknown gas \"$gas\" in `well_mixed_vmr`; RRTMGP's gas names \
                 are $(sort(collect(keys(idx_gases)))).",
            )
            vmr_wm[idx_gases[gas]] = FT(val)
        end
        vmr = VolumeMixingRatios.VmrGM(
            dev(profile.vmr_h2o),
            dev(profile.vmr_o3),
            DA{FT}(vmr_wm),
        )
        # Row 1 (col_dry) is computed by `prepare_atmosphere!` inside
        # `update_fluxes!`; row 4 (rel_hum) feeds only aerosol optics, which
        # this path rejects above.
        layerdata = zeros(FT, 4, nlay, ncol)
        layerdata[2, :, :] .= profile.p_lay
        layerdata[3, :, :] .= profile.t_lay
        # `lon` is unused by RRTMGP but must match `lat`'s type; `lat` enables
        # latitude-dependent gravity in the column-amount computation.
        AtmosphericStates.AtmosphericState(
            dev(zeros(FT, ncol)),
            dev(profile.lat),
            DA{FT}(layerdata),
            dev(profile.p_lev),
            dev(profile.t_lev),
            dev(profile.t_sfc),
            vmr,
            nothing,
            nothing,
        )
    end

    sfc_emis = fill!(DA{FT}(undef, lookups.nbnd_lw, ncol), FT(surface_emissivity))
    bcs_lw = BCs.LwBCs(sfc_emis, nothing)
    cos_zen = fill!(DA{FT}(undef, ncol), FT(cos_zenith))
    toa = fill!(DA{FT}(undef, ncol), FT(toa_flux))
    alb_dir = fill!(DA{FT}(undef, lookups.nbnd_sw, ncol), FT(surface_albedo))
    alb_dif = fill!(DA{FT}(undef, lookups.nbnd_sw, ncol), FT(surface_albedo))
    bcs_sw = BCs.SwBCs(cos_zen, toa, alb_dir, nothing, alb_dif)

    solver =
        RRTMGPSolver(grid_params, method, params, bcs_lw, bcs_sw, as; lookups)
    update_fluxes!(solver)
    return _radiation_output(solver)
end
