#####
##### Standalone (classroom) front door
#####
#
# Convenience entry points for using RRTMGP on its own, without a host model. The
# gray-atmosphere path needs neither NCDatasets nor ClimaParams, so it runs in a
# fresh session after `using RRTMGP`.

import ..GrayUtils
import ..AtmosphericStates
import ..BCs
using ClimaComms

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
    # `net_flux`/`face_pressure` are domain views, so use the domain layer count.
    nlay = s.grid_params.nlay - Int(s.grid_params.isothermal_boundary_layer)
    hr = similar(net_flux(s), nlay, ncol)
    GrayUtils.compute_gray_heating_rate!(
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
A `NamedTuple` with the level fluxes `lw_up`, `lw_dn`, `lw_net`, `sw_up`, `sw_dn`,
`sw_net`, and `net` (each `(nlev, ncol)` [W/m²]), the layer `heating_rate`
`(nlay, ncol)` [K/s], and the underlying `solver`.

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
    return (;
        solver,
        lw_up = lw_flux_up(solver),
        lw_dn = lw_flux_dn(solver),
        lw_net = lw_flux_net(solver),
        sw_up = sw_flux_up(solver),
        sw_dn = sw_flux_dn(solver),
        sw_net = sw_flux_net(solver),
        net = net_flux(solver),
        heating_rate = heating_rate(solver),
    )
end
