#####
##### Idealized standard-atmosphere profiles for the Layer-3 standalone path
#####

"""
    AtmosphereProfile

A self-contained, host-side description of a clear-sky atmospheric column (or
`ncol` identical columns), as produced by [`standard_atmosphere`](@ref) and
consumed by [`solve`](@ref). All arrays are plain `Array{FT}` on the host;
[`solve`](@ref) moves them to the compute device.

# Fields
- `p_lay`, `t_lay`: layer-center pressures [Pa] and temperatures [K], `(nlay, ncol)`.
- `p_lev`, `t_lev`: level pressures [Pa] and temperatures [K], `(nlay + 1, ncol)`.
- `z_lev`: level altitudes [m], `(nlay + 1, ncol)`.
- `t_sfc`: surface temperature [K], `(ncol,)`.
- `lat`: latitude [degrees], `(ncol,)`.
- `vmr_h2o`, `vmr_o3`: water-vapor and ozone volume mixing ratios, `(nlay, ncol)`.
- `well_mixed_vmr`: `Dict` of global-mean volume mixing ratios for the
  well-mixed gases, keyed by the RRTMGP gas names (see `gas_names_sw`); gases
  not listed are taken as zero.
"""
struct AtmosphereProfile{A1, A2, WM}
    p_lay::A2
    p_lev::A2
    t_lay::A2
    t_lev::A2
    z_lev::A2
    t_sfc::A1
    lat::A1
    vmr_h2o::A2
    vmr_o3::A2
    well_mixed_vmr::WM
end

# Idealized per-kind parameters: surface temperature [K], tropopause height
# [m], tropospheric lapse rate [K/m], stratospheric inverse lapse rate [K/m],
# surface and stratospheric-floor water-vapor vmr, default latitude [degrees].
# The values are idealized climatological choices (inspired by the AFGL
# reference-atmosphere climatology, but analytic rather than tabulated).
const _STANDARD_ATMOSPHERES = Dict(
    :tropical => (;
        t_sfc = 300.0,
        z_trop = 17.0e3,
        Γ_trop = 6.5e-3,
        Γ_strat = 2.2e-3,
        vmr_h2o_sfc = 2.3e-2,
        lat = 0.0,
    ),
    :midlatitude_summer => (;
        t_sfc = 294.0,
        z_trop = 13.0e3,
        Γ_trop = 6.5e-3,
        Γ_strat = 2.0e-3,
        vmr_h2o_sfc = 1.4e-2,
        lat = 45.0,
    ),
    :subarctic_winter => (;
        t_sfc = 257.0,
        z_trop = 9.0e3,
        Γ_trop = 5.0e-3,
        Γ_strat = 1.5e-3,
        vmr_h2o_sfc = 1.6e-3,
        lat = 65.0,
    ),
)

# Two-segment temperature profile: constant tropospheric lapse Γ_trop up to
# z_trop, then linear warming Γ_strat above (an idealized stratosphere).
function _standard_T(z, prm)
    z ≤ prm.z_trop && return prm.t_sfc - prm.Γ_trop * z
    return (prm.t_sfc - prm.Γ_trop * prm.z_trop) + prm.Γ_strat * (z - prm.z_trop)
end

# Hydrostatic pressure for the two-segment profile (exact closed forms for a
# constant-lapse-rate ideal-gas atmosphere: p ∝ T^(g/(R_d·Γ)), decreasing
# along −Γ and increasing along +Γ).
function _standard_p(z, prm, p_sfc, grav, R_d)
    T_trop = prm.t_sfc - prm.Γ_trop * prm.z_trop
    if z ≤ prm.z_trop
        return p_sfc * (_standard_T(z, prm) / prm.t_sfc)^(grav / (R_d * prm.Γ_trop))
    end
    p_trop = p_sfc * (T_trop / prm.t_sfc)^(grav / (R_d * prm.Γ_trop))
    return p_trop * (_standard_T(z, prm) / T_trop)^(-grav / (R_d * prm.Γ_strat))
end

# Idealized water vapor: exponential decay with a 2 km scale height above the
# surface value, floored at a stratospheric background of 4 ppmv.
_standard_vmr_h2o(z, prm) = max(prm.vmr_h2o_sfc * exp(-z / 2.0e3), 4.0e-6)

# Idealized ozone layer: log-pressure Gaussian peaking at ~7.5 ppmv near
# p = 1200 Pa (≈30 km), over a 30 ppbv tropospheric background.
_standard_vmr_o3(p) = 3.0e-8 + 7.5e-6 * exp(-(log(p / 1.2e3))^2 / (2 * 1.2^2))

"""
    standard_atmosphere(FT; kind = :midlatitude_summer, nlay = 60, ncol = 1,
                        z_top = 45.0e3, p_sfc = 101325.0,
                        params = default_parameters(FT))

Build an idealized clear-sky [`AtmosphereProfile`](@ref) on `nlay` layers with
levels uniformly spaced in altitude from the surface to `z_top` [m]:

- temperature: a constant tropospheric lapse rate up to an idealized
  tropopause, linear warming above (per-`kind` parameters);
- pressure: the exact hydrostatic profile for that temperature structure;
- water vapor: exponential decay (2 km scale height) to a 4 ppmv
  stratospheric floor; ozone: a log-pressure Gaussian layer peaking near
  30 km; well-mixed gases at present-day global means (CO₂ 420 ppmv,
  CH₄ 1.9 ppmv, N₂O 0.34 ppmv, CO 0.1 ppmv, O₂ 0.209, N₂ 0.781).

`kind` selects the idealized climatology: `:tropical`,
`:midlatitude_summer` (default), or `:subarctic_winter`. The profiles are
analytic and idealized — made for teaching and testing, not for reproducing
tabulated reference atmospheres. All `ncol` columns are identical.
"""
function standard_atmosphere(
    ::Type{FT};
    kind::Symbol = :midlatitude_summer,
    nlay::Int = 60,
    ncol::Int = 1,
    z_top = 45.0e3,
    p_sfc = 101325.0,
    params = default_parameters(FT),
) where {FT <: AbstractFloat}
    haskey(_STANDARD_ATMOSPHERES, kind) || error(
        "unknown standard-atmosphere kind `$(repr(kind))`; available kinds: " *
        "$(sort(collect(keys(_STANDARD_ATMOSPHERES))))",
    )
    prm = _STANDARD_ATMOSPHERES[kind]
    grav = Float64(RP.grav(params))
    R_d = Float64(RP.R_d(params))

    nlev = nlay + 1
    z_lev1 = range(0.0, Float64(z_top); length = nlev)
    z_lay1 = (z_lev1[1:(end - 1)] .+ z_lev1[2:end]) ./ 2

    p_lay1 = _standard_p.(z_lay1, Ref(prm), Float64(p_sfc), grav, R_d)
    p_lev1 = _standard_p.(z_lev1, Ref(prm), Float64(p_sfc), grav, R_d)

    tocols(v) = FT.(repeat(reshape(collect(v), :, 1), 1, ncol))
    return AtmosphereProfile(
        tocols(p_lay1),
        tocols(p_lev1),
        tocols(_standard_T.(z_lay1, Ref(prm))),
        tocols(_standard_T.(z_lev1, Ref(prm))),
        tocols(z_lev1),
        fill(FT(prm.t_sfc), ncol),
        fill(FT(prm.lat), ncol),
        tocols(_standard_vmr_h2o.(z_lay1, Ref(prm))),
        tocols(_standard_vmr_o3.(p_lay1)),
        Dict(
            "co2" => FT(420.0e-6),
            "ch4" => FT(1.9e-6),
            "n2o" => FT(0.34e-6),
            "co" => FT(1.0e-7),
            "o2" => FT(0.209),
            "n2" => FT(0.781),
        ),
    )
end
