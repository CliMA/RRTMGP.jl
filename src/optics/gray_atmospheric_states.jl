
"""
    AbstractGrayOpticalThickness

Abstract type for the gray-radiation optical-thickness parameterizations:
[`GrayOpticalThicknessSchneider2004`](@ref) and
[`GrayOpticalThicknessOGorman2008`](@ref).
"""
abstract type AbstractGrayOpticalThickness end

"""
    GrayOpticalThicknessSchneider2004{FT} <: AbstractGrayOpticalThickness
    GrayOpticalThicknessSchneider2004(FT; α = 3.5, te = 300, tt = 200, Δt = 60)

Parameters of the semi-grey optical-thickness profile of [schneider2004](@cite).
The vertical longwave optical depth scales with pressure as ``\\widehat{\\tau}(p) =
d_0(\\phi)\\,(p/p_0)^\\alpha``, so the exponent `α` controls how the absorber is
distributed with height (`α = 1` for a well-mixed absorber, larger `α` for one
concentrated near the surface). The optical thickness at the surface,
``d_0(\\phi)``, is set by the latitude-dependent radiative-equilibrium
temperature built from `te`, `tt`, and `Δt`.

# Fields
- `α`: Pressure exponent of the optical-depth profile.
- `te`: Global-mean surface temperature [K].
- `tt`: Temperature at the top of the atmosphere [K].
- `Δt`: Equator-to-pole surface temperature difference [K].
"""
struct GrayOpticalThicknessSchneider2004{FT} <: AbstractGrayOpticalThickness
    α::FT
    te::FT
    tt::FT
    Δt::FT
end
Adapt.@adapt_structure GrayOpticalThicknessSchneider2004

GrayOpticalThicknessSchneider2004(
    ::Type{FT};
    α = 3.5,
    te = 300,
    tt = 200,
    Δt = 60,
) where {FT} =
    GrayOpticalThicknessSchneider2004{FT}(FT(α), FT(te), FT(tt), FT(Δt))

"""
    GrayOpticalThicknessOGorman2008{FT} <: AbstractGrayOpticalThickness
    GrayOpticalThicknessOGorman2008(FT; α = 1.0, fₗ = 0.2, τₑ = 7.2, τₚ = 1.8, τ₀ = 0.22)

Parameters of the grey optical-thickness profile of [frierson2006](@cite) and
[ogorman2008](@cite). The
vertical longwave optical depth blends a linear and a quartic dependence on the
normalized pressure ``\\sigma = p/p_0``, ``\\widehat{\\tau} \\propto f_\\ell\\,\\sigma +
(1 - f_\\ell)\\,\\sigma^4``, so that it stays finite aloft while thickening toward
the surface, with an equator-to-pole contrast set by `τₑ` and `τₚ`. The shortwave
optical depth scales as ``\\sigma^2`` with amplitude `τ₀`. This is the default
optical-thickness model for [`solve_gray`](@ref RRTMGP.solve_gray).

# Fields
- `α`: Overall scaling factor for the longwave optical depth.
- `fₗ`: Weight of the linear (versus quartic) pressure dependence.
- `τₑ`: Longwave optical thickness at the equator.
- `τₚ`: Longwave optical thickness at the poles.
- `τ₀`: Shortwave optical thickness amplitude.
"""
struct GrayOpticalThicknessOGorman2008{FT} <: AbstractGrayOpticalThickness
    α::FT
    fₗ::FT
    τₑ::FT
    τₚ::FT
    τ₀::FT
end
Adapt.@adapt_structure GrayOpticalThicknessOGorman2008

GrayOpticalThicknessOGorman2008(
    ::Type{FT};
    α = 1.0,
    fₗ = 0.2,
    τₑ = 7.2,
    τₚ = 1.8,
    τ₀ = 0.22,
) where {FT <: AbstractFloat} =
    GrayOpticalThicknessOGorman2008{FT}(FT(α), FT(fₗ), FT(τₑ), FT(τₚ), FT(τ₀))


"""
    GrayAtmosphericState{FT, FTA1D, FTA2D, OTP} <: AbstractAtmosphericState

Atmospheric conditions, used to compute optical properties with the gray atmosphere
approximation.

# Fields
- `lat`: Latitude in degrees for each column `(ncol,)`.
- `p_lay`: Layer pressures [Pa, mb] `(nlay, ncol)`.
- `p_lev`: Level pressures [Pa, mb] `(nlay+1, ncol)`.
- `t_lay`: Layer temperatures [K] `(nlay, ncol)`.
- `t_lev`: Level temperatures [K] `(nlay+1, ncol)`.
- `z_lev`: Level altitudes [m] `(nlay+1, ncol)`.
- `t_sfc`: Surface temperatures [K] `(ncol)`.
- `otp`: Optical thickness parameters.
"""
struct GrayAtmosphericState{
    FT <: AbstractFloat,
    FTA1D <: AbstractArray{FT, 1},
    FTA2D <: AbstractArray{FT, 2},
    OTP <: AbstractGrayOpticalThickness,
} <: AbstractAtmosphericState
    lat::FTA1D
    p_lay::FTA2D
    p_lev::FTA2D
    t_lay::FTA2D
    t_lev::FTA2D
    z_lev::FTA2D
    t_sfc::FTA1D
    otp::OTP
end
Adapt.@adapt_structure GrayAtmosphericState
# Number of layers
@inline get_nlay(as::GrayAtmosphericState) = size(as.p_lay, 1)
# Number of columns
@inline get_ncol(as::GrayAtmosphericState) = size(as.p_lay, 2)
# Number of layers and columns
@inline get_dims(as::GrayAtmosphericState) = size(as.p_lay)

# view of layer pressures [Pa, mb]
@inline getview_p_lay(as::GrayAtmosphericState) = as.p_lay
@inline getview_p_lay(as::GrayAtmosphericState, gcol) =
    @inbounds view(as.p_lay, :, gcol)

# view of layer temperatures [K]
@inline getview_t_lay(as::GrayAtmosphericState) = as.t_lay
@inline getview_t_lay(as::GrayAtmosphericState, gcol) =
    @inbounds view(as.t_lay, :, gcol)

"""
    setup_gray_as_pr_grid(
        context::ClimaComms.AbstractCommsContext,
        nlay::Int,
        lat,
        p0,
        pe,
        otp::AbstractGrayOpticalThickness,
        param_set,
        ::Type{DA},
        step = "linear",
    )

Build a [`GrayAtmosphericState`](@ref) on a pressure grid from surface
pressure `p0` to top-of-atmosphere pressure `pe`, with the analytic
temperature profile of Schneider (2004), J. Atmos. Sci. 61 (12), 1317-1340
(doi: [10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2](https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2)).
One column is built per element of `lat`, with arrays of type `DA`.
"""
function setup_gray_as_pr_grid(
    context::ClimaComms.AbstractCommsContext,
    nlay::Int,
    lat::FTA1D,
    p0::FT,
    pe::FT,
    otp::AbstractGrayOpticalThickness,
    param_set::RP.ARP,
    ::Type{DA},
    step = "linear",
) where {FT <: AbstractFloat, FTA1D <: AbstractArray{FT, 1}, DA}
    nlev = Int(nlay + 1)
    ncol = length(lat)
    p_lay = DA{FT}(undef, nlay, ncol) # layer mean pressure
    p_lev = DA{FT}(undef, nlev, ncol) # level pressure
    t_lay = DA{FT}(undef, nlay, ncol) # layer mean temperature
    t_lev = DA{FT}(undef, nlev, ncol) # level temperature
    z_lev = DA{FT}(undef, nlev, ncol) # level altitude
    t_sfc = DA{FT}(undef, ncol)       # surface temperature
    efac = log(p0 / pe) / nlay       # multiplication factor for each pressure step
    d0 = DA{FT}(undef, ncol)       # optical depth (function of latitude)
    Δp = (p0 - pe) / nlay          # Δp for linear pressure distribution
    te = FT(300)                   # global mean surface temperature (K)
    tt = FT(200)                   # skin temp at top of atmosphere (K)
    Δt = FT(60)
    α = FT(3.5)                   # lapse rate of radiative equillibrium
    τ₀ = FT(0.22)
    r_d = RP.R_d(param_set)
    grav_ = RP.grav(param_set)
    args = (
        p_lev,
        p_lay,
        t_lev,
        t_lay,
        z_lev,
        t_sfc,
        lat,
        d0,
        efac,
        p0,
        pe,
        Δp,
        te,
        tt,
        Δt,
        α,
        τ₀,
        r_d,
        grav_,
        nlay,
    )
    device = ClimaComms.device(context)
    setup_gray_as_pr_grid!(device, ncol, args...)
    #------------------------------------------------
    return GrayAtmosphericState{
        eltype(t_sfc),
        typeof(t_sfc),
        typeof(p_lev),
        typeof(otp),
    }(
        lat,
        p_lay,
        p_lev,
        t_lay,
        t_lev,
        z_lev,
        t_sfc,
        otp,
    )
end

# This functions sets up a model temperature and pressure 
# distributions for a gray atmosphere based on an altitude grid
# with internal GLL point distribution within each cell
# see Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
# https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2

function setup_gray_as_pr_grid!(
    device::ClimaComms.AbstractCPUDevice,
    ncol,
    args...,
)
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            setup_gray_as_pr_grid_kernel!(args..., gcol)
        end
    end
end

function setup_gray_as_pr_grid_kernel!(
    p_lev::FTA2D,
    p_lay::FTA2D,
    t_lev::FTA2D,
    t_lay::FTA2D,
    z_lev::FTA2D,
    t_sfc::FTA1D,
    lat::FTA1D,
    d0::FTA1D,
    efac::FT,
    p0::FT,
    pe::FT,
    Δp::FT,
    te::FT,
    tt::FT,
    Δt::FT,
    α::FT,
    τ₀::FT,
    r_d::FT,
    grav_::FT,
    nlay::Int,
    gcol::Int,
) where {
    FT <: AbstractFloat,
    FTA1D <: AbstractArray{FT, 1},
    FTA2D <: AbstractArray{FT, 2},
}
    ts = te + Δt * (FT(1) / FT(3) - sin(lat[gcol] / FT(180) * FT(π))^2) # surface temp at a given latitude (K)
    d0[gcol] = FT((ts / tt)^FT(4) - FT(1)) # optical depth
    nlev = nlay + 1

    #---bot_at_1------------------------------
    p_lev[1, gcol] = p0
    t_lev[1, gcol] = tt * (FT(1) + d0[gcol] * (p_lev[1, gcol] / p0)^α)^FT(0.25)
    z_lev[1, gcol] = FT(0)

    @inbounds for ilay in 1:nlay
        #                if step == "linear"
        p_lev[ilay + 1, gcol] = p_lev[ilay, gcol] - Δp
        #                else
        #                    p_lev[ilay+1, gcol] = p_lev[ilay, gcol] * exp(-efac)
        #                end
        p_lay[ilay, gcol] =
            (p_lev[ilay, gcol] + p_lev[ilay + 1, gcol]) * FT(0.5)

        t_lev[ilay + 1, gcol] =
            tt * (FT(1) + d0[gcol] * (p_lev[ilay + 1, gcol] / p0)^α)^FT(0.25)
        t_lay[ilay, gcol] =
            tt * (FT(1) + d0[gcol] * (p_lay[ilay, gcol] / p0)^α)^FT(0.25)

        H = r_d * t_lay[ilay, gcol] / grav_
        Δz_lay = H * log(p_lev[ilay, gcol] / p_lev[ilay + 1, gcol])
        z_lev[ilay + 1, gcol] = Δz_lay + z_lev[ilay, gcol]
    end
    t_sfc[gcol] = t_lev[1, gcol]
    return nothing
end
