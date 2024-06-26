
abstract type AbstractGrayOpticalThickness end

"""
    GrayOpticalThicknessSchneider2004{FT} <: AbstractGrayOpticalThickness

Optical thickness function parameters from Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
DOI: https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayOpticalThicknessSchneider2004{FT} <: AbstractGrayOpticalThickness
    "scaling height ratio"
    α::FT
    "global mean surface temperature (K)"
    te::FT
    "temp at top of atmosphere (K)"
    tt::FT
    "Δt (K)"
    Δt::FT
end
Adapt.@adapt_structure GrayOpticalThicknessSchneider2004

GrayOpticalThicknessSchneider2004(::Type{FT}) where {FT} =
    GrayOpticalThicknessSchneider2004{FT}(FT(3.5), FT(300), FT(200), FT(60))

"""
    GrayOpticalThicknessOGorman2008{FT} <: AbstractGrayOpticalThickness

Optical thickness function parameters from O'Gorman 2008, Journal of Climate Vol 21, Page(s): 3815–3832.
DOI: https://doi.org/10.1175/2007JCLI2065.1

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayOpticalThicknessOGorman2008{FT} <: AbstractGrayOpticalThickness
    "scaling factor"
    α::FT
    "fₗ"
    fₗ::FT
    "longwave optical thickness at equator"
    τₑ::FT
    "longwave optical thickness at poles"
    τₚ::FT
    "optical thickness for shortwave radiation"
    τ₀::FT
end
Adapt.@adapt_structure GrayOpticalThicknessOGorman2008

GrayOpticalThicknessOGorman2008(::Type{FT}) where {FT <: AbstractFloat} =
    GrayOpticalThicknessOGorman2008{FT}(FT(1.0), FT(0.2), FT(7.2), FT(1.8), FT(0.22))


"""
    GrayAtmosphericState{FT,FTA1D,FTA2D} <: 
        AbstractAtmosphericState{FT,FTA1D}

Atmospheric conditions, used to compute optical properties with the gray atmosphere approximation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayAtmosphericState{
    FT <: AbstractFloat,
    FTA1D <: AbstractArray{FT, 1},
    FTA2D <: AbstractArray{FT, 2},
    OTP <: AbstractGrayOpticalThickness,
} <: AbstractAtmosphericState
    "latitude, in degrees, for each column; `(ncol,)`"
    lat::FTA1D
    "Layer pressures `[Pa, mb]`; `(nlay, ncol)`"
    p_lay::FTA2D
    "Level pressures `[Pa, mb]`; `(nlay+1, ncol)`"
    p_lev::FTA2D
    "Layer temperatures `[K]`; `(nlay, ncol)`"
    t_lay::FTA2D
    "Level temperatures `[K]`; `(nlay+1,ncol)`"
    t_lev::FTA2D
    "Level Altitude `[m]`; `(nlay+1,ncol)`"
    z_lev::FTA2D
    "Surface temperatures `[K]`; `(ncol)`"
    t_sfc::FTA1D
    "optical thickness parameters"
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
@inline getview_p_lay(as::GrayAtmosphericState, gcol) = @inbounds view(as.p_lay, :, gcol)

# view of layer temperatures [K]
@inline getview_t_lay(as::GrayAtmosphericState) = as.t_lay
@inline getview_t_lay(as::GrayAtmosphericState, gcol) = @inbounds view(as.t_lay, :, gcol)

# This functions sets up a model temperature and pressure 
# distributions for a gray atmosphere based on a pressure grid
# see Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
# https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2

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
    args = (p_lev, p_lay, t_lev, t_lay, z_lev, t_sfc, lat, d0, efac, p0, pe, Δp, te, tt, Δt, α, τ₀, r_d, grav_, nlay)
    device = ClimaComms.device(context)
    setup_gray_as_pr_grid!(device, ncol, args...)
    #------------------------------------------------
    return GrayAtmosphericState{eltype(t_sfc), typeof(t_sfc), typeof(p_lev), typeof(otp)}(
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

function setup_gray_as_pr_grid!(device::ClimaComms.AbstractCPUDevice, ncol, args...)
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
) where {FT <: AbstractFloat, FTA1D <: AbstractArray{FT, 1}, FTA2D <: AbstractArray{FT, 2}}
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
        p_lay[ilay, gcol] = (p_lev[ilay, gcol] + p_lev[ilay + 1, gcol]) * FT(0.5)

        t_lev[ilay + 1, gcol] = tt * (FT(1) + d0[gcol] * (p_lev[ilay + 1, gcol] / p0)^α)^FT(0.25)
        t_lay[ilay, gcol] = tt * (FT(1) + d0[gcol] * (p_lay[ilay, gcol] / p0)^α)^FT(0.25)

        H = r_d * t_lay[ilay, gcol] / grav_
        Δz_lay = H * log(p_lev[ilay, gcol] / p_lev[ilay + 1, gcol])
        z_lev[ilay + 1, gcol] = Δz_lay + z_lev[ilay, gcol]
    end
    t_sfc[gcol] = t_lev[1, gcol]
    return nothing
end
