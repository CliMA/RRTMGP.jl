
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
} <: AbstractAtmosphericState{FT, FTA1D}
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
    "Number of layers."
    nlay::Int
    "Number of columns."
    ncol::Int
end
Adapt.@adapt_structure GrayAtmosphericState
#---------------------------------------------------------------
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
    if device isa ClimaComms.CUDADevice
        max_threads = 256
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        @cuda threads = (tx) blocks = (bx) setup_gray_as_pr_grid_CUDA!(ncol, args...)
    else
        @inbounds begin
            ClimaComms.@threaded device for gcol in 1:ncol
                setup_gray_as_pr_grid_kernel!(args..., gcol)
            end
        end
    end
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
        nlay,
        ncol,
    )
end

function setup_gray_as_pr_grid_CUDA!(ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        setup_gray_as_pr_grid_kernel!(args..., gcol)
    end
    return nothing
end
# This functions sets up a model temperature and pressure 
# distributions for a gray atmosphere based on an altitude grid
# with internal GLL point distribution within each cell
# see Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
# https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2

function setup_gray_as_alt_grid(
    lat::FTA1D,
    p0::FT,
    zend::FT,
    ncls::Int,
    poly_order::Int,
    param_set::RP.ARP,
    ::Type{DA},
) where {FT <: AbstractFloat, FTA1D <: AbstractArray{FT, 1}, DA}
    ncol = length(lat)
    nlev = Int(ncls + 1 + (poly_order - 1) * ncls)
    nlay = Int(nlev - 1)
    qm1 = poly_order + 1
    zst = FT(0)
    te = FT(300)                   # global mean surface temperature (K)
    tt = FT(200)                   # skin temp at top of atmosphere (K)
    Δt = FT(60)
    α = FT(3.5)           # lapse rate of radiative equillibrium
    τ₀ = FT(0.22)

    p_lay = DA{FT}(undef, nlay, ncol)
    p_lev = DA{FT}(undef, nlev, ncol)
    t_lay = DA{FT}(undef, nlay, ncol)
    t_lev = DA{FT}(undef, nlev, ncol)
    z_lev = DA{FT}(undef, nlev, ncol)
    t_sfc = DA{FT}(undef, ncol)       # surface temperature
    d0 = DA{FT}(undef, ncol)       # optical depth (function of latitude)

    pts, _ = GaussQuadrature.legendre(FT, poly_order + 1, GaussQuadrature.both)
    pts = (pts .+ FT(1)) ./ FT(2)
    Δz_cell = zend / ncls

    for icol in 1:ncol
        #---------bot_at_1---------------------------------------
        z_lev[1, icol] = zst
        z_lev[nlev, icol] = zend

        for icls in 1:ncls
            for j in 1:(poly_order - 1)
                z_lev[(icls - 1) * poly_order + 1 + j, icol] =
                    z_lev[(icls - 1) * poly_order + 1, icol] + pts[1 + j] * Δz_cell
            end
            z_lev[icls * poly_order + 1, icol] = zst + icls * Δz_cell
        end

        ts = te + Δt * (FT(1 / 3) - sin(lat[icol] / FT(180) * FT(π))^2) # surface temp at a given latitude (K)
        d0[icol] = FT((ts / tt)^4 - 1) # optical depth

        p_lev[1, icol] = p0
        t_lev[1, icol] = tt * (1 + d0[icol] * (p_lev[icol, 1] / p0)^α)^FT(0.25)

        for ilay in 1:nlay
            H = RP.R_d(param_set) * t_lev[ilay, icol] / RP.grav(param_set)
            Δz_lay = abs(z_lev[ilay + 1, icol] - z_lev[ilay, icol])
            p_lev[ilay + 1, icol] = p_lev[ilay, icol] * exp(-Δz_lay / H)
            t_lev[ilay + 1, icol] = tt * (1 + d0[icol] * (p_lev[ilay + 1, icol] / p0)^α)^FT(0.25)

            t_lay[ilay, icol] = FT(0.5) * (t_lev[ilay, icol] + t_lev[ilay + 1, icol])
            p_lay[ilay, icol] = FT(0.5) * (p_lev[ilay, icol] + p_lev[ilay + 1, icol])
        end
        #--------------------------------------------------------
        t_sfc[icol] = t_lev[1, icol]
    end
    return GrayAtmosphericState{FT, DA{FT, 1}, DA{FT, 2}}(
        p_lay,
        p_lev,
        t_lay,
        t_lev,
        z_lev,
        t_sfc,
        α,
        τ₀,
        d0,
        nlay,
        ncol,
    )
end

#-------------------------------------------------------------------------
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

    for ilay in 1:nlay
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

    #---------------------------------
end
#-------------------------------------------------------------------------
