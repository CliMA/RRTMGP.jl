module GrayAtmosphericStates

using ..Device: array_type
using DocStringExtensions

using GaussQuadrature
using CLIMAParameters
using CLIMAParameters.Planet: grav, R_d

export GrayAtmosphericState, setup_gray_as_pr_grid, setup_gray_as_alt_grid

"""
    GrayAtmosphericState{FT}

Atmospheric conditions, used to compute optical properties with the gray atmosphere approximation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GrayAtmosphericState{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
    B<:Bool,
}
    "Layer pressures `[Pa, mb]`; (`ncol,nlay`)"
    p_lay::FTA2D
    "Level pressures `[Pa, mb]`; (`ncol,nlay+1`)"
    p_lev::FTA2D
    "Layer temperatures `[K]`; (`ncol,nlay`)"
    t_lay::FTA2D
    "Level temperatures `[K]`; (`ncol,nlay+1`)"
    t_lev::FTA2D
    "Level Altitude `[m]`; (`ncol,nlay+1`)"
    z_lev::FTA2D
    "Δz `[m]`; (`ncol,nlay`)"
    Δz_lay::FTA2D
    "Surface temperatures `[K]`; (`ncol`)"
    t_sfc::FTA1D
    "lapse rate"
    α::FT
    "optical thickness parameter"
    d0::FTA1D
    "Mesh orientation"
    top_at_1::B
    "Number of layers."
    nlay::I
    "Number of columns."
    ncol::I
end

# This functions sets up a model temperature and pressure 
# distributions for a gray atmosphere based on a pressure grid
# see Schneider 2004, J. Atmos. Sci. (2004) 61 (12): 1317–1340.
# https://doi.org/10.1175/1520-0469(2004)061<1317:TTATTS>2.0.CO;2

function setup_gray_as_pr_grid(
    nlay::Int,
    lat::FTA1D,
    p0::FT,
    pe::FT,
    param_set::AbstractEarthParameterSet,
    top_at_1::Bool,
    ::Type{DA},
    step = "linear",
) where {FT<:AbstractFloat,FTA1D<:AbstractArray{FT,1},DA}
    nlev = Int(nlay + 1)
    ncol = length(lat)
    p_lay = DA{FT}(undef, ncol, nlay) # layer mean pressure
    p_lev = DA{FT}(undef, ncol, nlev) # level pressure
    t_lay = DA{FT}(undef, ncol, nlay) # layer mean temperature
    t_lev = DA{FT}(undef, ncol, nlev) # level temperature
    z_lev = DA{FT}(undef, ncol, nlev) # level altitude
    Δz_lay = DA{FT}(undef, ncol, nlay) # layer thickness
    efac = log(p0 / pe) / nlay       # multiplication factor for each pressure step
    Δp = (p0 - pe) / nlay          # Δp for linear pressure distribution
    te = FT(300)                   # global mean surface temperature (K)
    tt = FT(200)                   # skin temp at top of atmosphere (K)
    Δt = FT(60)
    α = FT(3.5)                   # lapse rate of radiative equillibrium
    t_sfc = DA{FT}(undef, ncol)       # surface temperature
    d0 = DA{FT}(undef, ncol)       # optical depth (function of latitude)
    for icol = 1:ncol
        ts = te + Δt * (FT(1) / FT(3) - sin(lat[icol])^2) # surface temp at a given latitude (K)
        d0[icol] = FT((ts / tt)^FT(4) - FT(1)) # optical depth

        if top_at_1
            #---top_at_1------------------------------
            p_lev[icol, nlev] = p0
            t_lev[icol, nlev] =
                tt * (FT(1) + d0[icol] * (p_lev[icol, nlev] / p0)^α)^FT(0.25)
            z_lev[icol, nlev] = 0

            for ilay = nlay:-1:1
                if step == "linear"
                    p_lev[icol, ilay] = p_lev[icol, ilay+1] - Δp
                else
                    p_lev[icol, ilay] = p_lev[icol, ilay+1] * exp(-efac)
                end

                p_lay[icol, ilay] =
                    (p_lev[icol, ilay] + p_lev[icol, ilay+1]) * FT(0.5)

                t_lev[icol, ilay] =
                    tt *
                    (FT(1) + d0[icol] * (p_lev[icol, ilay] / p0)^α)^FT(0.25)
                t_lay[icol, ilay] =
                    tt *
                    (FT(1) + d0[icol] * (p_lay[icol, ilay] / p0)^α)^FT(0.25)

                H =
                    FT(R_d(param_set)) * t_lay[icol, ilay] /
                    FT((grav(param_set)))
                Δz_lay[icol, ilay] =
                    H * log(p_lev[icol, ilay+1] / p_lev[icol, ilay])
                z_lev[icol, ilay] = Δz_lay[icol, ilay] + z_lev[icol, ilay+1]
            end
        else
            #---bot_at_1------------------------------
            p_lev[icol, 1] = p0
            t_lev[icol, 1] =
                tt * (FT(1) + d0[icol] * (p_lev[icol, 1] / p0)^α)^FT(0.25)
            z_lev[icol, 1] = FT(0)

            for ilay = 1:nlay
                if step == "linear"
                    p_lev[icol, ilay+1] = p_lev[icol, ilay] - Δp
                else
                    p_lev[icol, ilay+1] = p_lev[icol, ilay] * exp(-efac)
                end
                p_lay[icol, ilay] =
                    (p_lev[icol, ilay] + p_lev[icol, ilay+1]) * FT(0.5)

                t_lev[icol, ilay+1] =
                    tt *
                    (FT(1) + d0[icol] * (p_lev[icol, ilay+1] / p0)^α)^FT(0.25)
                t_lay[icol, ilay] =
                    tt *
                    (FT(1) + d0[icol] * (p_lay[icol, ilay] / p0)^α)^FT(0.25)

                H =
                    FT(R_d(param_set)) * t_lay[icol, ilay] /
                    FT((grav(param_set)))
                Δz_lay[icol, ilay] =
                    H * log(p_lev[icol, ilay] / p_lev[icol, ilay+1])
                z_lev[icol, ilay+1] = Δz_lay[icol, ilay] + z_lev[icol, ilay]
            end
        end
        t_sfc[icol] = (top_at_1 ? t_lev[icol, end] : t_lev[icol, 1])
    end
    #---------------------------------
    @show top_at_1
    return GrayAtmosphericState{FT,DA{FT,1},DA{FT,2},Int,typeof(top_at_1)}(
        p_lay,
        p_lev,
        t_lay,
        t_lev,
        z_lev,
        Δz_lay,
        t_sfc,
        α,
        d0,
        top_at_1,
        nlay,
        ncol,
    )
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
    param_set::AbstractEarthParameterSet,
    top_at_1::Bool,
    ::Type{DA},
) where {FT<:AbstractFloat,FTA1D<:AbstractArray{FT,1},DA}
    ncol = length(lat)
    nlev = Int(ncls + 1 + (poly_order - 1) * ncls)
    nlay = Int(nlev - 1)
    qm1 = poly_order + 1
    zst = FT(0)
    te = FT(300)                   # global mean surface temperature (K)
    tt = FT(200)                   # skin temp at top of atmosphere (K)
    Δt = FT(60)
    α = FT(3.5)           # lapse rate of radiative equillibrium

    p_lay = DA{FT}(undef, ncol, nlay)
    p_lev = DA{FT}(undef, ncol, nlev)
    t_lay = DA{FT}(undef, ncol, nlay)
    t_lev = DA{FT}(undef, ncol, nlev)
    z_lev = DA{FT}(undef, ncol, nlev)
    Δz_lay = DA{FT}(undef, ncol, nlay)
    t_sfc = DA{FT}(undef, ncol)       # surface temperature
    d0 = DA{FT}(undef, ncol)       # optical depth (function of latitude)

    pts, _ = GaussQuadrature.legendre(FT, poly_order + 1, GaussQuadrature.both)
    pts = (pts .+ 1.0) ./ 2.0
    Δz_cell = zend / ncls

    for icol = 1:ncol
        if top_at_1
            #---------top_at_1--------------------------------------
            z_lev[icol, nlev] = zst
            z_lev[icol, 1] = zend

            for icls = 1:ncls
                for j = 1:poly_order-1
                    z_lev[icol, (icls - 1) * poly_order+1+j] =
                        z_lev[icol, (icls-1)*poly_order+1] - pts[1+j] * Δz_cell
                end
                z_lev[icol, icls*poly_order+1] = zend - icls * Δz_cell
            end

            for ilay = 1:nlay
                Δz_lay[icol, ilay] = z_lev[icol, ilay] - z_lev[icol, ilay+1]
            end

            ts = te + Δt * (FT(1 / 3) - sin(lat[icol])^2) # surface temp at a given latitude (K)
            d0[icol] = FT((ts / tt)^4 - 1) # optical depth

            p_lev[icol, nlev] = p0
            t_lev[icol, nlev] =
                tt * (1 + d0[icol] * (p_lev[icol, nlev] / p0)^α)^0.25

            for ilay = nlay:-1:1
                H =
                    FT(R_d(param_set)) * t_lev[icol, ilay+1] /
                    FT(grav(param_set))
                p_lev[icol, ilay] =
                    p_lev[icol, ilay+1] * exp(-Δz_lay[icol, ilay] / H)
                t_lev[icol, ilay] =
                    tt * (1 + d0[icol] * (p_lev[icol, ilay] / p0)^α)^FT(0.25)
                t_lay[icol, ilay] =
                    FT(0.5) * (t_lev[icol, ilay] + t_lev[icol, ilay+1])
                p_lay[icol, ilay] =
                    FT(0.5) * (p_lev[icol, ilay] + p_lev[icol, ilay+1])
            end
            #--------------------------------------------------------
        else
            #---------bot_at_1---------------------------------------
            z_lev[icol, 1] = zst
            z_lev[icol, nlev] = zend

            for icls = 1:ncls
                for j = 1:poly_order-1
                    z_lev[icol, (icls - 1) * poly_order+1+j] =
                        z_lev[icol, (icls-1)*poly_order+1] + pts[1+j] * Δz_cell
                end
                z_lev[icol, icls*poly_order+1] = zst + icls * Δz_cell
            end

            for ilay = 1:nlay
                Δz_lay[icol, ilay] = z_lev[icol, ilay+1] - z_lev[icol, ilay]
            end

            ts = te + Δt * (FT(1 / 3) - sin(lat[icol])^2) # surface temp at a given latitude (K)
            d0[icol] = FT((ts / tt)^4 - 1) # optical depth

            p_lev[icol, 1] = p0
            t_lev[icol, 1] = tt * (1 + d0[icol] * (p_lev[icol, 1] / p0)^α)^0.25

            for ilay = 1:nlay
                H = FT(R_d(param_set)) * t_lev[icol, ilay] / FT(grav(param_set))

                p_lev[icol, ilay+1] =
                    p_lev[icol, ilay] * exp(-Δz_lay[icol, ilay] / H)
                t_lev[icol, ilay+1] =
                    tt * (1 + d0[icol] * (p_lev[icol, ilay+1] / p0)^α)^FT(0.25)

                t_lay[icol, ilay] =
                    FT(0.5) * (t_lev[icol, ilay] + t_lev[icol, ilay+1])
                p_lay[icol, ilay] =
                    FT(0.5) * (p_lev[icol, ilay] + p_lev[icol, ilay+1])
            end
            #--------------------------------------------------------
        end

        t_sfc[icol] = (top_at_1 ? t_lev[icol, end] : t_lev[icol, 1])
    end
    @show top_at_1
    return GrayAtmosphericState{FT,DA{FT,1},DA{FT,2},Int,typeof(top_at_1)}(
        p_lay,
        p_lev,
        t_lay,
        t_lev,
        z_lev,
        Δz_lay,
        t_sfc,
        α,
        d0,
        top_at_1,
        nlay,
        ncol,
    )
end


end
