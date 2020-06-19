using Test
using RRTMGP
using NCDatasets
using ProgressMeter
using TimerOutputs
const to = TimerOutput()
using RRTMGP.OpticalProps
using RRTMGP.Utilities
using RRTMGP.GasOptics
using RRTMGP.GasConcentrations
using RRTMGP.RadiativeBoundaryConditions
using RRTMGP.RTESolver
using RRTMGP.Fluxes
using RRTMGP.AtmosphericStates
using RRTMGP.SourceFunctions
using RRTMGP.AngularDiscretizations
using RRTMGP.GrayAtmos
using RRTMGP.MeshOrientations

using CLIMAParameters
import CLIMAParameters.Planet: molmass_dryair, molmass_water, cp_d, grav
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
# overriding CLIMAParameters as different precision is needed by RRTMGP

"""
Example program to demonstrate the calculation of longwave radiative fluxes in a model gray atmosphere.
"""
function gray_atmos_lw(
    gray_rrtmgp::GrayRRTMGP{FT,I},
    gray_rte::GrayRTE{FT},
) where {FT<:AbstractFloat,I<:Int}

    angle_disc = GaussQuadrature(FT, 1)
    fluxes = gray_rte.fluxes
    hr_lay = gray_rte.hr_lay
    nlay = size(hr_lay, 2)

    as = gray_rrtmgp.as
    optical_props = gray_rrtmgp.op
    source = gray_rrtmgp.src
    bcs = gray_rrtmgp.bcs
    #**********************************************************
    gas_optics_gray_atmos!(as, optical_props, source)
    #*************************************************************************
    fill!(fluxes.flux_up, 0)
    fill!(fluxes.flux_dn, 0)
    fill!(fluxes.flux_net, 0)
    rte_lw_gray_solve!(
        fluxes,
        optical_props,
        as.mesh_orientation,
        bcs,
        source,
        angle_disc,
    )
    #*************************************************************************
    for i = 1:nlay
        hr_lay[1, i] =
            FT(grav(param_set)) *
            (fluxes.flux_net[1, i+1] - fluxes.flux_net[1, i]) /
            (as.p_lev[i+1] - as.p_lev[i]) / FT(cp_d(param_set))
    end
    #*************************************************************************
    return nothing
end

function gray_atmos_lw_equil(optical_props_constructor)
    FT = Float64
    I = Int
    lat = FT(0)     # latitude
    p0 = FT(100000) # surface pressure (Pa)
    pe = FT(9000) # TOA pressure (Pa)
    ncol, nlay, nbnd, ngpt = 1, 60, 1, 1
    nlev = nlay + 1
    tb = FT(320)
    tstep = 6.0 # in hours
    Δt = FT(60 * 60 * tstep)
    ndays = 365 * 16 #1 #16
    nsteps = ndays * (24 / tstep)
    temp_toler = 0.1        # in Kelvin
    flux_grad_toler = 1e-5  # tolerance for flux gradient

    sfc_emis = Array{FT}(undef, nbnd, ncol)
    sfc_emis[1, 1] = 1.0

    println("Running longwave test for gray atmosphere model")

    gray_rrtmgp =
        GrayRRTMGP(nlay, lat, p0, pe, optical_props_constructor, sfc_emis, "lw")
    gray_rte = GrayRTE(ncol, nlay, nlev, FT)

    fluxes, hr_lay = gray_rte.fluxes, gray_rte.hr_lay

    gray_as = gray_rrtmgp.as
    optical_props = gray_rrtmgp.op
    source = gray_rrtmgp.src

    # generate initial temp distrib
    gray_as.t_lev .= tb
    gray_as.t_lev[1, :] .= range(FT(200), tb, length = nlay + 1)
    for i = 1:nlay
        gray_as.t_lay[1, i] =
            0.5 * (gray_as.t_lev[1, i] + gray_as.t_lev[1, i+1])
    end

    T_ex_lev = zeros(FT, nlev)
    flux_grad = Array{FT}(undef, nlay)

    for i = 1:nsteps
        gray_atmos_lw(gray_rrtmgp, gray_rte)
        for j = 1:nlay
            gray_as.t_lay[1, j] = gray_as.t_lay[1, j] + Δt * hr_lay[1, j]
        end

        tlay_to_tlev!(
            gray_as.p_lay,
            gray_as.p_lev,
            gray_as.t_lay,
            gray_as.t_lev,
        )

        T_ex_lev[:] .=
            (
                (fluxes.flux_dn[1, :] .+ (fluxes.flux_net[1, :] ./ 2.0)) /
                Stefan()
            ) .^ 0.25
        flux_grad[:] .= fluxes.flux_net[1, 2:end] - fluxes.flux_net[1, 1:end-1]

        if maximum(abs.(flux_grad)) < flux_grad_toler
            println("Integration time = $(i/(24.0/tstep) / 365.0) years")
            break
        end
        #----------------------------------------------------
    end
    error = maximum(abs.(T_ex_lev[:] - gray_as.t_lev[1, :]))
    @test error < temp_toler
    #--------------------------------------------------------
end

gray_atmos_lw_equil(OneScalar)
