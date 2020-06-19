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
using Printf

using CLIMAParameters
import CLIMAParameters.Planet: molmass_dryair, molmass_water, cp_d, grav
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
# overriding CLIMAParameters as different precision is needed by RRTMGP

using Plots

"""
Example program to demonstrate the calculation of longwave radiative fluxes in a model gray atmosphere.
"""
function gray_atmos_lw(
    gray_rrtmgp::GrayRRTMGP{FT,I},
    angle_disc::AbstractAngularDiscretization{FT,I},
) where {FT<:AbstractFloat,I<:Int}

    fluxes = gray_rrtmgp.rte.fluxes
    hr_lay = gray_rrtmgp.rte.hr_lay
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
    if typeof(gray_rrtmgp.op) == OneScalar{FT,I}
        rte_lw_gray_solve!(
            fluxes,
            optical_props,
            as.mesh_orientation,
            bcs,
            source,
            angle_disc,
        )
    else
        rte_lw_2stream_gray_solve!(
            fluxes,
            optical_props,
            as.mesh_orientation,
            bcs,
            source,
            angle_disc,
        )
    end
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

    println("Running longwave test for gray atmosphere model - $(optical_props_constructor)")

    gray_rte = GrayRTE(ncol, nlay, nlev, FT)
    gray_rrtmgp = GrayRRTMGP(
        nlay,
        lat,
        p0,
        pe,
        optical_props_constructor,
        sfc_emis,
        gray_rte,
        "lw",
    )

    fluxes, hr_lay = gray_rrtmgp.rte.fluxes, gray_rrtmgp.rte.hr_lay
    gray_as = gray_rrtmgp.as
    optical_props = gray_rrtmgp.op
    source = gray_rrtmgp.src
    angle_disc = GaussQuadrature(FT, 1)
    #----------------------------------------------------
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
        gray_atmos_lw(gray_rrtmgp, angle_disc)
        for j = 1:nlay
            gray_as.t_lay[1, j] = gray_as.t_lay[1, j] + Δt * hr_lay[1, j]
        end

        tlay_to_tlev!(
            gray_as.p_lay,
            gray_as.p_lev,
            gray_as.t_lay,
            gray_as.t_lev,
        )
        for ilev = 1:nlev
            T_ex_lev[ilev] =
                (
                    (
                        fluxes.flux_dn[1, ilev] +
                        (fluxes.flux_net[1, ilev] / 2.0)
                    ) / Stefan()
                )^0.25
        end
        for ilay = 1:nlay
            flux_grad[ilay] =
                abs(fluxes.flux_net[1, ilay+1] - fluxes.flux_net[1, ilay])
        end

        if maximum(abs.(flux_grad)) < flux_grad_toler
            println("Integration time = $(i/(24.0/tstep) / 365.0) years")
            break
        end
        #----------------------------------------------------
    end
    error = maximum(abs.(T_ex_lev[:] .- gray_as.t_lev[1, :]))

    @show error
    @test error < temp_toler
    #--------------------------------------------------------
    return nothing
end

function gray_atmos_lw_comparison(optical_props_constructor)
    FT = Float64
    I = Int
    lat = FT(0)     # latitude
    p0 = FT(100000) # surface pressure (Pa)
    pe = FT(9000) # TOA pressure (Pa)
    ncol, nbnd, ngpt = 1, 1, 1
    tb = FT(320)
    tstep = 6.0 # in hours
    Δt = FT(60 * 60 * tstep)
    ndays = 365 * 16 #1 #16
    nsteps = ndays * (24 / tstep)
    zend = FT(15000)
    ncls = 10
    poly_ord = 5
    angle_disc = GaussQuadrature(FT, 1)

    sfc_emis = Array{FT}(undef, nbnd, ncol)
    sfc_emis[1, 1] = 1.0

    println("Running longwave comparison, for pressure grid vs altitude grid, test for gray atmosphere model")
    #----pressure grid
    nlay = 60
    nlev = nlay + 1
    gray_rte_pr_grd = GrayRTE(ncol, nlay, nlev, FT)
    gray_rrtmgp_pr_grd = GrayRRTMGP(
        nlay,
        lat,
        p0,
        pe,
        optical_props_constructor,
        sfc_emis,
        gray_rte_pr_grd,
        "lw",
    )

    gray_atmos_lw(gray_rrtmgp_pr_grd, angle_disc)

    #----altitude grid
    nlev = ncls + 1 + (poly_ord - 1) * ncls
    nlay = nlev - 1
    gray_rte_alt_grd = GrayRTE(ncol, nlay, nlev, FT)
    gray_rrtmgp_alt_grd = GrayRRTMGP(
        lat,
        p0,
        zend,
        ncls,
        poly_ord,
        optical_props_constructor,
        sfc_emis,
        gray_rte_alt_grd,
        "lw",
    )
    ncol = size(gray_rrtmgp_alt_grd.as.p_lay, 1)
    nlay = size(gray_rrtmgp_alt_grd.as.p_lay, 2)
    nlev = nlay + 1

    gray_atmos_lw(gray_rrtmgp_alt_grd, angle_disc)

    #--------------------------------------------------
    case = "testing_flux_up_comparison"
    plot(
        gray_rrtmgp_pr_grd.as.p_lev[:],
        gray_rte_pr_grd.fluxes.flux_up[:],
        label = "pressure grid",
    )
    plot!(
        gray_rrtmgp_alt_grd.as.p_lev[:],
        gray_rte_alt_grd.fluxes.flux_up[:],
        label = "DG grid",
        markershapes = [:circle],
        title = "p lev vs flux_up",
        xlabel = "p lev (Pa)",
        ylabel = "flux (K)",
        legend = :topleft,
    )
    out_dir = "output"
    mkpath(out_dir)
    savefig(joinpath(out_dir, case * ".png"))

    case = "testing_flux_dn_comparison"
    plot(
        gray_rrtmgp_pr_grd.as.p_lev[:],
        gray_rte_pr_grd.fluxes.flux_dn[:],
        label = "pressure grid",
    )
    plot!(
        gray_rrtmgp_alt_grd.as.p_lev[:],
        gray_rte_alt_grd.fluxes.flux_dn[:],
        label = "DG grid",
        markershapes = [:circle],
        title = "p lev vs flux_dn",
        xlabel = "p lev (Pa)",
        ylabel = "flux (K)",
        legend = :topleft,
    )
    out_dir = "output"
    mkpath(out_dir)
    savefig(joinpath(out_dir, case * ".png"))

    @show gray_rte_pr_grd.fluxes.flux_up[:]
    @show gray_rte_alt_grd.fluxes.flux_up[:]

end
#------------------------------------------------------------------------------
function gray_atmos_sw(
    gray_rrtmgp::GrayRRTMGP{FT,I},
    gray_rte::GrayRTE{FT},
    angle_disc::AbstractAngularDiscretization{FT,I},
) where {FT<:AbstractFloat,I<:Int}

    return nothing
end

function gray_atmos_sw_test(optical_props_constructor)
    FT = Float64
    I = Int
    lat = FT(0)     # latitude
    p0 = FT(100000) # surface pressure (Pa)
    pe = FT(9000) # TOA pressure (Pa)
    ncol, nlay, nbnd, ngpt = 1, 60, 1, 1
    nlev = nlay + 1
    tb = FT(320)


    return nothing
end

#------------------------------------------------------------------------------
for i = 1:10
    @time gray_atmos_lw_equil(OneScalar)
end
println("------------------------------")
for i = 1:10
    @time gray_atmos_lw_equil(TwoStream)
end

#gray_atmos_lw_comparison(OneScalar)
#gray_atmos_sw_test(TwoStream)
