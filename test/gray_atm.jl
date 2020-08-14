using Test
using CUDA
using KernelAbstractions
using RRTMGP.Device: array_type
using RRTMGP
using RRTMGP.GrayFluxes
using RRTMGP.GrayAtmos
using RRTMGP.GrayRTESolver
using RRTMGP.GrayOptics

using CLIMAParameters
import CLIMAParameters.Planet: cp_d, grav, R_d
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
# overriding CLIMAParameters as different precision is needed by RRTMGP

using Plots

DA = array_type()

"""
Example program to demonstrate the calculation of longwave radiative fluxes in a model gray atmosphere.
"""
function gray_atmos_lw(
    gray_rrtmgp::GrayRRTMGP{FT,I,FTA1D,FTA2D,FTA3D,B},
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    B<:Bool,
}
    # setting references
    flux_up = gray_rrtmgp.flux.flux_up   # upward flux
    flux_dn = gray_rrtmgp.flux.flux_dn   # downward flux
    flux_net = gray_rrtmgp.flux.flux_net # net flux
    hr_lay = gray_rrtmgp.flux.hr_lay     # heating rate
    nlay = size(hr_lay, 2)               # number of layers per column

    as = gray_rrtmgp.as            # gray atmospheric state
    optical_props = gray_rrtmgp.op # optical properties
    source = gray_rrtmgp.src       # Planck sources

    # computing optical properties
    gas_optics_gray_atmos!(as, optical_props, source)

    # solving radiative transfer equation
    if typeof(gray_rrtmgp.op) == GrayOneScalar{FT,FTA2D}
        rte_lw_noscat_gray_solve!(gray_rrtmgp) # no-scattering solver
    else
        rte_lw_2stream_gray_solve!(gray_rrtmgp) # 2-stream solver
    end

    # computing heating rate
    p_lev = gray_rrtmgp.as.p_lev
    for icol = 1:as.ncol
        for ilay = 1:nlay
            hr_lay[icol, ilay] =
                FT(grav(param_set)) *
                (flux_net[icol, ilay+1] - flux_net[icol, ilay]) /
                (p_lev[icol, ilay+1] - p_lev[icol, ilay]) / FT(cp_d(param_set))
        end
    end

    return nothing
end

function gray_atmos_lw_equil(
    ::Type{OPC},
    ::Type{FT},
    ::Type{I},
    ::Type{DA},
) where {OPC<:AbstractGrayOpticalProps,FT<:AbstractFloat,I<:Int,DA}
    ncol = 9                          # number of columns
    nlay = 60                         # number of layers
    p0 = FT(100000)                   # surface pressure (Pa)
    pe = FT(9000)                     # TOA pressure (Pa)
    nbnd, ngpt = 1, 1                 # # of nbands/g-points (=1 for gray radiation)
    nlev = nlay + 1                   # # of layers
    tb = FT(320)                      # surface temperature
    tstep = 6.0                       # timestep in hours
    Δt = FT(60 * 60 * tstep)          # timestep in seconds
    ndays = 365 * 40                  # # of simulation days
    nsteps = ndays * (24 / tstep)     # number of timesteps
    temp_toler = FT(0.1)              # tolerance for temperature (Kelvin)
    flux_grad_toler = FT(1e-5)        # tolerance for flux gradient
    n_gauss_angles = I(1)             # for non-scattering calculation
    sfc_emis = Array{FT}(undef, ncol) # surface emissivity
    sfc_emis .= FT(1.0)
    top_at_1 = false                  # Top-of-atmos at pt# 1 (true/false)

    if ncol == 1
        lat = DA{FT}([0])             # latitude
    else
        lat = DA{FT}(range(FT(-π / 2), FT(π / 2), length = ncol)) # latitude
    end

    println("Running longwave test for gray atmosphere model - $(OPC)")
    gray_flux = GrayFlux(ncol, nlay, nlev, FT, DA)

    gray_rrtmgp = GrayRRTMGP(
        nlay,
        lat,
        p0,
        pe,
        OPC,
        sfc_emis,
        gray_flux,
        n_gauss_angles,
        param_set,
        top_at_1,
        "lw",
        DA,
    )

    flux_up = gray_rrtmgp.flux.flux_up
    flux_dn = gray_rrtmgp.flux.flux_dn
    flux_net = gray_rrtmgp.flux.flux_net
    hr_lay = gray_rrtmgp.flux.hr_lay
    gray_as = gray_rrtmgp.as
    t_lay = gray_as.t_lay
    p_lay = gray_as.p_lay
    t_lev = gray_as.t_lev
    p_lev = gray_as.p_lev
    optical_props = gray_rrtmgp.op
    source = gray_rrtmgp.src
    #----------------------------------------------------
    # generate initial temp distrib
    gray_as.t_lev .= tb
    for icol = 1:ncol
        gray_as.t_lev[icol, :] .= range(FT(200), tb, length = nlay + 1)
        for ilay = 1:nlay
            t_lay[icol, ilay] = 0.5 * (t_lev[icol, ilay] + t_lev[icol, ilay+1])
        end
    end

    T_ex_lev = Array{FT}(undef, ncol, nlev)
    flux_grad = Array{FT}(undef, ncol, nlay)
    for i = 1:nsteps
        gray_atmos_lw(gray_rrtmgp)
        # updating t_lay based on heating rate
        for icol = 1:ncol, ilay = 1:nlay
            t_lay[icol, ilay] += Δt * hr_lay[icol, ilay]
        end

        # computing t_lev from t_lay
        tlay_to_tlev!(p_lay, p_lev, t_lay, t_lev)

        for icol = 1:ncol, ilev = 1:nlev
            T_ex_lev[icol, ilev] =
                (
                    (flux_dn[icol, ilev] + (flux_net[icol, ilev] / FT(2))) /
                    FT(Stefan())
                )^FT(0.25)
        end
        for icol = 1:ncol, ilay = 1:nlay
            flux_grad[icol, ilay] =
                abs(flux_net[icol, ilay+1] - flux_net[icol, ilay])
        end
        if maximum(flux_grad) < flux_grad_toler
            println("step # = $i")
            println("Integration time = $(FT(i)/FT(24.0/tstep) / FT(365.0)) years")
            println("maximum(flux_grad) = $(maximum(flux_grad, dims=2))")
            break
        end
        #----------------------------------------------------
    end
    t_error = maximum(abs.(T_ex_lev .- gray_as.t_lev), dims = 2)
    println("*************************************************")
    if maximum(flux_grad) > flux_grad_toler
        println("Integration time = $(FT(nsteps)/FT(24.0/tstep) / FT(365.0)) years")
    end
    println("maximum(flux_grad) = $(maximum(flux_grad, dims=2))")
    println("*************************************************")
    println("t_error = $(t_error)")
    println("*************************************************")

    @test maximum(t_error) < temp_toler
    #--------------------------------------------------------
    return nothing
end

function gray_atmos_lw_comparison(
    optical_props_constructor,
    ncol::Int,
    ::Type{FT},
    ::Type{I},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}
    p0 = FT(100000)                  # surface pressure (Pa)
    pe = FT(9000)                    # TOA pressure (Pa)
    nbnd, ngpt = 1, 1
    tb = FT(320)
    tstep = FT(6.0)                  # time step in hours
    Δt = FT(60 * 60 * tstep)         # time step in seconds
    ndays = 365 * 16                 # total integration time in days
    nsteps = ndays * (24 / tstep)    # total # of time steps
    zend = FT(15000)                 # upper limit for altitude (m)
    ncls = 10
    poly_ord = 5                     # polynomial order for altitude grid
    top_at_1 = false                 # is top of atmos at pt# 1?

    n_gauss_angles = I(1) # for non-scattering calculation
    sfc_emis = Array{FT}(undef, ncol)
    sfc_emis .= 1.0

    println("Running longwave comparison, for pressure grid vs altitude grid, test for gray atmosphere model")
    if ncol == 1
        lat = DA{FT}([0]) # latitude
    else
        lat = DA{FT}(range(FT(-π / 2), FT(π / 2), length = ncol)) # latitude
    end
    #----pressure grid
    nlay = 60
    nlev = nlay + 1
    gray_flux_pr_grd = GrayFlux(ncol, nlay, nlev, FT, DA)
    gray_rrtmgp_pr_grd = GrayRRTMGP(
        nlay,
        lat,
        p0,
        pe,
        optical_props_constructor,
        sfc_emis,
        gray_flux_pr_grd,
        n_gauss_angles,
        param_set,
        top_at_1,
        "lw",
        DA,
    )
    @show gray_rrtmgp_pr_grd.as.top_at_1
    @time gray_atmos_lw(gray_rrtmgp_pr_grd)
    #----altitude grid
    #=
        nlev = ncls + 1 + (poly_ord - 1) * ncls
        nlay = nlev - 1
        gray_flux_alt_grd = GrayFlux(ncol, nlay, nlev, FT, DA)
        gray_rrtmgp_alt_grd = GrayRRTMGP(
            lat,
            p0,
            zend,
            ncls,
            poly_ord,
            optical_props_constructor,
            sfc_emis,
            gray_flux_alt_grd,
            n_gauss_angles,
            param_set,
            top_at_1,
            "lw",
            DA,
        )
        ncol = size(gray_rrtmgp_alt_grd.as.p_lay, 1)
        nlay = size(gray_rrtmgp_alt_grd.as.p_lay, 2)
        nlev = nlay + 1

        gray_atmos_lw(gray_rrtmgp_alt_grd)
        #--------------------------------------------------
        for icol = 1:ncol
            lat_str = "lat_" * string(lat[icol] * 180.0 / π) * "_"
            lat_tit = "lat = " * string(lat[icol] * 180.0 / π) * " deg"

            case = "testing_flux_up_comparison_" * lat_str
            plot(
                gray_rrtmgp_pr_grd.as.p_lev[icol, :],
                gray_rrtmgp_pr_grd.flux.flux_up[icol, :],
                label = "pressure grid",
            )
            plot!(
                gray_rrtmgp_alt_grd.as.p_lev[icol, :],
                gray_rrtmgp_alt_grd.flux.flux_up[icol, :],
                label = "DG grid",
                markershapes = [:circle],
                title = "p lev vs flux_up" * lat_tit,
                xlabel = "p lev (Pa)",
                ylabel = "flux (K)",
                legend = :topleft,
            )
            out_dir = "output"
            mkpath(out_dir)
            savefig(joinpath(out_dir, case * ".png"))

            case = "testing_flux_dn_comparison_" * lat_str
            plot(
                gray_rrtmgp_pr_grd.as.p_lev[icol, :],
                gray_rrtmgp_pr_grd.flux.flux_dn[icol, :],
                label = "pressure grid",
            )
            plot!(
                gray_rrtmgp_alt_grd.as.p_lev[icol, :],
                gray_rrtmgp_alt_grd.flux.flux_dn[icol, :],
                label = "DG grid",
                markershapes = [:circle],
                title = "p lev vs flux_dn" * lat_tit,
                xlabel = "p lev (Pa)",
                ylabel = "flux (K)",
                legend = :topleft,
            )
            out_dir = "output"
            mkpath(out_dir)
            savefig(joinpath(out_dir, case * ".png"))
        end
    =#

    for icol in (1,)  #ncol
        println("********************")
        @show lat[icol] * 180 / π
        println("********************")
        println("------flux_up comparison- $(optical_props_constructor)------------------------")
        @show gray_rrtmgp_pr_grd.flux.flux_up[icol, :]
        #        @show gray_rrtmgp_alt_grd.flux.flux_up[icol, :]
        println("------flux_dn comparison--$(optical_props_constructor)-----------------------")
        @show gray_rrtmgp_pr_grd.flux.flux_dn[icol, :]
        #        @show gray_rrtmgp_alt_grd.flux.flux_dn[icol, :]
    end
    return nothing
end

#------------------------------------------------------------------------------
#@time gray_atmos_lw_equil(GrayOneScalar, Float64, Int, DA)
#@time gray_atmos_lw_equil(GrayTwoStream, Float64, Int, DA)
#for i = 1:2
#    @time gray_atmos_lw_equil(GrayOneScalar, Float64, Int, DA)
#end
#println("------------------------------")
#for i = 1:2
#    @time gray_atmos_lw_equil(GrayTwoStream, Float64, Int, DA)
#end

#println("CUDA.has_cuda_gpu = $(CUDA.has_cuda_gpu())")
@time gray_atmos_lw_comparison(GrayOneScalar, Int(512), Float64, Int, DA)
#gray_atmos_lw_comparison(GrayTwoStream, Int(9), Float64, Int, DA)
