using Test
using RRTMGP.Device: array_type, array_device
using KernelAbstractions
using CUDA
using RRTMGP
using RRTMGP.GrayFluxes
using RRTMGP.GrayAtmos
using RRTMGP.GrayRTESolver
using RRTMGP.GrayOptics
using RRTMGP.GrayUtils

using CLIMAParameters
import CLIMAParameters.Planet: cp_d, grav, R_d
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
# overriding CLIMAParameters as different precision is needed by RRTMGP

#using Plots

DA = array_type()

"""
Example program to demonstrate the calculation of longwave radiative fluxes in a model gray atmosphere.
"""
function gray_atmos_lw_equil(
    ::Type{OPC},
    ncol::Int,
    ::Type{FT},
    ::Type{I},
    ::Type{DA},
) where {OPC<:AbstractGrayOpticalProps,FT<:AbstractFloat,I<:Int,DA}
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
    max_threads = Int(256)            # maximum number of threads for KA kernels

    if ncol == 1
        lat = DA{FT}([0])             # latitude
    else
        lat = DA{FT}(range(FT(-π / 2), FT(π / 2), length = ncol)) # latitude
    end

    println("Running longwave test for gray atmosphere model - $(OPC); ncol = $ncol; DA = $DA")
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

    flux = gray_rrtmgp.flux
    flux_up = gray_rrtmgp.flux.flux_up
    flux_dn = gray_rrtmgp.flux.flux_dn
    flux_net = gray_rrtmgp.flux.flux_net
    gray_as = gray_rrtmgp.as
    t_lay = gray_as.t_lay
    p_lay = gray_as.p_lay
    t_lev = gray_as.t_lev
    p_lev = gray_as.p_lev
    optical_props = gray_rrtmgp.op
    source = gray_rrtmgp.src
    #----------------------------------------------------
    hr_lay = DA{FT}(undef, nlay, ncol)
    T_ex_lev = DA{FT}(undef, nlev, ncol)
    flux_grad = DA{FT}(undef, nlay, ncol)
    flux_grad_err = FT(0)

    for i = 1:nsteps
        # calling the long wave gray radiation solver
        gray_atmos_lw!(gray_rrtmgp, max_threads = max_threads)
        # computing heating rate
        compute_gray_heating_rate!(gray_as, flux, hr_lay, param_set)
        # updating t_lay and t_lev based on heating rate
        update_profile_lw!(
            gray_as,
            flux,
            hr_lay,
            flux_grad,
            T_ex_lev,
            Δt,
            nlay,
            nlev,
            ncol,
        )
        #----------------------------------------
        flux_grad_err = maximum(flux_grad)
        if flux_grad_err < flux_grad_toler
            break
        end
        #----------------------------------------------------
    end

    t_error = maximum(abs.(T_ex_lev .- gray_as.t_lev))
    println("*************************************************")
    println("Integration time = $(FT(nsteps)/FT(24.0/tstep) / FT(365.0)) years")
    println("t_error = $(t_error); flux_grad_err = $(flux_grad_err)")
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
    gray_atmos_lw!(gray_rrtmgp_pr_grd)
    #----altitude grid
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
    ncol = size(gray_rrtmgp_alt_grd.as.p_lay, 2)
    nlay = size(gray_rrtmgp_alt_grd.as.p_lay, 1)
    nlev = nlay + 1

    gray_atmos_lw!(gray_rrtmgp_alt_grd)
    #--------------------------------------------------
    for icol = 1:ncol
        lat_str = "lat_" * string(lat[icol] * 180.0 / π) * "_"
        lat_tit = "lat = " * string(lat[icol] * 180.0 / π) * " deg"

        case = "flux_up_" * string(optical_props_constructor) * "_" * lat_str
        plot(
            gray_rrtmgp_pr_grd.as.p_lev[:, icol],
            gray_rrtmgp_pr_grd.flux.flux_up[:, icol],
            label = "pressure grid",
        )
        plot!(
            gray_rrtmgp_alt_grd.as.p_lev[:, icol],
            gray_rrtmgp_alt_grd.flux.flux_up[:, icol],
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

        case = "flux_dn_" * string(optical_props_constructor) * "_" * lat_str
        plot(
            gray_rrtmgp_pr_grd.as.p_lev[:, icol],
            gray_rrtmgp_pr_grd.flux.flux_dn[:, icol],
            label = "pressure grid",
        )
        plot!(
            gray_rrtmgp_alt_grd.as.p_lev[:, icol],
            gray_rrtmgp_alt_grd.flux.flux_dn[:, icol],
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
    return nothing
end
#------------------------------------------------------------------------------
if DA == CuArray
    @time gray_atmos_lw_equil(GrayOneScalar, Int(4096), Float64, Int, DA)
    @time gray_atmos_lw_equil(GrayTwoStream, Int(4096), Float64, Int, DA)
else
    @time gray_atmos_lw_equil(GrayOneScalar, Int(9), Float64, Int, DA)
    @time gray_atmos_lw_equil(GrayTwoStream, Int(9), Float64, Int, DA)
end
#=
# for visual verification
@time gray_atmos_lw_comparison(GrayOneScalar, Int(3), Float64, Int, DA)
@time gray_atmos_lw_comparison(GrayTwoStream, Int(3), Float64, Int, DA)
=#
