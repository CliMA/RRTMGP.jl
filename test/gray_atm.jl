using Test
using RRTMGP.Device: array_type, array_device
using KernelAbstractions
using CUDA
using RRTMGP
using RRTMGP.AngularDiscretizations
using RRTMGP.Fluxes
using RRTMGP.RTE
using RRTMGP.GrayRTESolver
using RRTMGP.Optics
using RRTMGP.GrayUtils
using RRTMGP.AtmosphericStates
using RRTMGP.Sources
using RRTMGP.BCs

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
    opc::Symbol,
    ncol::Int,
    ::Type{FT},
    ::Type{I},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}
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
    inc_flux = nothing                # incoming flux
    max_threads = Int(256)            # maximum number of threads for KA kernels

    if ncol == 1
        lat = DA{FT}([0])             # latitude
    else
        lat = DA{FT}(range(FT(-π / 2), FT(π / 2), length = ncol)) # latitude
    end

    println("Running longwave test for gray atmosphere model - $opc; ncol = $ncol; DA = $DA")
    as = setup_gray_as_pr_grid(nlay, lat, p0, pe, param_set, DA)
    op = init_optical_props(opc, FT, DA, ncol, nlay)
    src_lw = source_func_longwave(FT, ncol, nlay, ngpt, opc, DA)
    bcs_lw = LwBCs(DA{FT,1}(sfc_emis), inc_flux)
    angle_disc = AngularDiscretization(opc, FT, n_gauss_angles, DA)
    fluxb_lw = nothing #FluxLW(ncol, nlay, FT, DA) 
    flux_lw = FluxLW(ncol, nlay, FT, DA)

    slv = Solver{
        FT,
        I,
        DA{FT,1},
        DA{FT,2},
        typeof(as),
        typeof(op),
        typeof(src_lw),
        Nothing,
        typeof(bcs_lw),
        Nothing,
        typeof(angle_disc),
        typeof(fluxb_lw),
        Nothing,
        typeof(flux_lw),
        Nothing,
    }(
        as,
        op,
        src_lw,
        nothing,
        bcs_lw,
        nothing,
        angle_disc,
        fluxb_lw,
        nothing,
        flux_lw,
        nothing,
    )

    flux = slv.flux_lw
    flux_up = slv.flux_lw.flux_up
    flux_dn = slv.flux_lw.flux_dn
    flux_net = slv.flux_lw.flux_net
    gray_as = slv.as
    t_lay = gray_as.t_lay
    p_lay = gray_as.p_lay
    t_lev = gray_as.t_lev
    p_lev = gray_as.p_lev
    optical_props = slv.op
    source = slv.src_lw
    #----------------------------------------------------
    hr_lay = DA{FT}(undef, nlay, ncol)
    T_ex_lev = DA{FT}(undef, nlev, ncol)
    flux_grad = DA{FT}(undef, nlay, ncol)
    flux_grad_err = FT(0)

    for i = 1:nsteps
        # calling the long wave gray radiation solver
        solve_lw!(slv, max_threads = max_threads)
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
end
#------------------------------------------------------------------------------

function gray_atmos_sw_test(
    opc::Symbol,
    ncol::Int,
    ::Type{FT},
    ::Type{I},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}
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
    n_gauss_angles = I(1)             # for non-scattering calculation
    sfc_emis = Array{FT}(undef, ncol) # surface emissivity
    sfc_emis .= FT(1.0)
    max_threads = Int(256)            # maximum number of threads for KA kernels

    toa_flux = Array{FT}(undef, ncol) # top of atmosphere flux
    sfc_alb_direct = Array{FT}(undef, ncol) # surface albedo (direct)
    sfc_alb_diffuse = Array{FT}(undef, ncol) # surface albedo (diffuse)
    inc_flux_diffuse = nothing
    zenith = Array{FT}(undef, ncol) # cosecant of solar zenith angle

    top_at_1 = false                          # Top-of-atmos at pt# 1 (true/false)

    if ncol == 1
        lat = DA{FT}([0])                     # latitude
    else
        lat = DA{FT}(range(FT(-π / 2), FT(π / 2), length = ncol)) # latitude
    end

    toa_flux .= FT(1407.679)
    sfc_alb_direct .= FT(0.1)
    sfc_alb_diffuse .= FT(0.1)
    zenith .= FT(1.66)               # corresponding to ~52.95 deg zenith angle    

    as = setup_gray_as_pr_grid(nlay, lat, p0, pe, param_set, DA) # init gray atmos state
    op = init_optical_props(opc, FT, DA, ncol, nlay) # init optical properties
    src_sw = source_func_shortwave(FT, ncol, nlay, opc, DA)
    bcs_sw = SwBCs(
        DA,
        FT,
        toa_flux,
        sfc_alb_direct,
        sfc_alb_diffuse,
        zenith,
        inc_flux_diffuse,
    )
    fluxb_sw = nothing
    flux_sw = FluxSWNoScat(ncol, nlay, FT, DA)
    println("Running shortwave test for gray atmosphere model - $(opc); ncol = $ncol; DA = $DA")

    slv = Solver{
        FT,
        I,
        DA{FT,1},
        DA{FT,2},
        typeof(as),
        typeof(op),
        Nothing,
        typeof(src_sw),
        Nothing,
        typeof(bcs_sw),
        Nothing,
        Nothing,
        typeof(fluxb_sw),
        Nothing,
        typeof(flux_sw),
    }(
        as,
        op,
        nothing,
        src_sw,
        nothing,
        bcs_sw,
        nothing,
        nothing,
        fluxb_sw,
        nothing,
        flux_sw,
    )
    solve_sw!(slv, max_threads = max_threads)

    # testing with exact solution
    ot_tot = sum(slv.op.τ[:, 1]) * zenith[1]
    exact = toa_flux[1] * zenith[1] * exp(-ot_tot)

    rel_toler = FT(0.001)
    rel_error = abs(slv.flux_sw.flux_dn_dir[1] - exact) / exact
    println("relative error = $rel_error")
    @test rel_error < rel_toler
end

if DA == CuArray
    @time gray_atmos_lw_equil(:OneScalar, Int(4096), Float64, Int, DA)
    @time gray_atmos_lw_equil(:TwoStream, Int(4096), Float64, Int, DA)
else
    @time gray_atmos_lw_equil(:OneScalar, Int(9), Float64, Int, DA)
    @time gray_atmos_lw_equil(:TwoStream, Int(9), Float64, Int, DA)
end

gray_atmos_sw_test(:OneScalar, Int(1), Float64, Int, DA)
gray_atmos_sw_test(:TwoStream, Int(1), Float64, Int, DA)
