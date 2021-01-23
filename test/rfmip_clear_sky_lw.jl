using RRTMGP
using RRTMGP.Device: array_type, array_device
using RRTMGP.Vmrs
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.AngularDiscretizations
using RRTMGP.RTE
using RRTMGP.GrayRTESolver

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using NCDatasets

include("read_rfmip_clear_sky.jl")
#---------------------------------------------------------------
function lw_rfmip(
    opc::Symbol,
    ::Type{FT},
    ::Type{I},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}

    homefolder = pwd()
    lw_file = joinpath(
        homefolder,
        "data",
        "rrtmgp",
        "data",
        "rrtmgp-data-lw-g256-2018-12-04.nc",
    )

    FTA1D = DA{FT,1}
    FTA2D = DA{FT,2}
    max_threads = Int(256)
    exp_no = 1
    n_gauss_angles = 1

    # reading longwave lookup data
    ds_lw = Dataset(lw_file, "r")
    lookup_lw, idx_gases = LookUpLW(ds_lw, I, FT, DA)
    close(ds_lw)
    # reading rfmip data to atmospheric state
    lw_input_file = joinpath(
        homefolder,
        "data",
        "examples",
        "rfmip-clear-sky",
        "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc",
    )

    ds_lw_in = Dataset(lw_input_file, "r")

    as = setup_rfmip_as(
        ds_lw_in,
        idx_gases,
        exp_no,
        lookup_lw,
        FT,
        DA,
        max_threads,
    )
    close(ds_lw_in)

    ncol, nlay, ngpt = as.ncol, as.nlay, lookup_lw.n_gpt
    nlev = nlay + 1
    op = init_optical_props(opc, FT, DA, ncol, nlay)                     # allocating optical properties object
    src_lw = source_func_longwave(FT, ncol, nlay, ngpt, opc, DA)        # allocating longwave source function object
    bcs_lw = LwBCs{FT,typeof(as.sfc_emis),Nothing}(as.sfc_emis, nothing) # setting up boundary conditions
    ang_disc = AngularDiscretization(opc, FT, n_gauss_angles, DA)        # initializing Angular discretization
    fluxb_lw = FluxLW(ncol, nlay, FT, DA)                                # flux storage for bandwise calculations
    flux_lw = FluxLW(ncol, nlay, FT, DA)                                 # longwave fluxes
    # initializing RTE solver
    slv = Solver{
        FT,
        I,
        FTA1D,
        FTA2D,
        typeof(as),
        typeof(op),
        typeof(src_lw),
        Nothing,
        typeof(bcs_lw),
        Nothing,
        typeof(ang_disc),
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
        ang_disc,
        fluxb_lw,
        nothing,
        flux_lw,
        nothing,
    )
    solve_lw!(slv, lookup_lw, max_threads = max_threads)

    for i = 1:10
        @time solve_lw!(slv, lookup_lw, max_threads = max_threads)
    end


    # reading comparison data

    flux_up_comp_file = joinpath(
        homefolder,
        "data",
        "examples",
        "rfmip-clear-sky",
        "rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
    )


    flux_dn_comp_file = joinpath(
        homefolder,
        "data",
        "examples",
        "rfmip-clear-sky",
        "rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc",
    )

    flip_ind = nlev:-1:1

    ds_flux_up = Dataset(flux_up_comp_file, "r")
    comp_flux_up = ds_flux_up["rlu"][:][flip_ind, :, exp_no]
    close(ds_flux_up)

    ds_flux_dn = Dataset(flux_dn_comp_file, "r")
    comp_flux_dn = ds_flux_dn["rld"][:][flip_ind, :, exp_no]
    close(ds_flux_dn)

    max_err_flux_up = maximum(abs.(slv.flux_lw.flux_up - comp_flux_up))
    max_err_flux_dn = maximum(abs.(slv.flux_lw.flux_dn - comp_flux_dn))
    println(" $opc; max_err_flux_up = $max_err_flux_up")
    println(" $opc; max_err_flux_dn = $max_err_flux_dn")
end


lw_rfmip(:OneScalar, Float64, Int, array_type())
lw_rfmip(:TwoSream, Float64, Int, array_type())
