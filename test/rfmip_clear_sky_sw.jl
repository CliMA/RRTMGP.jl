using RRTMGP
using RRTMGP.Device: array_type, array_device
using RRTMGP.LookUpTables
using RRTMGP.Vmrs
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using NCDatasets

include("read_rfmip_clear_sky.jl")
#---------------------------------------------------------------
function sw_rfmip(
    opc::Symbol,
    ::Type{FT},
    ::Type{I},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}

    homefolder = pwd()
    sw_file = joinpath(
        homefolder,
        "data",
        "rrtmgp",
        "data",
        "rrtmgp-data-sw-g224-2018-12-04.nc",
    )

    FTA1D = DA{FT,1}
    FTA2D = DA{FT,2}
    max_threads = Int(256)
    exp_no = 1
    # reading shortwave lookup data
    ds_sw = Dataset(sw_file, "r")
    lookup_sw, idx_gases = LookUpSW(ds_sw, I, FT, DA)
    close(ds_sw)
    # reading rfmip data to atmospheric state
    sw_input_file = joinpath(
        homefolder,
        "data",
        "examples",
        "rfmip-clear-sky",
        "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc",
    )

    ds_sw_in = Dataset(sw_input_file, "r")

    as = setup_rfmip_as(
        ds_sw_in,
        idx_gases,
        exp_no,
        lookup_sw,
        FT,
        DA,
        max_threads,
    )
    close(ds_sw_in)

    ncol, nlay, ngpt = as.ncol, as.nlay, lookup_sw.n_gpt
    nlev = nlay + 1
    op = init_optical_props(opc, FT, DA, ncol, nlay)                # allocating optical properties object
    src_sw = source_func_shortwave(FT, ncol, nlay, opc, DA)        # allocating longwave source function object
    # setting up boundary conditions
    toa_flux = deepcopy(as.irrad)
    sfc_alb_direct = deepcopy(as.sfc_alb)
    zenith = deepcopy(as.zenith)
    inc_flux_diffuse = nothing
    sfc_alb_diffuse = nothing
    #    bcs_sw = SwBCs(DA, FT, toa_flux, sfc_alb_direct, 
    #                   sfc_alb_diffuse, zenith, inc_flux_diffuse)

    ang_disc = nothing        # initializing Angular discretization
    if opc == OneScalar
        fluxb_sw = FluxSWNoScat(ncol, nlay, FT, DA)                 # flux storage for bandwise calculations
        flux_sw = FluxSWNoScat(ncol, nlay, FT, DA)                 # longwave fluxes
    else
        fluxb_sw = FluxSW2Str(ncol, nlay, FT, DA)                 # flux storage for bandwise calculations
        flux_sw = FluxSW2Str(ncol, nlay, FT, DA)                 # longwave fluxes
    end

    #igpt = 1

    #compute_optical_props!(as, lookup_sw, op, igpt, max_threads = max_threads)
end

sw_rfmip(:OneScalar, Float64, Int, array_type())
