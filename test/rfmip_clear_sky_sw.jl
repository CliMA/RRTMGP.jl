using Test
using Pkg.Artifacts
using NCDatasets
import ClimaComms

using RRTMGP
using RRTMGP.LookUpTables
using RRTMGP.Vmrs
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.Fluxes
using RRTMGP.BCs
using RRTMGP.RTE
using RRTMGP.RTESolver
import RRTMGP.Parameters.RRTMGPParameters
import ClimaParams as CP
# overriding some parameters to match with RRTMGP FORTRAN code

include("reference_files.jl")
include("read_rfmip_clear_sky.jl")
#---------------------------------------------------------------
function sw_rfmip(context, ::Type{OPC}, ::Type{SRC}, ::Type{VMR}, ::Type{FT}) where {FT <: AbstractFloat, OPC, SRC, VMR}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)

    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    opc = Symbol(OPC)
    sw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :sw) # sw lookup tables
    sw_input_file = get_ref_filename(:atmos_state, :clearsky)      # clear-sky atmos state
    # reference data files for comparison
    flux_up_file = get_ref_filename(:comparison, :clearsky, λ = :sw, flux_up_dn = :flux_up, opc = :TwoStream)
    flux_dn_file = get_ref_filename(:comparison, :clearsky, λ = :sw, flux_up_dn = :flux_dn, opc = :TwoStream)

    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    max_threads = Int(256)
    exp_no = 1

    # reading shortwave lookup data
    ds_sw = Dataset(sw_file, "r")
    lookup_sw, idx_gases = LookUpSW(ds_sw, FT, DA)
    close(ds_sw)
    # reading rfmip data to atmospheric state
    ds_sw_in = Dataset(sw_input_file, "r")

    (as, _, sfc_alb_direct, cos_zenith, toa_flux, usecol) =
        setup_rfmip_as(context, ds_sw_in, idx_gases, exp_no, lookup_sw, FT, VMR, max_threads, param_set)
    close(ds_sw_in)

    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    op = OPC(FT, ncol, nlay, DA)            # allocating optical properties object
    #src_sw = SRC(FT, DA, nlay, ncol)        # allocating longwave source function object
    src_sw = source_func_shortwave(FT, ncol, nlay, opc, DA)        # allocating longwave source function object

    # setting up boundary conditions
    inc_flux_diffuse = nothing
    sfc_alb_diffuse = FTA2D(deepcopy(sfc_alb_direct))
    bcs_sw = SwBCs{FT, FTA1D, Nothing, FTA2D}(cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    fluxb_sw = FluxSW(ncol, nlay, FT, DA) # flux storage for bandwise calculations
    flux_sw = FluxSW(ncol, nlay, FT, DA)  # shortwave fluxes for band calculations

    # initializing RTE solver
    slv = Solver(context, as, op, nothing, src_sw, nothing, bcs_sw, nothing, fluxb_sw, nothing, flux_sw)
    #--------------------------------------------------
    solve_sw!(slv, max_threads, lookup_sw)

    # reading comparison data
    flip_ind = nlev:-1:1

    ds_flux_up = Dataset(flux_up_file, "r")
    comp_flux_up = Array(ds_flux_up["rsu"])[flip_ind, :, exp_no]
    close(ds_flux_up)

    ds_flux_dn = Dataset(flux_dn_file, "r")
    comp_flux_dn = Array(ds_flux_dn["rsd"])[flip_ind, :, exp_no]
    close(ds_flux_dn)

    flux_up = Array(slv.flux_sw.flux_up)
    flux_dn = Array(slv.flux_sw.flux_dn)

    for i in 1:ncol
        if usecol[i] == 0
            flux_up[:, i] .= FT(0)
            flux_dn[:, i] .= FT(0)
        end
    end

    max_err_flux_up = maximum(abs.(flux_up .- comp_flux_up))
    max_err_flux_dn = maximum(abs.(flux_dn .- comp_flux_dn))

    color2 = :cyan

    printstyled(
        "Stand-alone clear-sky shortwave test with ncol = $ncol, nlev = $nlev, OPC = $opc, FT = $FT\n",
        color = color2,
    )
    printstyled("device = $device\n\n", color = color2)
    println("L∞ error in flux_up = $max_err_flux_up")
    println("L∞ error in flux_dn = $max_err_flux_dn\n")


    toler = FT(0.001)
    @test maximum(abs.(flux_up .- comp_flux_up)) ≤ toler
    @test maximum(abs.(flux_dn .- comp_flux_dn)) ≤ toler
    return nothing
end

sw_rfmip(ClimaComms.context(), TwoStream, SourceSW2Str, VmrGM, Float64) # two-stream solver should be used for the short-wave problem
#sw_rfmip(OneScalar, VmrGM, Float64, Int, array_type()) # this only computes flux_dn_dir
