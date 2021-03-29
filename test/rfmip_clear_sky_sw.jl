using Test
using Pkg.Artifacts
using NCDatasets

using RRTMGP
using RRTMGP.Device: array_type, array_device
using RRTMGP.LookUpTables
using RRTMGP.Vmrs
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.Fluxes
using RRTMGP.BCs
using RRTMGP.RTE
using RRTMGP.RTESolver

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
# overriding some parameters to match with RRTMGP FORTRAN code
CLIMAParameters.Planet.grav(::EarthParameterSet) = 9.80665
CLIMAParameters.Planet.molmass_dryair(::EarthParameterSet) = 0.028964
CLIMAParameters.Planet.molmass_water(::EarthParameterSet) = 0.018016

include("reference_files.jl")
include("read_rfmip_clear_sky.jl")
#---------------------------------------------------------------
function sw_rfmip(
    opc::Symbol,
    ::Type{FT},
    ::Type{I},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}

    sw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :sw) # sw lookup tables
    sw_input_file = get_ref_filename(:atmos_state, :clearsky)      # clear-sky atmos state
    # reference data files for comparison
    flux_up_file = get_ref_filename(
        :comparison,
        :clearsky,
        λ = :sw,
        flux_up_dn = :flux_up,
        opc = :TwoStream,
    )
    flux_dn_file = get_ref_filename(
        :comparison,
        :clearsky,
        λ = :sw,
        flux_up_dn = :flux_dn,
        opc = :TwoStream,
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
    ds_sw_in = Dataset(sw_input_file, "r")

    (as, _, sfc_alb_direct, zenith, toa_flux) = setup_rfmip_as(
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
    op = init_optical_props(opc, FT, DA, ncol, nlay)               # allocating optical properties object
    src_sw = source_func_shortwave(FT, ncol, nlay, opc, DA)        # allocating longwave source function object

    # setting up boundary conditions
    inc_flux_diffuse = nothing
    sfc_alb_diffuse = nothing
    bcs_sw = SwBCs(
        DA,
        FT,
        zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )

    ang_disc = nothing        # initializing Angular discretization
    fluxb_sw = init_flux_sw(ncol, nlay, FT, DA, opc) # flux storage for bandwise calculations
    flux_sw = init_flux_sw(ncol, nlay, FT, DA, opc)  # shortwave fluxes for band calculations
    # initializing RTE solver
    slv = Solver{
        FT,
        I,
        FTA1D,
        FTA2D,
        typeof(as),
        typeof(op),
        Nothing,
        typeof(src_sw),
        Nothing,
        typeof(bcs_sw),
        typeof(ang_disc),
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
        ang_disc,
        nothing,
        fluxb_sw,
        nothing,
        flux_sw,
    )
    # block out coluumns with zenith > π/2
    usecol = BitArray(undef, ncol)
    usecol .= 1
    for i = 1:ncol
        if bcs_sw.zenith[i] > FT(π) / 2 - 2 * eps(FT)
            bcs_sw.toa_flux[i] = FT(0)
            usecol[i] = 0
        end
    end
    #--------------------------------------------------

    solve_sw!(slv, lookup_sw, max_threads = max_threads)
    for i = 1:10
        @time solve_sw!(slv, lookup_sw, max_threads = max_threads)
    end

    # reading comparison data

    flip_ind = nlev:-1:1

    ds_flux_up = Dataset(flux_up_file, "r")
    comp_flux_up = ds_flux_up["rsu"][:][flip_ind, :, exp_no]
    close(ds_flux_up)

    ds_flux_dn = Dataset(flux_dn_file, "r")
    comp_flux_dn = ds_flux_dn["rsd"][:][flip_ind, :, exp_no]
    close(ds_flux_dn)

    flux_up = Array(slv.flux_sw.flux_up)
    flux_dn = Array(slv.flux_sw.flux_dn)

    for i = 1:ncol
        if usecol[i] == 0
            flux_up[:, i] .= FT(0)
            flux_dn[:, i] .= FT(0)
        end
    end

    max_err_flux_up = maximum(abs.(flux_up .- comp_flux_up))
    max_err_flux_dn = maximum(abs.(flux_dn .- comp_flux_dn))

    println("=======================================")
    println("Clear-sky shortwave test, opc  = $opc")
    println("max_err_flux_up = $max_err_flux_up")
    println("max_err_flux_dn = $max_err_flux_dn")

    toler = FT(0.001)
    @test maximum(abs.(flux_up .- comp_flux_up)) ≤ toler
    @test maximum(abs.(flux_dn .- comp_flux_dn)) ≤ toler
end

sw_rfmip(:TwoStream, Float64, Int, array_type()) # two-stream solver should be used for the short-wave problem
#sw_rfmip(:OneScalar, Float64, Int, array_type()) # this only computes flux_dn_dir
