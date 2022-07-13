using Test
using Pkg.Artifacts
using NCDatasets

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
using RRTMGP.RTESolver

import CLIMAParameters as CP
include(joinpath(pkgdir(RRTMGP), "parameters", "create_parameters.jl"))
# overriding some parameters to match with RRTMGP FORTRAN code
overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
param_set = create_insolation_parameters(Float64, overrides)

include("reference_files.jl")
include("read_all_sky.jl")

DA = array_type()

function all_sky(
    ::Type{OPC},
    ::Type{FT},
    ::Type{I},
    ::Type{DA};
    use_lut::Bool = true,
) where {FT <: AbstractFloat, I <: Int, DA, OPC}
    opc = Symbol(OPC)
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    max_threads = Int(256)
    n_gauss_angles = 1
    ncol = 128 # repeats col#1 128 time per RRTMGP example

    lw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :lw)             # lw lookup tables for gas optics
    lw_cld_file = get_ref_filename(:lookup_tables, :cloudysky, λ = :lw)        # lw cloud lookup tables
    sw_file = get_ref_filename(:lookup_tables, :clearsky, λ = :sw)             # sw lookup tables for gas optics
    sw_cld_file = get_ref_filename(:lookup_tables, :cloudysky, λ = :sw)        # lw cloud lookup tables

    input_file = get_ref_filename(:atmos_state, :cloudysky)                    # all-sky atmos state
    if use_lut
        flux_file = get_ref_filename(:comparison, :cloudysky, lut_pade = :lut) # flux files for comparison (LUT)
    else
        flux_file = get_ref_filename(:comparison, :cloudysky, lut_pade = :pade) # flux files for comparison (pade)
    end

    #reading longwave gas optics lookup data
    ds_lw = Dataset(lw_file, "r")
    lookup_lw, idx_gases = LookUpLW(ds_lw, I, FT, DA)
    close(ds_lw)
    # reading longwave cloud lookup data
    ds_lw_cld = Dataset(lw_cld_file, "r")
    lookup_lw_cld = LookUpCld(ds_lw_cld, I, FT, DA, use_lut)
    close(ds_lw_cld)
    #reading shortwave gas optics lookup data
    ds_sw = Dataset(sw_file, "r")
    lookup_sw, idx_gases = LookUpSW(ds_sw, I, FT, DA)
    close(ds_sw)
    # reading longwave cloud lookup data
    ds_sw_cld = Dataset(sw_cld_file, "r")
    lookup_sw_cld = LookUpCld(ds_sw_cld, I, FT, DA, use_lut)
    close(ds_sw_cld)
    # reading input file 
    ds_in = Dataset(input_file, "r")
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, zenith, toa_flux, bot_at_1 = setup_allsky_as(
        ds_in,
        idx_gases,
        lookup_lw,
        lookup_sw,
        lookup_lw_cld,
        lookup_sw_cld,
        use_lut,
        ncol,
        FT,
        DA,
        max_threads,
    )
    close(ds_in)
    ncol, nlay = as.ncol, as.nlay
    nlev = nlay + 1
    #---reading comparison files -----------------------------------
    if bot_at_1
        flip_ind = 1:nlev
    else
        flip_ind = nlev:-1:1
    end

    # Setting up longwave problem---------------------------------------
    ngpt_lw = lookup_lw.n_gpt
    op = OPC(FT, ncol, nlay, DA) # allocating optical properties object
    src_lw = source_func_longwave(param_set, FT, ncol, nlay, opc, DA)   # allocating longwave source function object
    bcs_lw = LwBCs{FT, typeof(sfc_emis), Nothing}(sfc_emis, nothing)    # setting up boundary conditions
    fluxb_lw = FluxLW(ncol, nlay, FT, DA)                             # flux storage for bandwise calculations
    flux_lw = FluxLW(ncol, nlay, FT, DA)                              # longwave fluxes

    # Setting up shortwave problem---------------------------------------
    ngpt_sw = lookup_sw.n_gpt
    src_sw = source_func_shortwave(FT, ncol, nlay, opc, DA)        # allocating longwave source function object

    # setting up boundary conditions
    inc_flux_diffuse = nothing

    bcs_sw = SwBCs(zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)

    fluxb_sw = FluxSW(ncol, nlay, FT, DA) # flux storage for bandwise calculations
    flux_sw = FluxSW(ncol, nlay, FT, DA)  # shortwave fluxes for band calculations    
    #-------------------------------------------------------------------
    # initializing RTE solver
    slv = Solver(as, op, src_lw, src_sw, bcs_lw, bcs_sw, fluxb_lw, fluxb_sw, flux_lw, flux_sw)
    #------calling solvers
    println("calling longwave solver; ncol = $ncol")
    @time solve_lw!(slv, max_threads, lookup_lw, lookup_lw_cld)
    println("calling shortwave solver; ncol = $ncol")
    @time solve_sw!(slv, max_threads, lookup_sw, lookup_sw_cld)
    #-------------
    # comparison
    method = use_lut ? "Lookup Table Interpolation method" : "PADE method"
    ds_comp = Dataset(flux_file, "r")
    comp_flux_up_lw = Array(transpose(ds_comp["lw_flux_up"][:][1:ncol, flip_ind]))
    comp_flux_dn_lw = Array(transpose(ds_comp["lw_flux_dn"][:][1:ncol, flip_ind]))
    comp_flux_up_sw = Array(transpose(ds_comp["sw_flux_up"][:][1:ncol, flip_ind]))
    comp_flux_dn_sw = Array(transpose(ds_comp["sw_flux_dn"][:][1:ncol, flip_ind]))
    close(ds_comp)

    flux_up_lw = Array(slv.flux_lw.flux_up)
    flux_dn_lw = Array(slv.flux_lw.flux_dn)

    max_err_flux_up_lw = maximum(abs.(flux_up_lw .- comp_flux_up_lw))
    max_err_flux_dn_lw = maximum(abs.(flux_dn_lw .- comp_flux_dn_lw))

    println("=======================================")
    println("Cloudy-sky longwave test - $opc")
    println(method)
    println("max_err_flux_up_lw = $max_err_flux_up_lw")
    println("max_err_flux_dn_lw = $max_err_flux_dn_lw")

    flux_up_sw = Array(slv.flux_sw.flux_up)
    flux_dn_sw = Array(slv.flux_sw.flux_dn)
    flux_dn_dir_sw = Array(slv.flux_sw.flux_dn_dir)

    max_err_flux_up_sw = maximum(abs.(flux_up_sw .- comp_flux_up_sw))
    max_err_flux_dn_sw = maximum(abs.(flux_dn_sw .- comp_flux_dn_sw))

    println("Cloudy-sky shortwave test - $opc")
    println(method)
    println("max_err_flux_up_sw = $max_err_flux_up_sw")
    println("max_err_flux_dn_sw = $max_err_flux_dn_sw")
    println("=======================================")

    toler = FT(1e-5)

    @test max_err_flux_up_lw < toler
    @test max_err_flux_dn_lw < toler

    @test max_err_flux_up_sw < toler
    @test max_err_flux_dn_sw < toler
end

@testset "Cloudy (all-sky, Two-stream calculations using lookup table method" begin
    @time all_sky(TwoStream, Float64, Int, DA, use_lut = true)
end
@testset "Cloudy (all-sky), Two-stream calculations using Pade method" begin
    @time all_sky(TwoStream, Float64, Int, DA, use_lut = false)
end
