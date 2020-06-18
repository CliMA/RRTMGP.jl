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

using CLIMAParameters
import CLIMAParameters.Planet: molmass_dryair, molmass_water
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
# overriding CLIMAParameters as different precision is needed by RRTMGP
CLIMAParameters.Planet.molmass_dryair(::EarthParameterSet) = 0.028964
CLIMAParameters.Planet.molmass_water(::EarthParameterSet) = 0.018016

@static if haspkg("Plots")
    using Plots
    const export_plots = true
else
    const export_plots = false
end

#-----------
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# Only run the per-grid-point tests
pgp_only = false
include("data_set_files.jl")
datafolder = RRTMGP.data_folder_rrtmgp()
println("datafolder = $datafolder")
ds = dataset_dict(data_files_dict(datafolder, "lw"))

#-----------

include(joinpath("PostProcessing.jl"))
include(joinpath("..", "ReadInputData", "ReadInputs.jl"))
include(joinpath("..", "ReadInputData", "LoadCoefficients.jl"))

"""
Example program to demonstrate the calculation of longwave radiative fluxes in clear, aerosol-free skies.
  The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
  The large problem (1800 profiles) is divided into blocks

Program is invoked as rrtmgp_rfmip_lw [block_size input_file  coefficient_file upflux_file downflux_file]
  All arguments are optional but need to be specified in order.

RTE shortwave driver

RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  Here we're just reporting broadband fluxes

RRTMGP's gas optics class needs to be initialized with data read from a netCDF files

Optical properties of the atmosphere as array of values
   In the longwave we include only absorption optical depth (_1scl)
   Shortwave calculations use optical depth, single-scattering albedo, asymmetry parameter (_2str)
"""
function rfmip_clear_sky_lw_ga(ds, optical_props_constructor)

    FT = Float64
    I = Int
    deg_to_rad = acos(-FT(1)) / FT(180)
    angle_disc = GaussQuadrature(FT, 1)

    ncol, nlay, nexp = read_size(ds[:rfmip])
    #**********************************************************
    ncol, nexp = 1, 1
    #**********************************************************
    forcing_index = 1
    block_size = ncol * nexp #8

    # How big is the problem? Does it fit into blocks of the size we've specified?
    @assert 1 <= forcing_index <= 3
    #
    # Identify the set of gases used in the calculation based on the forcing index
    #   A gas might have a different name in the k-distribution than in the files
    #   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
    #
    kdist_gas_names = determine_gas_names(ds[:k_dist], forcing_index)
    # --------------------------------------------------
    #
    # Prepare data for use in rte+rrtmgp
    #
    #
    # Allocation on assignment within reading routines
    #
    p_lay_all, p_lev_all, t_lay_all, t_lev_all = read_and_block_pt(ds[:rfmip])
    #Restricting to col #1****************************************************
    p_lay_all = p_lay_all[1:1, :]
    p_lev_all = p_lev_all[1:1, :]
    t_lay_all = t_lay_all[1:1, :]
    t_lev_all = t_lev_all[1:1, :]
    #*************************************************************************
    # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
    #
    top_at_1 = p_lay_all[1, 1] < p_lay_all[1, nlay]
    #
    # Read the gas concentrations and surface properties
    #
    gas_conc = read_and_block_gases_ty(ds[:rfmip], kdist_gas_names)
    #Restricting to col #1****************************************************
    gas_conc = GasConcs{FT,I}(
        gas_conc.gas_names,
        gas_conc.concs[:, 1:1, :],
        ncol,
        nlay,
    )
    #*************************************************************************
    sfc_emis_all, t_sfc_all = read_and_block_lw_bc(ds[:rfmip])
    #Restricting to col #1****************************************************
    sfc_emis_all = sfc_emis_all[1:1]
    t_sfc_all = t_sfc_all[1:1]
    #*************************************************************************
    # Read k-distribution information:
    k_dist = load_and_init(ds[:k_dist], FT, gas_conc.gas_names)
    @assert source_is_internal(k_dist)

    nbnd = get_nband(k_dist.optical_props)
    ngpt = get_ngpt(k_dist.optical_props)
    #
    # RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
    #   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
    #   This introduces an error but shows input sanitizing.
    #
    if top_at_1
        p_lev_all[:, 1] .= get_press_min(k_dist.ref) + eps(FT)
    else
        p_lev_all[:, nlay+1] .= get_press_min(k_dist.ref) + eps(FT)
    end
    #
    # RTE will fail if passed solar zenith angles greater than 90 degree. We replace any with
    #   nighttime columns with a default solar zenith angle. We'll mask these out later, of
    #   course, but this gives us more work and so a better measure of timing.
    #

    #
    # Allocate space for output fluxes (accessed via pointers in FluxesBroadBand),
    #   gas optical properties, and source functions. The %alloc() routines carry along
    #   the spectral discretization from the k-distribution.
    #
    flux_up = Array{FT}(undef, block_size, nlay + 1)
    flux_dn = Array{FT}(undef, block_size, nlay + 1)

    sfc_emis_spec = Array{FT}(undef, nbnd, block_size)

    optical_props =
        optical_props_constructor(k_dist.optical_props, block_size, nlay, ngpt)
    source = SourceFuncLongWave(block_size, nlay, k_dist.optical_props)
    #
    # Loop over blocks
    #
    fluxes = FluxesBroadBand(FT, (size(flux_up, 1), size(flux_up, 2)))

    local as

    for icol = 1:block_size
        for ibnd = 1:nbnd
            sfc_emis_spec[ibnd, icol] = sfc_emis_all[icol]
        end
    end
    p_lay = p_lay_all[:, :]
    p_lev = p_lev_all[:, :]
    t_lay = t_lay_all[:, :]
    t_lev = t_lev_all[:, :]
    t_sfc = t_sfc_all[:]

    as = AtmosphericState(
        gas_conc,
        p_lay,
        p_lev,
        t_lay,
        t_lev,
        k_dist.ref,
        param_set,
        nothing,
        t_sfc,
    )

    fluxes.flux_up .= FT(0)
    fluxes.flux_dn .= FT(0)

    gas_optics!(k_dist, as, optical_props, source)

    bcs = LongwaveBCs(sfc_emis_spec)

    rte_lw!(fluxes, optical_props, as.mesh_orientation, bcs, source, angle_disc)

    flux_up[:, :] .= fluxes.flux_up
    flux_dn[:, :] .= fluxes.flux_dn

    if export_plots
        case = "clearsky_lw_" * string(optical_props_constructor)
        heating_rate, z =
            compute_heating_rate(fluxes.flux_up, fluxes.flux_dn, as)
        plot(
            heating_rate,
            z,
            title = "Clear sky longwave heating rates",
            xlabel = "heating rate",
            ylabel = "pressure",
        )
        out_dir = "output"
        mkpath(out_dir)
        savefig(joinpath(out_dir, case * ".png"))
    end

    # reshaping the flux_up and flux_dn arrays for comparison with Fortran code.
    flux_up = reshape_for_comparison(flux_up, nlay, ncol, nexp)
    flux_dn = reshape_for_comparison(flux_dn, nlay, ncol, nexp)

    rel_err_up, rel_err_dn = rel_err(flux_up, flux_dn, rlu_ref, rld_ref)

    println("******************************************************")
    println("rfmip_clear_sky_lw_1c -> $optical_props_constructor \n")
    println("Max rel_err_flux_up = $rel_err_up; Max rel_err_flux_dn = $rel_err_dn")
    println("******************************************************")

    if optical_props_constructor == TwoStream
        @test rel_err_up < FT(0.01)
        @test rel_err_dn < FT(0.15)
    else
        @test rel_err_up < FT(1e-6)
        @test rel_err_dn < FT(1e-6)
    end

    return nothing
end
#-----------------------------------------------------------
println("OneScalar, normal version")
rfmip_clear_sky_lw_ga(ds, OneScalar)
println("============================================")
println("TwoStream, normal version")
rfmip_clear_sky_lw_ga(ds, TwoStream)
println("============================================")
#-----------------------------------------------------------
