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

include(joinpath("PostProcessing.jl"))
include(joinpath("..", "ReadInputData", "ReadInputs.jl"))
include(joinpath("..", "ReadInputData", "LoadCoefficients.jl"))

"""
Example program to demonstrate the calculation of shortwave radiative fluxes in clear, aerosol-free skies.
  The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
  The large problem (1800 profiles) is divided into blocks

Program is invoked as rrtmgp_rfmip_sw [block_size input_file  coefficient_file upflux_file downflux_file]
  All arguments are optional but need to be specified in order.

RTE shortwave driver

RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  Here we're just reporting broadband fluxes

RRTMGP's gas optics class needs to be initialized with data read from a netCDF files

modules for reading and writing files


"""
function rfmip_clear_sky_sw(ds, optical_props_constructor)

    FT = Float64
    I = Int
    deg_to_rad = acos(-FT(1)) / FT(180)

    ncol, nlay, nexp = read_size(ds[:rfmip])

    forcing_index = 1
    block_size = ncol * nexp #8

    # How big is the problem? Does it fit into blocks of the size we've specified?
    @assert mod(ncol * nexp, block_size) == 0 # number of columns must fit evenly into blocks
    @assert !(forcing_index < 1 || forcing_index > 3)

    #
    # Identify the set of gases used in the calculation based on the forcing index
    #   A gas might have a different name in the k-distribution than in the files
    #   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
    #
    kdist_gas_names = @timeit to "determine_gas_names" determine_gas_names(
        ds[:k_dist],
        forcing_index,
    )

    # --------------------------------------------------
    #
    # Prepare data for use in rte+rrtmgp
    #
    #
    # Allocation on assignment within reading routines
    #
    p_lay,
    p_lev,
    t_lay,
    t_lev = @timeit to "read_and_block_pt" read_and_block_pt(ds[:rfmip])
    #
    # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
    #

    top_at_1 = p_lay[1, 1] < p_lay[1, nlay]
    #
    # Read the gas concentrations and surface properties
    #
    gas_conc = @timeit to "read_and_block_gases_ty" read_and_block_gases_ty(
        ds[:rfmip],
        kdist_gas_names,
    )
    surface_albedo, total_solar_irradiance, solar_zenith_angle =
        @timeit to "read_and_block_sw_bc" read_and_block_sw_bc(ds[:rfmip],)

    # Read k-distribution information:
    k_dist = @timeit to "load_and_init" load_and_init(
        ds[:k_dist],
        FT,
        gas_conc.gas_names,
    )
    @assert source_is_external(k_dist)

    nbnd = get_nband(k_dist.optical_props)
    ngpt = get_ngpt(k_dist.optical_props)

    #
    # RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
    #   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
    #   This introduces an error but shows input sanitizing.
    #
    if top_at_1
        p_lev[:, 1] .= get_press_min(k_dist.ref) + eps(FT)
    else
        p_lev[:, nlay+1] .= get_press_min(k_dist.ref) + eps(FT)
    end

    toa_flux = zeros(FT, block_size, get_ngpt(k_dist.optical_props))
    def_tsi = zeros(FT, block_size)
    usecol = Array{Bool}(undef, block_size)
    #
    # RTE will fail if passed solar zenith angles greater than 90 degree. We replace any with
    #   nighttime columns with a default solar zenith angle. We'll mask these out later, of
    #   course, but this gives us more work and so a better measure of timing.
    #
    usecol[1:block_size] .=
        solar_zenith_angle[1:block_size] .< FT(90) - FT(2) * eps(FT(90))

    #
    # Allocate space for output fluxes (accessed via pointers in FluxesBroadBand),
    #   gas optical properties, and source functions. The %alloc() routines carry along
    #   the spectral discretization from the k-distribution.
    #
    flux_up = zeros(FT, block_size, nlay + 1)
    flux_dn = zeros(FT, block_size, nlay + 1)

    μ_0 = zeros(FT, block_size)
    sfc_alb_spec = zeros(FT, nbnd, block_size)
    optical_props =
        optical_props_constructor(k_dist.optical_props, block_size, nlay, ngpt)

    #
    # Loop over blocks
    #
    fluxes = FluxesBroadBand(FT, (size(flux_up, 1), size(flux_up, 2)), true)

    local as

    as = AtmosphericState(
        gas_conc,
        p_lay,
        p_lev,
        t_lay,
        nothing,
        k_dist.ref,
        param_set,
        nothing,
        nothing,
    )

    # Compute the optical properties of the atmosphere and the Planck source functions
    #    from pressures, temperatures, and gas concentrations...
    fluxes.flux_up .= FT(0)
    fluxes.flux_dn .= FT(0)

    @timeit to "gas_optics!" gas_optics!(k_dist, as, optical_props, true)

    check_extent(toa_flux, (as.ncol, ngpt), "toa_flux")
    toa_flux .= repeat(k_dist.solar_src', as.ncol)
    # Boundary conditions
    #   (This is partly to show how to keep work on GPUs using OpenACC in a host application)
    # What's the total solar irradiance assumed by RRTMGP?
    #
    def_tsi[1:block_size] = sum(toa_flux, dims = 2)
    #
    # Normalize incoming solar flux to match RFMIP specification
    #
    for igpt = 1:ngpt
        for icol = 1:block_size
            toa_flux[icol, igpt] =
                toa_flux[icol, igpt] * total_solar_irradiance[icol] /
                def_tsi[icol]
        end
    end
    #
    # Expand the spectrally-constant surface albedo to a per-band albedo for each column
    #
    for icol = 1:block_size
        for ibnd = 1:nbnd
            sfc_alb_spec[ibnd, icol] = surface_albedo[icol]
        end
    end
    #
    # Cosine of the solar zenith angle
    #
    for icol = 1:block_size
        μ_0[icol] =
            usecol[icol] ? cos(solar_zenith_angle[icol] * deg_to_rad) : FT(1)
    end

    #
    # ... and compute the spectrally-resolved fluxes, providing reduced values
    #    via FluxesBroadBand
    #
    bcs = ShortwaveBCs(toa_flux, sfc_alb_spec, sfc_alb_spec)

    @timeit to "rte_sw!" rte_sw!(
        fluxes,
        optical_props,
        as.mesh_orientation,
        bcs,
        μ_0,
    )


    flux_up[:, :] .= fluxes.flux_up
    flux_dn[:, :] .= fluxes.flux_dn
    #
    # Zero out fluxes for which the original solar zenith angle is > 90 degrees.
    #
    for icol = 1:block_size
        if !usecol[icol]
            flux_up[icol, :] .= FT(0)
            flux_dn[icol, :] .= FT(0)
        end
    end

    if export_plots
        case = "clearsky_sw_" * string(optical_props_constructor)
        heating_rate, z =
            compute_heating_rate(fluxes.flux_up, fluxes.flux_dn, as)
        plot(
            heating_rate,
            z,
            title = "Clear sky shortwave heating rates",
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

    # comparing with reference data
    rsu_ref = ds[:flx_up]["rsu"][:]
    rsd_ref = ds[:flx_dn]["rsd"][:]

    rel_err_up, rel_err_dn = rel_err(flux_up, flux_dn, rsu_ref, rsd_ref)

    println("******************************************************")
    println("rfmip_clear_sky_sw -> $optical_props_constructor \n")
    println("Max rel_err_flux_up = $rel_err_up; Max rel_err_flux_dn = $rel_err_dn")
    println("******************************************************")

    if optical_props_constructor == TwoStream
        @test rel_err_up < FT(1e-7)
        @test rel_err_dn < FT(1e-7)
    else
        @test rel_err_up < FT(1.1)
        @test rel_err_dn < FT(1.1)
    end
    @show to
    return nothing
end
