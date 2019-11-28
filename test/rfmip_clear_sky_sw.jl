using Test
using RRTMGP
using NCDatasets
using ProgressMeter
using TimerOutputs
const to = TimerOutput()
using RRTMGP.OpticalProps
using RRTMGP.FortranIntrinsics
using RRTMGP.ArrayUtilities
using RRTMGP.GasOptics
using RRTMGP.GasConcentrations
using RRTMGP.RTESolver
using RRTMGP.Fluxes


include(joinpath("..","ReadInputData","ReadInputs.jl"))
include(joinpath("..","ReadInputData","LoadCoefficients.jl"))

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
function rfmip_clear_sky_sw(ds, optical_props_constructor; compile_first=false)

  FT = Float64
  I  = Int
  deg_to_rad = acos(-FT(1))/FT(180)

  ncol, nlay, nexp = read_size(ds[:rfmip])

  forcing_index = 1
  block_size = 8

  # How big is the problem? Does it fit into blocks of the size we've specified?
  @assert mod(ncol*nexp, block_size) == 0 # number of columns must fit evenly into blocks
  nblocks = Int((ncol*nexp)/block_size)

  @assert !(forcing_index < 1 || forcing_index > 3)

  #
  # Identify the set of gases used in the calculation based on the forcing index
  #   A gas might have a different name in the k-distribution than in the files
  #   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
  #
  kdist_gas_names, rfmip_gas_games = @timeit to "determine_gas_names" determine_gas_names(ds[:k_dist], forcing_index)
  # print("Calculation uses RFMIP gases: ")
  # @show rfmip_gas_games

  # --------------------------------------------------
  #
  # Prepare data for use in rte+rrtmgp
  #
  #
  # Allocation on assignment within reading routines
  #
  p_lay, p_lev, t_lay, t_lev = @timeit to "read_and_block_pt" read_and_block_pt(ds[:rfmip], block_size)
  #
  # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
  #

  top_at_1 = p_lay[1, 1, 1] < p_lay[1, nlay, 1]

  #
  # Read the gas concentrations and surface properties
  #
  gas_conc_array = @timeit to "read_and_block_gases_ty" read_and_block_gases_ty(ds[:rfmip], block_size, kdist_gas_names, rfmip_gas_games)
  surface_albedo, total_solar_irradiance, solar_zenith_angle = @timeit to "read_and_block_sw_bc" read_and_block_sw_bc(ds[:rfmip], block_size)

  #
  # Read k-distribution information. load_and_init() reads data from netCDF and calls
  #   k_dist%init(); users might want to use their own reading methods
  #
  k_dist = @timeit to "load_and_init" load_and_init(ds[:k_dist], FT, gas_conc_array[1].gas_names)
  @assert source_is_external(k_dist)

  nbnd = get_nband(k_dist.optical_props)
  ngpt = get_ngpt(k_dist.optical_props)

  println("--------- Problem size:")
  @show ncol,nlay,nbnd,ngpt
  @show nblocks,nexp,block_size

  ps = ProblemSize(block_size, nlay, ngpt)

  toa_flux = zeros(FT, block_size, get_ngpt(k_dist.optical_props))
  def_tsi = zeros(FT, block_size)
  usecol = Array{Bool}(undef, block_size, nblocks)
  #
  # RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
  #   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
  #   This introduces an error but shows input sanitizing.
  #
  if top_at_1
    p_lev[:,1,:] .= get_press_min(k_dist) + eps(FT)
  else
    p_lev[:,nlay+1,:] .= get_press_min(k_dist) + eps(FT)
  end

  #
  # RTE will fail if passed solar zenith angles greater than 90 degree. We replace any with
  #   nighttime columns with a default solar zenith angle. We'll mask these out later, of
  #   course, but this gives us more work and so a better measure of timing.
  #
  for b = 1:nblocks
    usecol[1:block_size,b] .= solar_zenith_angle[1:block_size,b] .< FT(90) - FT(2) * spacing(FT(90))
  end

  #
  # Allocate space for output fluxes (accessed via pointers in FluxesBroadBand),
  #   gas optical properties, and source functions. The %alloc() routines carry along
  #   the spectral discretization from the k-distribution.
  #
  flux_up = zeros(FT, block_size, nlay+1, nblocks)
  flux_dn = zeros(FT, block_size, nlay+1, nblocks)

  μ_0 = zeros(FT, block_size)
  sfc_alb_spec = zeros(FT, nbnd,block_size)
  optical_props = optical_props_constructor(k_dist.optical_props, ps)

  #
  # Loop over blocks
  #
  fluxes = FluxesBroadBand(FT, (size(flux_up,1),size(flux_up,2)), true)

  b_tot = (compile_first ? 1 : nblocks)
  @showprogress 1 "Computing..." for b = 1:b_tot
    fluxes.flux_up .= FT(0)
    fluxes.flux_dn .= FT(0)
    #
    # Compute the optical properties of the atmosphere and the Planck source functions
    #    from pressures, temperatures, and gas concentrations...
    #

    @timeit to "gas_optics_ext!" gas_optics_ext!(k_dist,
                p_lay[:,:,b],
                p_lev[:,:,b],
                t_lay[:,:,b],
                gas_conc_array[b],
                optical_props,
                toa_flux,
                nothing,
                b==b_tot)
    # Boundary conditions
    #   (This is partly to show how to keep work on GPUs using OpenACC in a host application)
    # What's the total solar irradiance assumed by RRTMGP?
    #
    def_tsi[1:block_size] = sum(toa_flux, dims=2)
    #
    # Normalize incoming solar flux to match RFMIP specification
    #
    for igpt = 1:ngpt
      for icol = 1:block_size
        toa_flux[icol,igpt] = toa_flux[icol,igpt] * total_solar_irradiance[icol,b]/def_tsi[icol]
      end
    end
    #
    # Expand the spectrally-constant surface albedo to a per-band albedo for each column
    #
    for icol = 1:block_size
      for ibnd = 1:nbnd
        sfc_alb_spec[ibnd,icol] = surface_albedo[icol,b]
      end
    end
    #
    # Cosine of the solar zenith angle
    #
    for icol = 1:block_size
      μ_0[icol] = fmerge(cos(solar_zenith_angle[icol,b]*deg_to_rad), FT(1), usecol[icol,b])
    end

    #
    # ... and compute the spectrally-resolved fluxes, providing reduced values
    #    via FluxesBroadBand
    #

    @timeit to "rte_sw!" rte_sw!(optical_props,
            top_at_1,
            μ_0,
            toa_flux,
            sfc_alb_spec,
            sfc_alb_spec,
            fluxes)


    flux_up[:,:,b] .= fluxes.flux_up
    flux_dn[:,:,b] .= fluxes.flux_dn
    #
    # Zero out fluxes for which the original solar zenith angle is > 90 degrees.
    #
    for icol = 1:block_size
      if !usecol[icol,b]
        flux_up[icol,:,b] .= FT(0)
        flux_dn[icol,:,b] .= FT(0)
      end
    end

  end

  # reshaping the flux_up and flux_dn arrays for comparison with Fortran code.
  flux_up = reshape_for_comparison(flux_up, nlay, ncol, nexp)
  flux_dn = reshape_for_comparison(flux_dn, nlay, ncol, nexp)

  # comparing with reference data
  rsu_ref = ds[:flx_up]["rsu"][:]
  rsd_ref = ds[:flx_dn]["rsd"][:]

  diff_up = maximum( abs.( flux_up .- rsu_ref ) )
  diff_dn = maximum( abs.( flux_dn .- rsd_ref ) )

  diff_up_ulps = maximum( abs.( flux_up .- rsu_ref ) ./ eps.(rsu_ref) )
  diff_dn_ulps = maximum( abs.( flux_dn .- rsd_ref ) ./ eps.(rsd_ref) )

  # @show sqrt(1/eps(FT))
  # @show diff_up, diff_up_ulps, maximum(abs.(rsu_ref))
  # @show diff_dn, diff_dn_ulps, maximum(abs.(rsd_ref))

  if !compile_first
    if optical_props_constructor isa TwoStream
      @test diff_up_ulps < sqrt(1/(1e6eps(FT)))
      @test diff_dn_ulps < sqrt(1/(1e6eps(FT)))
    else
      @test diff_up_ulps < sqrt(1/(eps(FT))) # 1.6776966e7
      @test diff_dn_ulps < sqrt(1/(eps(FT))) # 1.6777158e7
    end
  end
  @show to
  return nothing
end
