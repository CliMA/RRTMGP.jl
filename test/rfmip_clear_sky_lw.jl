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
using RRTMGP.LoadCoefficients
using RRTMGP.RFMIPIO
using RRTMGP.SourceFunctions
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
function rfmip_clear_sky_lw(ds, optical_props_constructor; compile_first=false)

  FT = Float64
  I  = Int
  deg_to_rad = acos(-FT(1))/FT(180)
  n_quad_angles = I(1)

  ncol, nlay, nexp = read_size(ds[:rfmip])

  forcing_index = 1
  block_size = 8

  # How big is the problem? Does it fit into blocks of the size we've specified?
  @assert mod(ncol*nexp, block_size) == 0

  nblocks = Int((ncol*nexp)/block_size)
  # println("Doing $(nblocks) blocks of size $(block_size)")

  @assert 1 <= forcing_index <= 3

  # Identify the set of gases used in the calculation based on the forcing index
  #   A gas might have a different name in the k-distribution than in the files
  #   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')

  kdist_gas_names, rfmip_gas_games = determine_gas_names(ds[:k_dist], forcing_index)
  # print("Calculation uses RFMIP gases: $(rfmip_gas_games)")

  # Read atmospheric state from netCDF files
  p_lay, p_lev, t_lay, t_lev = read_and_block_pt(ds[:rfmip], block_size)

  # Read the gas concentrations and surface properties from netCDF files
  gas_conc_array = read_and_block_gases_ty(ds[:rfmip], block_size, kdist_gas_names, rfmip_gas_games)
  sfc_emis, sfc_t = read_and_block_lw_bc(ds[:rfmip], block_size)

  # Initialize atmospheric state, deallocate locals
  atmos_state = AtmosphericState(gas_conc_array,p_lay,p_lev,t_lay,t_lev)
  top_at_1 = atmos_state.top_at_1
  p_lay, t_lay, p_lev, gas_conc_array, t_lev = ntuple(i->nothing,5)

  # Read k-distribution information from netCDF files.
  k_dist = load_and_init(ds[:k_dist], FT, gas_conc_array[1].gas_name)

  @assert source_is_internal(k_dist)

  nbnd = get_nband(k_dist.optical_props)
  ngpt = get_ngpt(k_dist.optical_props)

  # RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
  #   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
  #   This introduces an error but shows input sanitizing.
  if top_at_1
    atmos_state.p_lev[:,1,:] .= get_press_min(k_dist) + eps(FT)
  else
    atmos_state.p_lev[:,nlay+1,:] .= get_press_min(k_dist) + eps(FT)
  end

  # RTE will fail if passed solar zenith angles greater than 90 degree. We replace any with
  #   nighttime columns with a default solar zenith angle. We'll mask these out later, of
  #   course, but this gives us more work and so a better measure of timing.

  # Allocate space for output fluxes (accessed via pointers in FluxesBroadBand),
  #   gas optical properties, and source functions. The %alloc() routines carry along
  #   the spectral discretization from the k-distribution.
  flux_up = Array{FT}(undef, block_size, nlay+1, nblocks)
  flux_dn = Array{FT}(undef, block_size, nlay+1, nblocks)

  sfc_emis_spec = Array{FT}(undef, nbnd,block_size)

  ps = ProblemSize(block_size, nlay, ngpt)
  optical_props = optical_props_constructor(k_dist.optical_props, ps)
  source = SourceFuncLW(block_size,nlay,k_dist.optical_props)

  fluxes = FluxesBroadBand(FT, (size(flux_up,1),size(flux_up,2)))

  # Loop over blocks
  for b = 1:(compile_first ? 1 : nblocks)
    fluxes.flux_up .= FT(0)
    fluxes.flux_dn .= FT(0)

    update_view!(atmos_state, b)

    for icol = 1:block_size
      for ibnd = 1:nbnd
        sfc_emis_spec[ibnd,icol] = sfc_emis[icol,b]
      end
    end

    gas_optics_int!(k_dist,
                    atmos_state,
                    sfc_t[:,  b],
                    optical_props,
                    source;
                    b=b)

    rte_lw!(optical_props,top_at_1,source,sfc_emis_spec,fluxes,nothing,n_quad_angles)

    flux_up[:,:,b] .= fluxes.flux_up
    flux_dn[:,:,b] .= fluxes.flux_dn

  end

  # reshaping the flux_up and flux_dn arrays for comparison with Fortran code.
  flux_up = reshape_for_comparison(flux_up, nlay, ncol, nexp)
  flux_dn = reshape_for_comparison(flux_dn, nlay, ncol, nexp)

  # comparing with reference data
  rlu_ref = ds[:flx_up]["rlu"][:]
  rld_ref = ds[:flx_dn]["rld"][:]

  diff_up = maximum( abs.( flux_up .- rlu_ref ) )
  diff_dn = maximum( abs.( flux_dn .- rld_ref ) )

  diff_up_ulps = maximum( abs.( flux_up .- rlu_ref ) ./ eps.(rlu_ref) )
  diff_dn_ulps = maximum( abs.( flux_dn .- rld_ref ) ./ eps.(rld_ref) )

  # @show sqrt(1/eps(FT))
  # @show diff_up, diff_up_ulps, maximum(abs.(rlu_ref))
  # @show diff_dn, diff_dn_ulps, maximum(abs.(rld_ref))

  if !compile_first
    if optical_props_constructor isa TwoStream
      @test diff_up_ulps < sqrt(1/(1e6eps(FT)))
      @test diff_dn_ulps < sqrt(1/(1e6eps(FT)))
    else
      @test diff_up_ulps < sqrt(1/(1e5eps(FT)))
      @test diff_dn_ulps < sqrt(1/(1e3eps(FT)))
    end
  end
  return nothing
end
