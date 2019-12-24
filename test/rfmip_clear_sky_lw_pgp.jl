using Test
using RRTMGP
using NCDatasets
using ProgressMeter
using TimerOutputs
const to = TimerOutput()
using RRTMGP.OpticalProps
using RRTMGP.FortranIntrinsics
using RRTMGP.Utilities
using RRTMGP.GasOptics
using RRTMGP.GasConcentrations
using RRTMGP.RadiativeBoundaryConditions
using RRTMGP.RTESolver
using RRTMGP.Fluxes
using RRTMGP.AtmosphericStates
using RRTMGP.SourceFunctions
using RRTMGP.AngularDiscretizations
@static if haspkg("Plots")
  using Plots
  const export_plots = true
else
  const export_plots = false
end

include(joinpath("PostProcessing.jl"))
include(joinpath("..","ReadInputData","ReadInputs.jl"))
include(joinpath("..","ReadInputData","LoadCoefficients.jl"))

convert_optical_props_pgp(op::AbstractOpticalPropsArry) =
  op isa TwoStream ? convert(Array{TwoStreamPGP},op) :
                     convert(Array{OneScalarPGP},op)
convert_optical_props(op) =
  eltype(op) <: TwoStreamPGP ? convert(TwoStream,op) :
                                convert(OneScalar,op)

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
function rfmip_clear_sky_lw_pgp(ds, optical_props_constructor; compile_first=false)

  FT = Float64
  I  = Int
  deg_to_rad = acos(-FT(1))/FT(180)
  angle_disc = GaussQuadrature(FT, 1)

  ncol, nlay, nexp = read_size(ds[:rfmip])

  forcing_index = 1
  block_size = 8

  # How big is the problem? Does it fit into blocks of the size we've specified?
  @assert mod(ncol*nexp, block_size) == 0
  nblocks = Int((ncol*nexp)/block_size)

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
  p_lay_all, p_lev_all, t_lay_all, t_lev_all = read_and_block_pt(ds[:rfmip], block_size)
  #
  # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
  #

  top_at_1 = p_lay_all[1, 1, 1] < p_lay_all[1, nlay, 1]

  #
  # Read the gas concentrations and surface properties
  #
  gas_conc_array = read_and_block_gases_ty(ds[:rfmip], block_size, kdist_gas_names)
  sfc_emis_all, t_sfc_all = read_and_block_lw_bc(ds[:rfmip], block_size)

  #
  # Read k-distribution information. load_and_init() reads data from netCDF and calls
  #   k_dist%init(); users might want to use their own reading methods
  #
  k_dist = load_and_init(ds[:k_dist], FT, gas_conc_array[1].gas_names)
  @assert source_is_internal(k_dist)

  nbnd = get_nband(k_dist.optical_props)
  ngpt = get_ngpt(k_dist.optical_props)

  #
  # RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
  #   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
  #   This introduces an error but shows input sanitizing.
  #
  if top_at_1
    p_lev_all[:,1,:] .= get_press_min(k_dist.ref) + eps(FT)
  else
    p_lev_all[:,nlay+1,:] .= get_press_min(k_dist.ref) + eps(FT)
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
  flux_up = Array{FT}(undef, block_size, nlay+1, nblocks)
  flux_dn = Array{FT}(undef, block_size, nlay+1, nblocks)

  sfc_emis_spec = Array{FT}(undef, nbnd,block_size)

  optical_props = optical_props_constructor(k_dist.optical_props, block_size, nlay, ngpt)
  optical_props = convert_optical_props_pgp(optical_props)
  optical_props = convert_optical_props(optical_props)

  source = SourceFuncLongWave(block_size,nlay,k_dist.optical_props)

  #
  # Loop over blocks
  #
  fluxes = FluxesBroadBand(FT, (size(flux_up,1),size(flux_up,2)))

  local as

  for b = 1:(compile_first ? 1 : nblocks)
    for icol = 1:block_size
      for ibnd = 1:nbnd
        sfc_emis_spec[ibnd,icol] = sfc_emis_all[icol,b]
      end
    end
    gas_conc = gas_conc_array[b]
    p_lay = p_lay_all[:,:,b]
    p_lev = p_lev_all[:,:,b]
    t_lay = t_lay_all[:,:,b]
    t_lev = t_lev_all[:,:,b]
    t_sfc = t_sfc_all[:,  b]
    as = AtmosphericState(gas_conc, p_lay, p_lev, t_lay, t_lev, k_dist.ref, nothing, t_sfc)
    as = convert(Array{AtmosphericStatePGP}, as)
    as = convert(AtmosphericState, as)

    fluxes.flux_up .= FT(0)
    fluxes.flux_dn .= FT(0)

    source = convert(Array{SourceFuncLongWavePGP}, source)
    source = convert(SourceFuncLongWave, source)

    gas_optics!(k_dist, as, optical_props, source)

    bcs = LongwaveBCs(sfc_emis_spec)

    rte_lw!(optical_props,
            as.top_at_1,
            source,
            bcs,
            fluxes,
            angle_disc)

    flux_up[:,:,b] .= fluxes.flux_up
    flux_dn[:,:,b] .= fluxes.flux_dn

  end

  if export_plots
    case = "clearsky_lw_"*string(optical_props_constructor)
    heating_rate, z = compute_heating_rate(fluxes.flux_up, fluxes.flux_dn, as)
    plot(heating_rate, z, title="Clear sky longwave heating rates",
                          xlabel="heating rate",
                          ylabel="pressure")
    out_dir = "output"
    mkpath(out_dir)
    savefig(joinpath(out_dir,case*".png"))
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
