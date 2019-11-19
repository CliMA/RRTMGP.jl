using Test
using RRTMGP
using NCDatasets
using RRTMGP.mo_optical_props
using RRTMGP.mo_rte_solver_kernels
using RRTMGP.fortran_intrinsics
using RRTMGP.mo_util_array
using RRTMGP.mo_gas_optics_rrtmgp
using RRTMGP.mo_gas_concentrations
using RRTMGP.mo_rte_lw
using RRTMGP.mo_fluxes
using RRTMGP.mo_load_coefficients
using RRTMGP.mo_rfmip_io
using RRTMGP.mo_source_functions
"""
Example program to demonstrate the calculation of shortwave radiative fluxes in clear, aerosol-free skies.
  The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
  The large problem (1800 profiles) is divided into blocks

Program is invoked as rrtmgp_rfmip_sw [block_size input_file  coefficient_file upflux_file downflux_file]
  All arguments are optional but need to be specified in order.

Optical properties of the atmosphere as array of values
   In the longwave we include only absorption optical depth (_1scl)
   Shortwave calculations use optical depth, single-scattering albedo, asymmetry parameter (_2str)

RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  Here we're just reporting broadband fluxes

RRTMGP's gas optics class needs to be initialized with data read from a netCDF files
"""
function rfmip_clear_sky_lw(ds, optical_props_constructor; compile_first=false)
  # character(len=132) :: flxdn_file, flxup_file
  # integer            :: nargs, ncol, nlay, nbnd, ngpt, nexp, nblocks, block_size, forcing_index
  # logical            :: top_at_1
  # integer            :: b, icol, ibnd, igpt
  # character(len=4)   :: block_size_char, forcing_index_char = '1'

  # character(len=32 ), dimension(:),   allocatable :: kdist_gas_names, rfmip_gas_games
  # real(FT), dimension(:,:,:),         allocatable :: p_lay, p_lev, t_lay, t_lev # block_size, nlay, nblocks
  # real(FT), dimension(:,:,:), target, allocatable :: flux_up, flux_dn
  # real(FT), dimension(:,:  ),         allocatable :: surface_albedo, total_solar_irradiance, solar_zenith_angle
  #                                                    # block_size, nblocks
  # real(FT), dimension(:,:  ),         allocatable :: sfc_alb_spec # nbnd, block_size; spectrally-resolved surface albedo
  #
  # Classes used by rte+rrtmgp
  #
  # type(ty_gas_optics_rrtmgp)                     :: k_dist
  # type(ty_optical_props_2str)                    :: optical_props
  # type(ty_fluxes_broadband)                      :: fluxes

  # real(FT), dimension(:,:), allocatable          :: toa_flux # block_size, ngpt
  # real(FT), dimension(:  ), allocatable          :: def_tsi, mu0    # block_size
  # logical , dimension(:,:), allocatable          :: usecol # block_size, nblocks
  #
  # ty_gas_concentration holds multiple columns; we make an array of these objects to
  #   leverage what we know about the input file
  #
  # type(ty_gas_concs), dimension(:), allocatable  :: gas_conc_array
  FT = Float64
  I  = Integer
  deg_to_rad = acos(-FT(1))/FT(180)
  n_quad_angles = Integer.(1)

  #
  # Code starts
  #   all arguments are optional
  #


  ncol, nlay, nexp = read_size(ds[:rfmip])

  forcing_index = 1
  block_size = 8

  #
  # How big is the problem? Does it fit into blocks of the size we've specified?
  #
  @assert mod(ncol*nexp, block_size) == 0

  nblocks = Int((ncol*nexp)/block_size)
  # println("Doing $(nblocks) blocks of size $(block_size)")

  @assert 1 <= forcing_index <= 3

  #
  # Identify the set of gases used in the calculation based on the forcing index
  #   A gas might have a different name in the k-distribution than in the files
  #   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
  #

  kdist_gas_names, rfmip_gas_games = determine_gas_names(ds[:k_dist], forcing_index)
  # print("Calculation uses RFMIP gases: ")
  # @show rfmip_gas_games

  # --------------------------------------------------
  #
  # Prepare data for use in rte+rrtmgp
  #
  #
  # Allocation on assignment within reading routines
  #
  p_lay, p_lev, t_lay, t_lev = read_and_block_pt(ds[:rfmip], block_size)
  #
  # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
  #

  top_at_1 = p_lay[1, 1, 1] < p_lay[1, nlay, 1]

  #
  # Read the gas concentrations and surface properties
  #
  gas_conc_array = read_and_block_gases_ty(ds[:rfmip], block_size, kdist_gas_names, rfmip_gas_games)
  sfc_emis, sfc_t = read_and_block_lw_bc(ds[:rfmip], block_size)

  #
  # Read k-distribution information. load_and_init() reads data from netCDF and calls
  #   k_dist%init(); users might want to use their own reading methods
  #
  k_dist = load_and_init(ds[:k_dist], gas_conc_array[1])

  @assert source_is_internal(k_dist)

  nbnd = get_nband(k_dist.optical_props)
  ngpt = get_ngpt(k_dist.optical_props)
  # @show nbnd
  # @show ngpt

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

  #
  # Allocate space for output fluxes (accessed via pointers in ty_fluxes_broadband),
  #   gas optical properties, and source functions. The %alloc() routines carry along
  #   the spectral discretization from the k-distribution.
  #
  flux_up = Array{FT}(undef, block_size, nlay+1, nblocks)
  flux_dn = Array{FT}(undef, block_size, nlay+1, nblocks)

  sfc_emis_spec = Array{FT}(undef, nbnd,block_size)

  ps = ProblemSize(block_size, nlay, ngpt)
  optical_props = optical_props_constructor(k_dist.optical_props, ps)
  source = ty_source_func_lw(block_size,nlay,k_dist.optical_props)

  #
  # Loop over blocks
  #
  fluxes = ty_fluxes_broadband(FT, (size(flux_up,1),size(flux_up,2)))

  for b = 1:(compile_first ? 1 : nblocks)
    fluxes.flux_up .= FT(0)
    fluxes.flux_dn .= FT(0)

    for icol = 1:block_size
      for ibnd = 1:nbnd
        sfc_emis_spec[ibnd,icol] = sfc_emis[icol,b]
      end
    end

    gas_optics_int!(k_dist,
                p_lay[:,:,b],
                p_lev[:,:,b],
                t_lay[:,:,b],
                sfc_t[:,  b],
                gas_conc_array[b],
                optical_props,
                source;
                tlev = t_lev[:,:,b])

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
    if optical_props_constructor isa ty_optical_props_2str
      @test diff_up_ulps < sqrt(1/(1e6eps(FT)))
      @test diff_dn_ulps < sqrt(1/(1e6eps(FT)))
    else
      @test diff_up_ulps < sqrt(1/(1e5eps(FT)))
      @test diff_dn_ulps < sqrt(1/(1e3eps(FT)))
    end
  end
  return nothing
end
