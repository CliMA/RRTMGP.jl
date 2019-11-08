using Test
using JRRTMGP
using NCDatasets
using JRRTMGP.mo_optical_props
using JRRTMGP.mo_rte_solver_kernels
using JRRTMGP.fortran_intrinsics
using JRRTMGP.mo_util_array
using JRRTMGP.mo_gas_optics_rrtmgp
using JRRTMGP.mo_gas_concentrations
using JRRTMGP.mo_rte_sw
using JRRTMGP.mo_fluxes
using JRRTMGP.mo_load_coefficients
using JRRTMGP.mo_rfmip_io

# This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
#
# Contacts: Robert Pincus and Eli Mlawer
# email:  rrtmgp@aer.com
#
# Copyright 2015-2018,  Atmospheric and Environmental Research and
# Regents of the University of Colorado.  All right reserved.
#
# Use and duplication is permitted under the terms of the
#    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
# -------------------------------------------------------------------------------------------------
#
# Example program to demonstrate the calculation of shortwave radiative fluxes in clear, aerosol-free skies.
#   The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
#   The large problem (1800 profiles) is divided into blocks
#
# Program is invoked as rrtmgp_rfmip_sw [block_size input_file  coefficient_file upflux_file downflux_file]
#   All arguments are optional but need to be specified in order.
#

function run_driver(datafolder, optical_props_constructor, nblocks_iterations=nothing)
  #
  # RTE shortwave driver
  #
  # RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  #   Here we're just reporting broadband fluxes
  #
  # modules for reading and writing files
  #
  # RRTMGP's gas optics class needs to be initialized with data read from a netCDF files
  #

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
  deg_to_rad = acos(-FT(1))/FT(180)

  # -------------------------------------------------------------------------------------------------
  #
  # Code starts
  #   all arguments are optional
  #

  rfmip_file = joinpath(datafolder, "examples","rfmip-clear-sky", "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc")
  kdist_file = joinpath(datafolder, "rrtmgp", "data", "rrtmgp-data-sw-g224-2018-12-04.nc")
  clear_sky_dir = joinpath(datafolder, "examples","rfmip-clear-sky")

  # sw_flx_up_for_res_file = joinpath("..","..","rte-rrtmgp", "examples","rfmip-clear-sky", "rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc") # Results from Fortran code
  # sw_flx_dn_for_res_file = joinpath("..","..","rte-rrtmgp", "examples","rfmip-clear-sky", "rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc") # Results from Fortran code
  sw_flx_up_for_res_file = joinpath(clear_sky_dir, "rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc") # Results from Fortran code
  sw_flx_dn_for_res_file = joinpath(clear_sky_dir, "rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc") # Results from Fortran code

  ds_sw_flx_up = Dataset(sw_flx_up_for_res_file, "r") # reading the NetCDF file in read only mode
  ds_sw_flx_dn = Dataset(sw_flx_dn_for_res_file, "r") # reading the NetCDF file in read only mode
  ds = Dataset(rfmip_file, "r") # reading the NetCDF file in read only mode
  ds_k_dist = Dataset(kdist_file, "r") # reading the NetCDF file in read only mode

  ncol, nlay, nexp = read_size(ds)

  forcing_index = 1
  block_size = 8

  #
  # How big is the problem? Does it fit into blocks of the size we've specified?
  #
  if mod(ncol*nexp, block_size) ≠ 0
    error("rrtmgp_rfmip_sw: number of columns does not fit evenly into blocks.")
  end
  nblocks = Int((ncol*nexp)/block_size)
  println("Doing $(nblocks) blocks of size $(block_size)")

  # TODO: Fix readme
  # read(forcing_index_char, "(i4)") forcing_index
  @assert !(forcing_index < 1 || forcing_index > 3)

  #
  # Identify the set of gases used in the calculation based on the forcing index
  #   A gas might have a different name in the k-distribution than in the files
  #   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
  #
  # names_in_kdist, names_in_file = determine_gas_names(ds_k_dist, forcing_index)
  kdist_gas_names, rfmip_gas_games = determine_gas_names(ds_k_dist, forcing_index)
  print("Calculation uses RFMIP gases: ")
  @show rfmip_gas_games

  # --------------------------------------------------
  #
  # Prepare data for use in rte+rrtmgp
  #
  #
  # Allocation on assignment within reading routines
  #
  p_lay, p_lev, t_lay, t_lev = read_and_block_pt(ds, block_size)
  #
  # Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
  #

  top_at_1 = p_lay[1, 1, 1] < p_lay[1, nlay, 1]

  #
  # Read the gas concentrations and surface properties
  #
  gas_conc_array = read_and_block_gases_ty(ds, block_size, kdist_gas_names, rfmip_gas_games)
  surface_albedo, total_solar_irradiance, solar_zenith_angle = read_and_block_sw_bc(ds, block_size)

  #
  # Read k-distribution information. load_and_init() reads data from netCDF and calls
  #   k_dist%init(); users might want to use their own reading methods
  #
  k_dist = load_and_init(ds_k_dist, gas_conc_array[1])
  !source_is_external(k_dist) && error("rrtmgp_rfmip_sw: k-distribution file is not SW")

  nbnd = get_nband(k_dist.optical_props)
  ngpt = get_ngpt(k_dist.optical_props)

  toa_flux = Array{FT}(undef, block_size, get_ngpt(k_dist.optical_props))
  def_tsi = Array{FT}(undef, block_size)
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
  # Allocate space for output fluxes (accessed via pointers in ty_fluxes_broadband),
  #   gas optical properties, and source functions. The %alloc() routines carry along
  #   the spectral discretization from the k-distribution.
  #
  flux_up = zeros(FT, block_size, nlay+1, nblocks)
  flux_dn = zeros(FT, block_size, nlay+1, nblocks)

  mu0 = Array{FT}(undef, block_size)
  sfc_alb_spec = Array{FT}(undef, nbnd,block_size)
  optical_props = optical_props_constructor(FT, Int)
  copy_and_alloc!(optical_props, block_size, nlay, k_dist.optical_props)

  #
  # Loop over blocks
  #
  fluxes = ty_fluxes_broadband(FT)

  nblocks_iterations==nothing && (nblocks_iterations = nblocks)
  for b = 1:nblocks_iterations
    @show b/nblocks
    fup = fluxes.flux_up = @view(flux_up[:,:,b])
    fdn = fluxes.flux_dn = @view(flux_dn[:,:,b])
    #
    # Compute the optical properties of the atmosphere and the Planck source functions
    #    from pressures, temperatures, and gas concentrations...
    #

    gas_optics_ext!(k_dist,
                p_lay[:,:,b],
                p_lev[:,:,b],
                t_lay[:,:,b],
                gas_conc_array[b],
                optical_props,
                toa_flux)
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
      mu0[icol] = fmerge(cos(solar_zenith_angle[icol,b]*deg_to_rad), FT(1), usecol[icol,b])
    end

    #
    # ... and compute the spectrally-resolved fluxes, providing reduced values
    #    via ty_fluxes_broadband
    #

    rte_sw!(optical_props,
            top_at_1,
            mu0,
            toa_flux,
            sfc_alb_spec,
            sfc_alb_spec,
            fluxes)


    #
    # Zero out fluxes for which the original solar zenith angle is > 90 degrees.
    #
    for icol = 1:block_size
      if !usecol[icol,b]
        flux_up[icol,:,b] .= FT(0)
        flux_dn[icol,:,b] .= FT(0)
      end
    end

    @assert fup === fluxes.flux_up
    @assert fdn === fluxes.flux_dn
  end

  # reshaping the flux_up and flux_dn arrays for comparison with Fortran code.
  flux_up = reshape_for_comparison(flux_up, nlay, ncol, nexp)
  flux_dn = reshape_for_comparison(flux_dn, nlay, ncol, nexp)

  rsu_for = ds_sw_flx_up["rsu"][:]
  rsd_for = ds_sw_flx_dn["rsd"][:]

  diff_up = maximum( abs.( flux_up .- rsu_for ) )
  diff_dn = maximum( abs.( flux_dn .- rsd_for ) )

  diff_up_ulps = maximum( abs.( flux_up .- rsu_for ) ./ eps.(rsu_for) )
  diff_dn_ulps = maximum( abs.( flux_dn .- rsd_for ) ./ eps.(rsd_for) )

  # @show sqrt(1/eps(FT))
  # @show diff_up, diff_up_ulps, maximum(abs.(rsu_for))
  # @show diff_dn, diff_dn_ulps, maximum(abs.(rsd_for))

  if optical_props_constructor isa ty_optical_props_1scl
    @test diff_up_ulps < sqrt(1/(eps(FT)))
    @test diff_dn_ulps < sqrt(1/(eps(FT)))
  else
    @test diff_up_ulps < sqrt(1/(1e6eps(FT)))
    @test diff_dn_ulps < sqrt(1/(1e6eps(FT)))
  end

  # flxdn_file = "rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f" * string(forcing_index) * "_gn.nc"
  # flxup_file = "rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f" * string(forcing_index) * "_gn.nc"
  # unblock_and_write!(trim(flxup_file), "rsu", flux_up)
  # unblock_and_write!(trim(flxdn_file), "rsd", flux_dn)

  close(ds_sw_flx_up)
  close(ds_sw_flx_dn)
  close(ds)
  close(ds_k_dist)
end

@testset "Shortwave driver ty_optical_props_1scl" begin
  datafolder = JRRTMGP.data_folder_rrtmgp()

  run_driver(datafolder, ty_optical_props_1scl)
end

@testset "Shortwave driver ty_optical_props_2str" begin
  datafolder = JRRTMGP.data_folder_rrtmgp()

  run_driver(datafolder, ty_optical_props_2str)
end
