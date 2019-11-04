using Test
using JRRTMGP
using NCDatasets
using JRRTMGP.mo_optical_props
using JRRTMGP.mo_simple_netcdf
using JRRTMGP.mo_rte_solver_kernels
using JRRTMGP.fortran_intrinsics
using JRRTMGP.mo_util_array
using JRRTMGP.mo_gas_optics_rrtmgp
using JRRTMGP.mo_gas_concentrations
using JRRTMGP.mo_rte_lw
using JRRTMGP.mo_fluxes
using JRRTMGP.mo_load_coefficients
using JRRTMGP.mo_rfmip_io
using JRRTMGP.mo_source_functions
using JRRTMGP.mo_cloud_optics
using JRRTMGP.mo_rte_sw
using JRRTMGP.mo_load_cloud_coefficients
using JRRTMGP.mo_garand_atmos_io

include("mo_test_files_io.jl")

function vmr_2d_to_1d!(gas_concs, gas_concs_garand, name)
  # use mo_gas_concentrations, only: ty_gas_concs

  # type(ty_gas_concs), intent(in)    :: gas_concs_garand
  # type(ty_gas_concs), intent(inout) :: gas_concs
  # character(len=*),   intent(in)    :: name
  # integer,            intent(in)    :: sz1, sz2

  # real(wp) :: tmp(sz1, sz2), tmp_col(sz2)

  tmp = get_vmr(gas_concs_garand, name)
  tmp_col = tmp[1, :]

  set_vmr!(gas_concs, name, tmp_col)
end

function run_driver(datafolder; λ_string="")
  @assert λ_string == "sw" || λ_string == "lw"

  # # ----------------------------------------------------------------------------------
  # # Variables
  # # ----------------------------------------------------------------------------------
  # # Arrays: dimensions (col, lay)
  # real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev
  # real(wp), dimension(:,:),   allocatable :: col_dry
  # real(wp), dimension(:,:),   allocatable :: temp_array

  # #
  # # Longwave only
  # #
  # real(wp), dimension(:,:),   allocatable :: t_lev
  # real(wp), dimension(:),     allocatable :: t_sfc
  # real(wp), dimension(:,:),   allocatable :: emis_sfc # First dimension is band
  # #
  # # Shortwave only
  # #
  # real(wp), dimension(:),     allocatable :: mu0
  # real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif # First dimension is band
  # #
  # # Source functions
  # #
  # #   Longwave
  # type(ty_source_func_lw), save               :: lw_sources
  # #   Shortwave
  # real(wp), dimension(:,:), allocatable, save :: toa_flux
  # #
  # # Clouds
  # #
  # real(wp), allocatable, dimension(:,:) :: lwp, iwp, rel, rei
  # logical,  allocatable, dimension(:,:) :: cloud_mask
  # #
  # # Output variables
  # #
  # real(wp), dimension(:,:), target,
  #                           allocatable :: flux_up, flux_dn, flux_dir
  # #
  # # Derived types from the RTE and RRTMGP libraries
  # #
  # type(ty_gas_optics_rrtmgp) :: k_dist
  # type(ty_cloud_optics)      :: cloud_optics_
  # type(ty_gas_concs)         :: gas_concs, gas_concs_garand, gas_concs_1col
  # class(ty_optical_props_arry),
  #                allocatable :: atmos, clouds
  # type(ty_fluxes_broadband)  :: fluxes

  # #
  # # Inputs to RRTMGP
  # #
  # logical :: top_at_1, is_sw, is_lw

  # integer  :: ncol, nlay, nbnd, ngpt
  # integer  :: icol, ilay, ibnd, iloop, igas
  # real(wp) :: rel_val, rei_val

  # character(len=8) :: char_input
  # integer  :: nUserArgs=0, nloops
  # logical :: use_luts = .true., write_fluxes = .true.
  # integer, parameter :: ngas = 8
  # character(len=3), dimension(ngas)
  #                    :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']

  # character(len=256) :: input_file, k_dist_file, cloud_optics_file
  # #
  # # Timing variables
  # #
  # integer(kind=8)              :: start, finish, start_all, finish_all, clock_rate
  # real(wp)                     :: avg
  # integer(kind=8), allocatable :: elapsed(:)

  #
  # Parse command line for any file names, block size
  #
  # rrtmgp_clouds rrtmgp-clouds.nc $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc  128 1
  # rrtmgp_clouds rrtmgp-clouds.nc $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc  128 1

  gas_names = ["h2o", "co2", "o3 ", "n2o", "co ", "ch4", "o2 ", "n2 "]
  ngas = length(gas_names)
  nloops = 1
  ncol = 128
  use_luts = true
  write_fluxes = true
  λ_string=="lw" && (n_g_points = "256")
  λ_string=="sw" && (n_g_points = "224")
  input_file = joinpath(datafolder, "examples", "all-sky", "rrtmgp-allsky.nc")
  reference_file = joinpath(datafolder, "examples", "all-sky", "ref", "rrtmgp-allsky.nc")
  k_dist_file = joinpath(datafolder, "rrtmgp", "data", "rrtmgp-data-$(λ_string)-g"*n_g_points*"-2018-12-04.nc")
  cloud_optics_file = joinpath(datafolder, "extensions", "cloud_optics", "rrtmgp-cloud-optics-coeffs-$(λ_string).nc")

  #
  # Read temperature, pressure, gas concentrations.
  #   Arrays are allocated as they are read
  #
  FT = Float64
  I = Int64
  ds_input = Dataset(input_file, "r")
  ds_ref = Dataset(reference_file, "r")
  ds_k_dist = Dataset(k_dist_file, "r")
  ds_cloud_optics = Dataset(cloud_optics_file, "r")
  p_lay, t_lay, p_lev, t_lev, gas_concs_garand, col_dry = read_atmos(ds_input, FT)
  deallocate!(col_dry)
  nlay = size(p_lay, 2)
  # For clouds we'll use the first column, repeated over and over
  gas_concs = ty_gas_concs(FT, ncol, nlay)
  gas_concs.gas_name = gas_names
  for igas = 1:ngas
    vmr_2d_to_1d!(gas_concs, gas_concs_garand, gas_names[igas])
  end
  #  If we trusted in Fortran allocate-on-assign we could skip the temp_array here
  temp_array = Array{FT}(undef, ncol, nlay)
  temp_array .= reshape(spread(p_lay[1,:], 1, ncol), ncol, nlay)
  p_lay = deepcopy(temp_array)
  temp_array = Array{FT}(undef, ncol, nlay)
  temp_array .= reshape(spread(t_lay[1,:], 1, ncol), ncol, nlay)
  t_lay = deepcopy(temp_array)
  temp_array = Array{FT}(undef, ncol, nlay+1)
  temp_array .= reshape(spread(p_lev[1,:], 1, ncol), ncol, nlay+1)
  p_lev = deepcopy(temp_array)
  temp_array = Array{FT}(undef, ncol, nlay+1)
  temp_array .= reshape(spread(t_lev[1,:], 1, ncol), ncol, nlay+1)
  t_lev = deepcopy(temp_array)

  # This puts pressure and temperature arrays on the GPU
  # load data into classes
  k_dist = load_and_init(ds_k_dist, gas_concs)
  is_sw = source_is_external(k_dist)
  is_lw = !is_sw
  #
  # Should also try with Pade calculations
  cloud_optics_ = ty_cloud_optics(FT,I)
  if use_luts
    load_cld_lutcoeff!(cloud_optics_, ds_cloud_optics)
  else
    load_cld_padecoeff!(cloud_optics_, ds_cloud_optics)
  end
  set_ice_roughness!(cloud_optics_, 2)

  #
  # Problem sizes
  #
  nbnd = get_nband(k_dist.optical_props)
  ngpt = get_ngpt(k_dist.optical_props)
  top_at_1 = p_lay[1, 1] < p_lay[1, nlay]


  # LW calculations neglect scattering; SW calculations use the 2-stream approximation
  #   Here we choose the right variant of optical_props.
  #
  if is_sw
    clouds = ty_optical_props_2str(FT, I)
    atmos = ty_optical_props_2str(FT, I)
  else
    atmos = ty_optical_props_1scl(FT, I)
    clouds = ty_optical_props_1scl(FT, I)
  end

  # Clouds optical props are defined by band
  init!(clouds, "Clouds", get_band_lims_wavenumber(k_dist.optical_props))

  #
  # Allocate arrays for the optical properties themselves.
  #
  copy_and_alloc!(atmos,  ncol, nlay, k_dist.optical_props)
  alloc!(clouds, ncol, nlay)

  #  Boundary conditions depending on whether the k-distribution being supplied
  #   is LW or SW
  if is_sw
    # toa_flux is thread private
    toa_flux = Array{FT}(undef, ncol, ngpt)
    #
    sfc_alb_dir = Array{FT}(undef, nbnd, ncol)
    sfc_alb_dif = Array{FT}(undef, nbnd, ncol)
    mu0 = Array{FT}(undef, ncol)
    # Ocean-ish values for no particular reason
    sfc_alb_dir .= FT(0.06)
    sfc_alb_dif .= FT(0.06)
    mu0 .= FT(.86)
  else
    # lw_sorces is thread private
    lw_sources = ty_source_func_lw(ncol, nlay, k_dist.optical_props)

    t_sfc = Array{FT}(undef, ncol)
    emis_sfc = Array{FT}(undef, nbnd, ncol)
    # Surface temperature
    t_sfc .= t_lev[1, fmerge(nlay+1, 1, top_at_1)]
    emis_sfc .= FT(0.98)
  end

  #
  # Fluxes
  #
  flux_up = Array{FT}(undef, ncol,nlay+1)
  flux_dn = Array{FT}(undef, ncol,nlay+1)

  if is_sw
    flux_dir = Array{FT}(undef, ncol,nlay+1)
  end
  #
  # Clouds
  #
  lwp = Array{FT}(undef, ncol,nlay)
  iwp = Array{FT}(undef, ncol,nlay)
  rel = Array{FT}(undef, ncol,nlay)
  rei = Array{FT}(undef, ncol,nlay)
  cloud_mask = Array{Bool}(undef, ncol,nlay)

  # Restrict clouds to troposphere (< 100 hPa = 100*100 Pa)
  #   and not very close to the ground
  rel_val = FT(0.5) * (get_min_radius_liq(cloud_optics_) + get_max_radius_liq(cloud_optics_))
  rei_val = FT(0.5) * (get_min_radius_ice(cloud_optics_) + get_max_radius_ice(cloud_optics_))
  for ilay=1:nlay
    for icol=1:ncol
      cloud_mask[icol,ilay] = p_lay[icol,ilay] < FT(100) * FT(100) && p_lay[icol,ilay] > FT(900)
      #
      # Ice and liquid will overlap in a few layers
      #
      lwp[icol,ilay] = fmerge(FT(10),  FT(0), cloud_mask[icol,ilay] && t_lay[icol,ilay] > FT(263))
      iwp[icol,ilay] = fmerge(FT(10),  FT(0), cloud_mask[icol,ilay] && t_lay[icol,ilay] < FT(273))
      rel[icol,ilay] = fmerge(rel_val, FT(0), lwp[icol,ilay] > FT(0))
      rei[icol,ilay] = fmerge(rei_val, FT(0), iwp[icol,ilay] > FT(0))
    end
  end
  fluxes = ty_fluxes_broadband(FT)
  #
  # Multiple iterations for big problem sizes, and to help identify data movement
  #   For CPUs we can introduce OpenMP threading over loop iterations
  #

  for iloop = 1:nloops
    cloud_optics!(cloud_optics_, lwp, iwp, rel, rei, clouds)
    #
    # Solvers
    #
    fluxes.flux_up = flux_up[:,:]
    fluxes.flux_dn = flux_dn[:,:]
    if is_lw
      gas_optics!(k_dist, p_lay, p_lev,
                  t_lay, t_sfc,
                  gas_concs,
                  atmos,
                  lw_sources;
                  tlev = t_lev)
      increment!(clouds, atmos)
      rte_lw!(atmos, top_at_1,
              lw_sources,
              emis_sfc,
              fluxes)
    else
      fluxes.flux_dn_dir = flux_dir[:,:]

      gas_optics!(k_dist, p_lay, p_lev,
                  t_lay,
                  gas_concs,
                  atmos,
                  toa_flux)
      delta_scale!(clouds)
      increment!(clouds, atmos)
      rte_sw!(atmos, top_at_1,
              mu0,   toa_flux,
              sfc_alb_dir, sfc_alb_dif,
              fluxes)
    end
  end


  # Compare with reference:
  ref_flux_up = Array{FT}(ds_ref[λ_string*"_flux_up"][:])
  ref_flux_dn = Array{FT}(ds_ref[λ_string*"_flux_dn"][:])

  diff_up = maximum( abs.( flux_up .- ref_flux_up ) )
  diff_dn = maximum( abs.( flux_dn .- ref_flux_dn ) )

  diff_up_ulps = maximum( abs.( flux_up .- ref_flux_up ) ./ eps.(ref_flux_up) )
  diff_dn_ulps = maximum( abs.( flux_dn .- ref_flux_dn ) ./ eps.(ref_flux_dn) )

  @show sqrt(1/eps(FT))
  @show diff_up, diff_up_ulps, maximum(abs.(ref_flux_up))
  @show diff_dn, diff_dn_ulps, maximum(abs.(ref_flux_dn))

  @test diff_up_ulps < sqrt(1/(1e6eps(FT)))
  @test diff_dn_ulps < sqrt(1/(1e6eps(FT)))

  close(ds_input)
  close(ds_ref)
  close(ds_k_dist)
  close(ds_cloud_optics)
end


@testset "All sky" begin
  datafolder = JRRTMGP.data_folder_rrtmgp()

  run_driver(datafolder; λ_string = "lw")
  run_driver(datafolder, λ_string = "sw")
end

