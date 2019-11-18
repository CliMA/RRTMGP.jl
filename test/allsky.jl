using Profile
using Test
using RRTMGP
using TimerOutputs
const to = TimerOutput()
using NCDatasets
using RRTMGP.mo_optical_props
using RRTMGP.mo_simple_netcdf
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
using RRTMGP.mo_cloud_optics
using RRTMGP.mo_rte_sw
using RRTMGP.mo_load_cloud_coefficients

include("mo_cloud_sampling.jl")
include("mo_test_files_io.jl")

function vmr_2d_to_1d!(gas_concs::ty_gas_concs{FT},
                       gas_concs_garand::ty_gas_concs{FT},
                       name::String,
                       sz1::I,
                       sz2::I) where {FT<:AbstractFloat,I<:Int}
  tmp = Array{FT}(undef, sz1, sz2)
  get_vmr!(tmp, gas_concs_garand, name)
  tmp_col = tmp[1, :]
  set_vmr!(gas_concs, name, tmp_col)
end

function all_sky(ds; use_luts=false, λ_string="", compile_first=false)
  @assert λ_string == "sw" || λ_string == "lw"
  k_dist_sym = Symbol(:k_dist,λ_string)
  cloud_optics_sym = Symbol(:cloud_optics,λ_string)

  # ----------------------------------------------------------------------------------
  # Variables
  # ----------------------------------------------------------------------------------
  # Arrays: dimensions (col, lay)
  # real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev
  # real(wp), dimension(:,:),   allocatable :: col_dry
  # real(wp), dimension(:,:),   allocatable :: temp_array

  #
  # Longwave only
  #
  # real(wp), dimension(:,:),   allocatable :: t_lev
  # real(wp), dimension(:),     allocatable :: t_sfc
  # real(wp), dimension(:,:),   allocatable :: emis_sfc # First dimension is band
  #
  # Shortwave only
  #
  # real(wp), dimension(:),     allocatable :: mu0
  # real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif # First dimension is band
  #
  # Source functions
  #
  #   Longwave
  # type(ty_source_func_lw), save               :: lw_sources
  #   Shortwave
  # real(wp), dimension(:,:), allocatable, save :: toa_flux
  #
  # Clouds
  #
  # real(wp), allocatable, dimension(:,:) :: lwp, iwp, rel, rei
  # logical,  allocatable, dimension(:,:) :: cloud_mask
  #
  # Output variables
  #
  # real(wp), dimension(:,:), target,
  #                           allocatable :: flux_up, flux_dn, flux_dir
  #
  # Derived types from the RTE and RRTMGP libraries
  #
  # type(ty_gas_optics_rrtmgp) :: k_dist
  # type(ty_cloud_optics)      :: cloud_optics_
  # type(ty_gas_concs)         :: gas_concs, gas_concs_garand, gas_concs_1col
  # class(ty_optical_props_arry),
  #                allocatable :: atmos, clouds
  # type(ty_fluxes_broadband)  :: fluxes

  #
  # Inputs to RRTMGP
  #
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
  #
  # Timing variables
  #
  # integer(kind=8)              :: start, finish, start_all, finish_all, clock_rate
  # real(wp)                     :: avg
  # integer(kind=8), allocatable :: elapsed(:)

  #
  # Parse command line for any file names, block size
  #
  # rrtmgp_clouds rrtmgp-clouds.nc $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc  128 1
  # rrtmgp_clouds rrtmgp-clouds.nc $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc  128 1

  gas_names = lowercase.(strip.(["h2o", "co2", "o3", "n2o", "co", "ch4", "o2", "n2"]))
  ngas = length(gas_names)
  nloops = 1
  ncol = 128
  write_fluxes = true

  #
  # Read temperature, pressure, gas concentrations.
  #   Arrays are allocated as they are read
  #
  FT = Float64
  I = Int64
  p_lay, t_lay, p_lev, t_lev, gas_concs_garand, col_dry = @timeit to "read_atmos" read_atmos(ds[:input], FT, gas_names)

  col_dry = nothing
  nlay = size(p_lay, 2)
  # For clouds we'll use the first column, repeated over and over
  gsc = GasConcSize(ncol, nlay, (ncol, nlay), ngas)
  gas_concs = ty_gas_concs(FT, gas_names, ncol, nlay, gsc)
  for igas = 1:ngas
    vmr_2d_to_1d!(gas_concs, gas_concs_garand, gas_names[igas], size(p_lay, 1), nlay)
  end

  #  If we trusted in Fortran allocate-on-assign we could skip the temp_array here
  temp_array = zeros(FT, ncol, nlay)
  temp_array .= spread_new(p_lay[1,:], 1, ncol)
  p_lay = deepcopy(temp_array)
  temp_array = zeros(FT, ncol, nlay)
  temp_array .= spread_new(t_lay[1,:], 1, ncol)
  t_lay = deepcopy(temp_array)
  temp_array = zeros(FT, ncol, nlay+1)
  temp_array .= spread_new(p_lev[1,:], 1, ncol)
  p_lev = deepcopy(temp_array)
  temp_array = zeros(FT, ncol, nlay+1)
  temp_array .= spread_new(t_lev[1,:], 1, ncol)
  t_lev = deepcopy(temp_array)

  # This puts pressure and temperature arrays on the GPU
  # load data into classes
  k_dist = load_and_init(ds[:k_dist], gas_concs)
  is_sw = source_is_external(k_dist)
  is_lw = !is_sw
  #
  # Should also try with Pade calculations
  if use_luts
    cloud_optics_ = ty_cloud_optics_new(FT,I)
    load_cld_lutcoeff!(cloud_optics_, ds[:cloud_optics])
  else
    cloud_optics_ = ty_cloud_optics(FT,I)
    load_cld_padecoeff!(cloud_optics_, ds[:cloud_optics])
  end
  set_ice_roughness!(cloud_optics_, 2)

  #
  # Problem sizes
  #
  nbnd = get_nband(k_dist.optical_props)
  ngpt = get_ngpt(k_dist.optical_props)
  top_at_1 = p_lay[1, 1] < p_lay[1, nlay]

  ps = ProblemSize(ncol, nlay, ngpt)

  # Clouds optical props are defined by band
  clouds_base = ty_optical_props_base("Clouds", get_band_lims_wavenumber(k_dist.optical_props))

  # LW calculations neglect scattering; SW calculations use the 2-stream approximation
  #   Here we choose the right variant of optical_props.
  if is_sw
    clouds = ty_optical_props_2str(clouds_base, ps)
    atmos = ty_optical_props_2str(k_dist.optical_props,ps)
  else
    clouds = ty_optical_props_1scl(clouds_base, ps)
    atmos = ty_optical_props_1scl(k_dist.optical_props, ps)
  end

  #
  # Allocate arrays for the optical properties themselves.
  #

  #  Boundary conditions depending on whether the k-distribution being supplied
  #   is LW or SW
  if is_sw
    # toa_flux is thread private
    toa_flux = zeros(FT, ncol, ngpt)
    #
    sfc_alb_dir = zeros(FT, nbnd, ncol)
    sfc_alb_dif = zeros(FT, nbnd, ncol)
    mu0 = zeros(FT, ncol)
    # Ocean-ish values for no particular reason
    sfc_alb_dir .= FT(0.06)
    sfc_alb_dif .= FT(0.06)
    mu0 .= FT(.86)
  else
    # lw_sorces is thread private
    lw_sources = ty_source_func_lw(ncol, nlay, k_dist.optical_props)

    t_sfc = zeros(FT, ncol)
    emis_sfc = zeros(FT, nbnd, ncol)
    # Surface temperature
    t_sfc .= t_lev[1, fmerge(nlay+1, 1, top_at_1)]
    emis_sfc .= FT(0.98)
  end

  #
  # Fluxes
  #
  flux_up = zeros(FT, ncol,nlay+1)
  flux_dn = zeros(FT, ncol,nlay+1)

  if is_sw
    flux_dir = zeros(FT, ncol,nlay+1)
  end
  #
  # Clouds
  #
  lwp = zeros(FT, ncol,nlay)
  iwp = zeros(FT, ncol,nlay)
  rel = zeros(FT, ncol,nlay)
  rei = zeros(FT, ncol,nlay)
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

  for iloop = 1:(compile_first ? 1 : nloops)
    cloud_optics!(cloud_optics_, lwp, iwp, rel, rei, clouds)
    #
    # Solvers
    #
    fluxes.flux_up = @view(flux_up[:,:])
    fluxes.flux_dn = @view(flux_dn[:,:])
    if is_lw
      gas_optics_int!(k_dist, p_lay, p_lev,
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

      gas_optics_ext!(k_dist, p_lay, p_lev,
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
  ref_flux_up = Array{FT}(ds[:ref][λ_string*"_flux_up"][:])
  ref_flux_dn = Array{FT}(ds[:ref][λ_string*"_flux_dn"][:])

  diff_up = maximum( abs.( flux_up .- ref_flux_up ) )
  diff_dn = maximum( abs.( flux_dn .- ref_flux_dn ) )

  diff_up_ulps = maximum( abs.( flux_up .- ref_flux_up ) ./ eps.(ref_flux_up) )
  diff_dn_ulps = maximum( abs.( flux_dn .- ref_flux_dn ) ./ eps.(ref_flux_dn) )

  # @show sqrt(1/eps(FT))
  # @show diff_up, diff_up_ulps, maximum(abs.(ref_flux_up))
  # @show diff_dn, diff_dn_ulps, maximum(abs.(ref_flux_dn))
  if !compile_first
    if use_luts
      @test diff_up_ulps < sqrt(1/(10eps(FT)))
      @test diff_dn_ulps < sqrt(1/(10eps(FT)))
    else
      # Need better test
      # @show flux_up
      # @show flux_dn
      # @show diff_up_ulps < sqrt(1/(eps(FT)))
      # @show diff_dn_ulps < sqrt(1/(eps(FT)))
    end
  end
  return nothing

end

