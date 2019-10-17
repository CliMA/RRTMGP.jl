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
using JRRTMGP.mo_rte_sw
using JRRTMGP.mo_source_functions
using JRRTMGP.mo_test_files_io
using JRRTMGP.mo_fluxes
using JRRTMGP.mo_load_coefficients
using DataDeps
using JRRTMGP.mo_rfmip_io

function run_test()

  # use mo_optical_props, only: ty_optical_props, ty_optical_props_arry,
  #                             ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  # use mo_source_functions,
  #                       only: ty_source_func_lw
  # use mo_rte_solver_kernels,
  #                       only: sw_two_stream, sw_source_2str,
  #                             lw_two_stream, lw_source_2str, lw_combine_sources, lw_source_noscat,
  #                             apply_BC

  # use mo_test_files_io, only: read_optical_prop_values, read_spectral_disc,
  #                             is_sw, is_lw, read_direction,
  #                             read_sw_bc, read_sw_solar_sources,
  #                             read_lw_bc, read_lw_Planck_sources,
  #                             write_sources
  # implicit none
  # ----------------------------------------------------------------------------------
  # integer :: ncol, nlay, ngpt

  # character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'

  # class(ty_optical_props_arry), allocatable :: atmos

  # # SW-specific
  # real(wp), dimension(:    ), allocatable :: mu0, tsi
  # real(wp), dimension(:,  :), allocatable :: sfc_alb_dir, sfc_alb_dif, toa_src
  # real(wp), dimension(:,  :), allocatable :: Rdir, Tdir, Tnoscat
  # real(wp), dimension(:,:,:), allocatable :: flux_dn_dir
  # real(wp)                                :: tsi_scaling
  # # LW-specific
  # real(wp), dimension(:,:  ), allocatable :: gamma1, gamma2
  # real(wp), dimension(:,:  ), allocatable :: sfc_emis, sfc_src, emis_sfc_bnd
  # real(wp), dimension(:    ), allocatable :: t_sfc
  # type(ty_source_func_lw),    target      :: sources
  # real(wp), dimension(:,:,:), pointer     :: lev_src_up, lev_src_dn
  # # Generic
  # real(wp), dimension(:,  :), allocatable :: Rdif, Tdif
  # real(wp), dimension(:,:,:), allocatable :: source_dn, source_up
  # real(wp), dimension(:,  :), allocatable :: source_sfc
  # real(wp), dimension(:,  :), allocatable :: lev_source

  # type(ty_optical_props) :: spectral_disc
  # integer :: i, j, k
  # logical :: do_sw, do_lw, top_at_1
  # ----------------------------------------------------------------------------------

  clear_sky_dir = joinpath("..", "src", "examples","rfmip-clear-sky")
  data_dir = joinpath("..", "src", "rrtmgp","data")

  ENV["DATADEPS_ALWAYS_ACCEPT"] = true

  register(DataDep("rrtmgp-lw-inputs-outputs-clear",
                   "Reference files for test_two_stream sources tests (long-wave)",
                   "https://caltech.box.com/shared/static/j1mksrrq0bfhsdaeupi9vri2jtjfepfv.gz",
                   "5f43b57fd021703a557e98161944b012dcaded22a48bc3437848a43219fb1b72",
                   post_fetch_method=JRRTMGP.extract_targz))

  datafolder = datadep"rrtmgp-lw-inputs-outputs-clear"
  datafile = joinpath(datafolder,"test","sources","ref","rrtmgp-lw-inputs-outputs-clear.nc")
  ds = Dataset(datafile, "r")

  FT, I = Float64, Int

  atmos = read_optical_prop_values(ds, FT, I)
  ncol = get_ncol(atmos)
  nlay = get_nlay(atmos)
  ngpt = get_ngpt(atmos)
  source_up = Array{FT}(undef, ncol, nlay, ngpt)
  source_dn = Array{FT}(undef, ncol, nlay, ngpt)
  source_sfc = Array{FT}(undef, ncol, ngpt)
  Rdif = Array{FT}(undef, ncol,nlay)
  Tdif = Array{FT}(undef, ncol,nlay)

  spectral_disc = ty_optical_props_1scl(FT,I)
  read_spectral_disc!(ds, spectral_disc, FT, I)

  top_at_1 = read_direction(ds)
  do_sw = is_sw(ds)
  do_lw = !do_sw
  if do_sw
    flux_dn_dir = Array{FT}(undef, ncol, nlay+1, ngpt)
    Rdir = Array{FT}(undef, ncol,nlay)
    Tdir = Array{FT}(undef, ncol,nlay)
    Tnoscat = Array{FT}(undef, ncol,nlay)
    mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif = read_sw_bc(ds, FT)
    toa_src = read_sw_solar_sources(ds)
    mu0 .= cos.(mu0 .* acos(-FT(1))/FT(180))
    flux_dn_dir = apply_BC(ncol, nlay, ngpt, top_at_1, toa_src, mu0)
  else
    gamma1 = Array{FT}(undef, ncol,nlay)
    gamma2 = Array{FT}(undef, ncol,nlay)
    sfc_emis = Array{FT}(undef, ncol, ngpt)
    lev_source = Array{FT}(undef, ncol,nlay+1)
    sources = ty_source_func_lw(FT, I)
    read_lw_Planck_sources!(ds, sources)
    t_sfc, emis_sfc_bnd = read_lw_bc(ds)
    for k = 1:ngpt
      sfc_emis[1:ncol, k] = emis_sfc_bnd[convert_gpt2band(spectral_disc, k), :]
    end
  end

  if atmos isa ty_optical_props_1scl
    if do_lw
      if top_at_1
        lev_src_up = sources.lev_source_dec
        lev_src_dn = sources.lev_source_inc
      else
        lev_src_up = sources.lev_source_inc
        lev_src_dn = sources.lev_source_dec
      end
      for k = 1:ngpt
        lw_source_noscat!(ncol, nlay,
                          sources.lay_source[:,:,k], lev_src_up[:,:,k], lev_src_dn[:,:,k],
                          # Diffusivity angle would be applied in solver
                          atmos.tau[:,:,k]*FT(1.66), exp(-atmos.tau[:,:,k]*FT(1.66)),
                          @view(source_dn[:,:,k]), @view(source_up[:,:,k]))
        # Copied from kernels
        source_sfc[:, k] .= sfc_emis[:,k] .* sources.sfc_source[:,k]
      end
    else
      error("Haven't implemented direct beam source for SW")
    end

  elseif atmos isa ty_optical_props_2str
    if do_lw
      for k = 1:ngpt
        lw_two_stream!(ncol, nlay,
                       atmos.tau[:,:,k], atmos.ssa[:,:,k], atmos.g[:,:,k],
                       gamma1, gamma2, Rdif, Tdif)
       #
       # RRTMGP provides source functions at each level using the spectral mapping
       #   of each adjacent layer. Combine these for two-stream calculations
       #
       lw_combine_sources!(ncol, nlay, top_at_1,
                           sources.lev_source_inc[:,:,k], sources.lev_source_dec[:,:,k],
                           lev_source)
       #
       # Source function for diffuse radiation
       #
        lw_source_2str!(ncol, nlay, top_at_1,
                        sfc_emis[:,k], sources.sfc_source[:,k],
                        sources.lay_source[:,:,k], lev_source,
                        gamma1, gamma2, Rdif, Tdif, atmos.tau[:,:,k],
                        @view(source_dn[:,:,k]), @view(source_up[:,:,k]), @view(source_sfc[:,k]))
      end
    else
      for k = 1:ngpt
        sw_two_stream!(ncol, nlay, mu0,
                       atmos.tau[:,:,k], atmos.ssa[:,:,k], atmos.g[:,:,k],
                       Rdif, Tdif, Rdir, Tdir, Tnoscat)
        sw_source_2str!(ncol, nlay, top_at_1, Rdir, Tdir, Tnoscat,
                            sfc_alb_dir[convert_gpt2band(spectral_disc, k), :],
                            @view(source_up[:,:,k]), @view(source_dn[:,:,k]), @view(source_sfc[:,k]), @view(flux_dn_dir[:,:,k]))
      end
    end
  elseif atmos isa ty_optical_props_nstr
    error("Haven't implemented source calculations for multi-stream inputs")
  end

  # test_data(source_up, "source_up")
  # test_data(source_dn, "source_dn")

end

@testset "test_two_stream" begin
  run_test()
end

