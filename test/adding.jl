using Test
using JRRTMGP

using mo_optical_props
using mo_rte_solver_kernels
using mo_test_files_io

function stop_on_err(msg)
  #
  # Print error message and stop
  #
  # character(len=*), intent(in) :: msg
  if len_trim(msg) > 0
    println(trim(msg))
    println("test_adding stopping")
    error("Done in stop_on_err in adding.jl")
  end
end
# ----------------------------------------------------------------------------------
# program test_adding
@testset "test_adding" begin
  # use mo_rte_kind,      only: wp, wl
  # use mo_optical_props, only: ty_optical_props, ty_optical_props_arry
  # use mo_rte_solver_kernels,
  #                       only: adding

  # use mo_test_files_io, only: is_sw, read_direction, read_spectral_disc,
  #                             read_sw_solar_sources, read_sw_bc, read_lw_bc,
  #                             read_two_stream, read_sources,
  #                             write_gpt_fluxes
  # implicit none
  # # ----------------------------------------------------------------------------------
  # integer :: ncol, nlay, ngpt
  # integer :: b, nBlocks, colS, colE
  # integer, parameter :: blockSize = 8
  blockSize = 8

  # character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'
  fileName = "rrtmgp-inputs-outputs.nc"

  # class(ty_optical_props_arry), allocatable :: atmos

  # real(wp), dimension(:,:,:), allocatable :: Rdif, Tdif, source_up, source_dn
  # real(wp), dimension(:,  :), allocatable :: source_sfc
  # real(wp), dimension(:,:,:), allocatable :: Rdir, Tdir, Tnoscat
  # real(wp), dimension(:    ), allocatable :: mu0, tsi, t_sfc
  # real(wp), dimension(:,  :), allocatable :: toa_src
  # real(wp), dimension(:,  :), allocatable :: sfc_alb_dir, sfc_alb_dif,
  #                                            sfc_alb_gpt
  # real(wp)                                :: tsi_scaling

  # real(wp), dimension(:,:,:), allocatable :: flux_up, flux_dn, flux_dn_dir

  # logical :: top_at_1, do_sw
  # integer :: i, j, k, ibnd, igpt
  # type(ty_optical_props) :: spectral_disc
  # # ----------------------------------------------------------------------------------
  do_sw = is_sw(fileName)
  top_at_1 = read_direction(fileName)
  spectral_disc = read_spectral_disc(fileName)
  source_up, source_dn, source_sfc = read_sources(fileName)
  if do_sw
    Rdif, Tdif, Rdir, Tdir, Tnoscat = read_two_stream(fileName)
    mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif = read_sw_bc(fileName)
    DT = eltype(Rdif)
    mu0[:] = cos(mu0[:] * acos(-DT(1))/DT(180))
  else
    Rdif, Tdif = read_two_stream(fileName)
    t_sfc, sfc_alb_dif = read_lw_bc(fileName)
    DT = eltype(Rdif)
    # Read emissivity; convert to surface albedo
    sfc_alb_dif[:,:] .= DT(1) .- sfc_alb_dif[:,:]
  end

  ncol = size(source_up, 1)
  nlay = size(source_up, 2)
  ngpt = size(source_up, 3)
  flux_up = Array{DT}(undef, ncol,nlay+1,ngpt)
  flux_dn = Array{DT}(undef, ncol,nlay+1,ngpt)
  sfc_alb_gpt = Array{DT}(ncol, ngpt)
  if top_at_1
    flux_dn[:,     1,:] = DT(0)
  else
    flux_dn[:,nlay+1,:] = DT(0)
  end
   # Expand surface albedos from bands to gpoints
  for igpt = 1:ngpt
    ibnd = spectral_disc.convert_gpt2band[igpt]
    sfc_alb_gpt[1:ncol,igpt] = sfc_alb_dif[ibnd,1:ncol]
  end

  for k = 1:ngpt
    adding!(ncol, nlay, logical[top_at_1, wl],
                sfc_alb_gpt[:,k],
                Rdif[:,:,k], Tdif[:,:,k],
                source_dn[:,:,k], source_up[:,:,k], source_sfc[:,k],
                flux_up[:,:,k], flux_dn[:,:,k])
  end

  if do_sw
    #
    # Compute direct beam for solar  - this is done in sources()
    #
    toa_src = read_sw_solar_sources(fileName)
    flux_dn_dir = Array{DT}(undef, ncol, nlay+1, ngpt)

    if top_at_1
      flux_dn_dir[:,    1,:]  .= toa_src[:,:] .* spread[mu0, 2, ngpt]
      for j = 1:nlay
        flux_dn_dir[:,j+1,:] .= Tnoscat[:,j,:] .* flux_dn_dir[:,j,:]
      end
    else
      flux_dn_dir[:,nlay+1,:] .= toa_src[:,:] .* spread[mu0, 2, ngpt]
      for j = nlay:-1:1
        flux_dn_dir[:,j,:] .= Tnoscat[:,j,:] .* flux_dn_dir[:,j+1,:]
      end
    end
    #
    # Add direct beam so flux_dn is total
    #
    flux_dn[:,:,:] .= flux_dn[:,:,:] .+ flux_dn_dir[:,:,:]
    write_gpt_fluxes!(fileName, flux_up, flux_dn, flux_dn_dir)
  else
    write_gpt_fluxes!(fileName, flux_up, flux_dn)
  end

end
