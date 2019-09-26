subroutine stop_on_err(msg)
  #
  # Print error message and stop
  #
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: msg
  if(len_trim(msg) > 0) then
    write (error_unit,*) trim(msg)
    write (error_unit,*) "test_adding stopping"
    stop
  end if
end subroutine
# ----------------------------------------------------------------------------------
program test_adding
  use mo_rte_kind,      only: wp, wl
  use mo_optical_props, only: ty_optical_props, ty_optical_props_arry
  use mo_rte_solver_kernels, &
                        only: adding

  use mo_test_files_io, only: is_sw, read_direction, read_spectral_disc, &
                              read_sw_solar_sources, read_sw_bc, read_lw_bc, &
                              read_two_stream, read_sources,  &
                              write_gpt_fluxes
  implicit none
  # ----------------------------------------------------------------------------------
  integer :: ncol, nlay, ngpt
  integer :: b, nBlocks, colS, colE
  integer, parameter :: blockSize = 8

  character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'

  class(ty_optical_props_arry), allocatable :: atmos

  real(wp), dimension(:,:,:), allocatable :: Rdif, Tdif, source_up, source_dn
  real(wp), dimension(:,  :), allocatable :: source_sfc
  real(wp), dimension(:,:,:), allocatable :: Rdir, Tdir, Tnoscat
  real(wp), dimension(:    ), allocatable :: mu0, tsi, t_sfc
  real(wp), dimension(:,  :), allocatable :: toa_src
  real(wp), dimension(:,  :), allocatable :: sfc_alb_dir, sfc_alb_dif, &
                                             sfc_alb_gpt
  real(wp)                                :: tsi_scaling

  real(wp), dimension(:,:,:), allocatable :: flux_up, flux_dn, flux_dn_dir

  logical :: top_at_1, do_sw
  integer :: i, j, k, ibnd, igpt
  type(ty_optical_props) :: spectral_disc
  # ----------------------------------------------------------------------------------
  do_sw = is_sw(fileName)
  call read_direction (fileName, top_at_1)
  call read_spectral_disc(fileName, spectral_disc)
  call read_sources(fileName, source_up, source_dn, source_sfc)
  if(do_sw) then
    call read_two_stream(fileName, Rdif, Tdif, Rdir, Tdir, Tnoscat)
    call read_sw_bc     (fileName, mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    mu0(:) = cos(mu0(:) * acos(-1._wp)/180.)
  else
    call read_two_stream(fileName, Rdif, Tdif)
    call read_lw_bc(fileName, t_sfc, sfc_alb_dif)
    # Read emissivity; convert to surface albedo
    sfc_alb_dif(:,:) = 1._wp - sfc_alb_dif(:,:)
  end if

  ncol = size(source_up, 1)
  nlay = size(source_up, 2)
  ngpt = size(source_up, 3)
  allocate(flux_up    (ncol,nlay+1,ngpt), &
           flux_dn    (ncol,nlay+1,ngpt))
  allocate(sfc_alb_gpt(ncol, ngpt))
  if(top_at_1) then
    flux_dn(:,     1,:)  = 0._wp
  else
    flux_dn(:,nlay+1,:)  = 0._wp
  end if
   # Expand surface albedos from bands to gpoints
  do igpt = 1, ngpt
    ibnd = spectral_disc%convert_gpt2band(igpt)
    sfc_alb_gpt(1:ncol,igpt) = sfc_alb_dif(ibnd,1:ncol)
  end do

  do k = 1, ngpt
    call adding(ncol, nlay, logical(top_at_1, wl), &
                sfc_alb_gpt(:,k),  &
                Rdif(:,:,k), Tdif(:,:,k),                            &
                source_dn(:,:,k), source_up(:,:,k), source_sfc(:,k), &
                flux_up(:,:,k), flux_dn(:,:,k))
  end do

  if(do_sw) then
    #
    # Compute direct beam for solar  - this is done in sources()
    #
    call read_sw_solar_sources(fileName, toa_src)
    allocate(flux_dn_dir(ncol, nlay+1, ngpt))

    if(top_at_1) then
      flux_dn_dir(:,     1,:)  = toa_src(:,:) * spread(mu0, 2, ngpt)
      do j = 1, nlay
        flux_dn_dir(:,j+1,:) = Tnoscat(:,j,:) * flux_dn_dir(:,j,:)
      end do
    else
      flux_dn_dir(:,nlay+1,:)  = toa_src(:,:) * spread(mu0, 2, ngpt)
      do j = nlay, 1, -1
        flux_dn_dir(:,j,:) = Tnoscat(:,j,:) * flux_dn_dir(:,j+1,:)
      end do
    end if
    #
    # Add direct beam so flux_dn is total
    #
    flux_dn(:,:,:) = flux_dn(:,:,:) + flux_dn_dir(:,:,:)
    call write_gpt_fluxes(fileName, flux_up, flux_dn, flux_dn_dir)
  else
    call write_gpt_fluxes(fileName, flux_up, flux_dn)
  end if

end program test_adding
