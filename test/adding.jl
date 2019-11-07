using JRRTMGP
using NCDatasets
using JRRTMGP.mo_optical_props
using JRRTMGP.mo_rte_solver_kernels
using JRRTMGP.mo_simple_netcdf
using JRRTMGP.fortran_intrinsics
using Test

function compare(ds, var_tuple)
  D = Dict()
  for (A, name) in var_tuple
    A_ds = get_array(ds, name, eltype(A))
    d = abs.(A_ds - A)
    # @show size(A_ds), length(size(d))
    @show name, max(d...)
    # D[name] = [sum(d, dims=i) for i in 1:length(size(d))]
    # D[name] = [sum(d, dims=[j for j in 1:length(size(d)) if i ≠ j]) for i in 1:length(size(d))]
    D[name] = false
  end
  return D
end


@testset "test_adding" begin
  # integer :: ncol, nlay, ngpt
  # integer :: b, nBlocks, colS, colE
  # integer, parameter :: blockSize = 8
  blockSize = 8

  # character(len=128) :: fileName = 'rrtmgp-inputs-outputs.nc'
  data_dir = "../../rte-rrtmgp-test-val/test/adding/ref/"
  fileName = data_dir*"rrtmgp-sw-inputs-outputs-clear.nc"

  # data_dir = "../rte-rrtmgp-test-val/test/adding/ref/"
  # fileName = data_dir*"rrtmgp-sw-inputs-outputs-clear.nc"

  # class(ty_optical_props_arry), allocatable :: atmos

  # real(FT), dimension(:,:,:), allocatable :: Rdif, Tdif, source_up, source_dn
  # real(FT), dimension(:,  :), allocatable :: source_sfc
  # real(FT), dimension(:,:,:), allocatable :: Rdir, Tdir, Tnoscat
  # real(FT), dimension(:    ), allocatable :: mu0, tsi, t_sfc
  # real(FT), dimension(:,  :), allocatable :: toa_src
  # real(FT), dimension(:,  :), allocatable :: sfc_alb_dir, sfc_alb_dif,
  #                                            sfc_alb_gpt
  # real(FT)                                :: tsi_scaling

  # real(FT), dimension(:,:,:), allocatable :: flux_up, flux_dn, flux_dn_dir

  # logical :: top_at_1, do_sw
  # integer :: i, j, k, ibnd, igpt
  # type(ty_optical_props) :: spectral_disc
  # # ----------------------------------------------------------------------------------
  ds = Dataset(fileName, "r")
  do_sw = haskey(ds, "solar_zenith_angle")
  # do_sw = var_exists(ncid, '')
  # top_at_1 = read_direction(fileName)
  top_at_1 = ds.attrib["top_at_1"]==1
  FT = Float64
  spectral_disc = read_spectral_disc(ds, FT)
  source_up, source_dn, source_sfc = read_sources(ds, FT)
  if do_sw
    Rdif, Tdif, Rdir, Tdir, Tnoscat = read_two_stream(ds, FT)
    mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif = read_sw_bc(ds, FT)
    mu0[:] .= cos.(mu0[:] * acos(-FT(1))/FT(180))
  else
    # Rdif, Tdif = read_two_stream(fileName)
  #   t_sfc, sfc_alb_dif = read_lw_bc(fileName)
  #   FT = eltype(Rdif)
  #   # Read emissivity; convert to surface albedo
  #   sfc_alb_dif[:,:] .= FT(1) .- sfc_alb_dif[:,:]
  end

  ncol = size(source_up, 1)
  nlay = size(source_up, 2)
  ngpt = size(source_up, 3)
  flux_up = Array{FT}(undef, ncol,nlay+1,ngpt)
  fill!(flux_up, 0)
  flux_dn = Array{FT}(undef, ncol,nlay+1,ngpt)
  fill!(flux_dn, 0)
  sfc_alb_gpt = Array{FT}(undef, ncol, ngpt)
  if top_at_1
    flux_dn[:,     1,:] .= FT(0)
  else
    flux_dn[:,nlay+1,:] .= FT(0)
  end
  # Expand surface albedos from bands to gpoints
  for igpt = 1:ngpt
    ibnd = convert_gpt2band(spectral_disc, igpt)
    sfc_alb_gpt[1:ncol,igpt] = sfc_alb_dif[ibnd,1:ncol]
  end

  for k = 1:ngpt
    adding!(ncol, nlay, top_at_1,
                sfc_alb_gpt[:,k],
                Rdif[:,:,k], Tdif[:,:,k],
                @view(source_dn[:,:,k]), @view(source_up[:,:,k]), @view(source_sfc[:,k]),
                flux_up[:,:,k], flux_dn[:,:,k])
  end

  if do_sw
    #
    # Compute direct beam for solar  - this is done in sources()
    #
    toa_src = read_sw_solar_sources(ds, FT)
    flux_dn_dir = Array{FT}(undef, ncol, nlay+1, ngpt)

    if top_at_1
      flux_dn_dir[:,    1,:]  .= toa_src[:,:] .* spread(mu0, 2, ngpt)
      for j = 1:nlay
        flux_dn_dir[:,j+1,:] .= Tnoscat[:,j,:] .* flux_dn_dir[:,j,:]
      end
    else
      flux_dn_dir[:,nlay+1,:] .= toa_src[:,:] .* spread(mu0, 2, ngpt)
      for j = nlay:-1:1
        flux_dn_dir[:,j,:] .= Tnoscat[:,j,:] .* flux_dn_dir[:,j+1,:]
      end
    end
    #
    # Add direct beam so flux_dn is total
    #
    flux_dn[:,:,:] .= flux_dn[:,:,:] .+ flux_dn_dir[:,:,:]
    # write_gpt_fluxes!(fileName, flux_up, flux_dn, flux_dn_dir)
  else
    # write_gpt_fluxes!(fileName, flux_up, flux_dn)
  end

  D = compare(ds, ((flux_dn,"gpt_flux_dn"),
                   (flux_up,"gpt_flux_up")))
  @show D
  close(ds)

end
