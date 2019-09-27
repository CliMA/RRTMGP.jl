using Test


get_array(ds, name, DT) = haskey(ds, name) ? convert(Array{DT}, ds[name][:]) : nothing
get_array(ds, name, DT, s) = haskey(ds, name) ? convert(Array{DT}, ds[name][:]) : zeros(s)
get_dim_size(ds, name) = ds.dim[name]

function read_sw_bc(ds, DT)
  ncol  = get_dim_size(ds, "col")
  nband = get_dim_size(ds, "band")
  mu0         =  get_array(ds, "mu0", DT, (ncol))
  tsi         =  get_array(ds, "tsi", DT, (ncol))
  sfc_alb_dir =  get_array(ds, "sfc_alb_dir", DT, (nband,ncol))
  sfc_alb_dif =  get_array(ds, "sfc_alb_dif", DT, (nband,ncol))

  tsi_scaling =  get_array(ds, "tsi_scaling", DT)
  return mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif
end

function read_sources(ds, DT)
  source_up =  get_array(ds, "source_up", DT)
  source_dn =  get_array(ds, "source_dn", DT)
  source_sfc = get_array(ds, "source_sfc", DT)
  return source_up, source_dn, source_sfc
end

function read_two_stream(ds, DT)
  Rdif = get_array(ds, "Rdif", DT)
  Tdif = get_array(ds, "Tdif", DT)
  Rdir = get_array(ds, "Rdir", DT)
  Tdir = get_array(ds, "Tdir", DT)
  Tnoscat = get_array(ds, "Tnoscat", DT)
  return Rdif, Tdif, Rdir, Tdir, Tnoscat
end

function read_spectral_disc(ds, DT)
    # character(len=*),       intent(in   ) :: fileName
    # class(ty_optical_props), intent(inout) :: spectral_disc

    # integer :: ncid
    # integer :: nband
    # integer,  dimension(:,:), allocatable :: band_lims_gpt
    # real(wp), dimension(:,:), allocatable :: band_lims_wvn

    # ! -------------------
    # if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
    #   call stop_on_err("read_spectral_disc: can't open file " // trim(fileName))

    # nband = length(size(ds["band"]))
    # nband = ds.dim["band"]
    # if (get_dim_size(ncid, 'pair') /= 2) &
    #   call stop_on_err("read_spectral_disc: pair dimension not 2 in file "//trim(fileName) )
    # @show ds["pair"][:]
    # @show ds["name"][:]

    band_lims_wvn = convert(Array{DT}, ds["band_lims_wvn"][:])
    band_lims_gpt = convert(Array{DT}, ds["band_lims_gpt"][:])
    return ty_optical_props_1scl("spectral_disc", band_lims_wvn, band_lims_gpt)

    # ncid = nf90_close(ncid)
end

using JRRTMGP
using NCDatasets
using JRRTMGP.mo_optical_props
using JRRTMGP.mo_rte_solver_kernels
# using JRRTMGP.mo_test_files_io

# function stop_on_err(msg)
#   #
#   # Print error message and stop
#   #
#   # character(len=*), intent(in) :: msg
#   if len_trim(msg) > 0
#     println(trim(msg))
#     println("test_adding stopping")
#     error("Done in stop_on_err in adding.jl")
#   end
# end
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
  data_dir = "../../rte-rrtmgp-test-val/test/adding/ref/"
  fileName = data_dir*"rrtmgp-sw-inputs-outputs-clear.nc"

  # data_dir = "../rte-rrtmgp-test-val/test/adding/ref/"
  # fileName = data_dir*"rrtmgp-sw-inputs-outputs-clear.nc"

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
  ds = Dataset(fileName, "r")
  do_sw = haskey(ds, "solar_zenith_angle")
  # do_sw = var_exists(ncid, '')
  @show do_sw
  # @show haskey(ds, "top_at_1")
  # top_at_1 = read_direction(fileName)
  top_at_1 = ds.attrib["top_at_1"]==1
  # @show top_at_1
  DT = Float64
  spectral_disc = read_spectral_disc(ds, DT)
  @show spectral_disc
  source_up, source_dn, source_sfc = read_sources(ds, DT)
  if do_sw
    Rdif, Tdif, Rdir, Tdir, Tnoscat = read_two_stream(ds, DT)
    mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif = read_sw_bc(ds, DT)
    mu0[:] .= cos.(mu0[:] * acos(-DT(1))/DT(180))
  else
  #   Rdif, Tdif = read_two_stream(fileName)
  #   t_sfc, sfc_alb_dif = read_lw_bc(fileName)
  #   DT = eltype(Rdif)
  #   # Read emissivity; convert to surface albedo
  #   sfc_alb_dif[:,:] .= DT(1) .- sfc_alb_dif[:,:]
  end

  # ncol = size(source_up, 1)
  # nlay = size(source_up, 2)
  # ngpt = size(source_up, 3)
  # flux_up = Array{DT}(undef, ncol,nlay+1,ngpt)
  # flux_dn = Array{DT}(undef, ncol,nlay+1,ngpt)
  # sfc_alb_gpt = Array{DT}(ncol, ngpt)
  # if top_at_1
  #   flux_dn[:,     1,:] .= DT(0)
  # else
  #   flux_dn[:,nlay+1,:] .= DT(0)
  # end
  #  # Expand surface albedos from bands to gpoints
  # for igpt = 1:ngpt
  #   ibnd = spectral_disc.convert_gpt2band[igpt]
  #   sfc_alb_gpt[1:ncol,igpt] = sfc_alb_dif[ibnd,1:ncol]
  # end

  # for k = 1:ngpt
  #   adding!(ncol, nlay, logical[top_at_1, wl],
  #               sfc_alb_gpt[:,k],
  #               Rdif[:,:,k], Tdif[:,:,k],
  #               source_dn[:,:,k], source_up[:,:,k], source_sfc[:,k],
  #               flux_up[:,:,k], flux_dn[:,:,k])
  # end

  # if do_sw
  #   #
  #   # Compute direct beam for solar  - this is done in sources()
  #   #
  #   toa_src = read_sw_solar_sources(fileName)
  #   flux_dn_dir = Array{DT}(undef, ncol, nlay+1, ngpt)

  #   if top_at_1
  #     flux_dn_dir[:,    1,:]  .= toa_src[:,:] .* spread[mu0, 2, ngpt]
  #     for j = 1:nlay
  #       flux_dn_dir[:,j+1,:] .= Tnoscat[:,j,:] .* flux_dn_dir[:,j,:]
  #     end
  #   else
  #     flux_dn_dir[:,nlay+1,:] .= toa_src[:,:] .* spread[mu0, 2, ngpt]
  #     for j = nlay:-1:1
  #       flux_dn_dir[:,j,:] .= Tnoscat[:,j,:] .* flux_dn_dir[:,j+1,:]
  #     end
  #   end
  #   #
  #   # Add direct beam so flux_dn is total
  #   #
  #   flux_dn[:,:,:] .= flux_dn[:,:,:] .+ flux_dn_dir[:,:,:]
  #   # write_gpt_fluxes!(fileName, flux_up, flux_dn, flux_dn_dir)
  # else
  #   # write_gpt_fluxes!(fileName, flux_up, flux_dn)
  # end
  close(ds)

end
