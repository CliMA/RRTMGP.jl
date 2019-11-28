using RRTMGP
using NCDatasets
using RRTMGP.OpticalProps
using RRTMGP.RTESolver
using RRTMGP.FortranIntrinsics
using Test

get_array(ds, name, FT) = haskey(ds, name) ? convert(Array{FT}, ds[name][:]) : nothing
get_array(ds, name, FT, s) = haskey(ds, name) ? convert(Array{FT}, ds[name][:]) : zeros(s)
function read_sw_solar_sources(ds, FT)
  ncol  = ds.dim["col"]
  toa_src = get_array(ds, "toa_src", FT, (ncol))
  return toa_src
end

function read_sw_bc(ds, FT)
  ncol  = ds.dim["col"]
  nband = ds.dim["band"]
  μ_0         =  get_array(ds, "μ_0", FT, (ncol))
  tsi         =  get_array(ds, "tsi", FT, (ncol))
  sfc_alb_dir =  get_array(ds, "sfc_alb_dir", FT, (nband,ncol))
  sfc_alb_dif =  get_array(ds, "sfc_alb_dif", FT, (nband,ncol))

  tsi_scaling =  get_array(ds, "tsi_scaling", FT)
  return μ_0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif
end

function read_sources(ds, FT)
  source_up =  get_array(ds, "source_up", FT)
  source_dn =  get_array(ds, "source_dn", FT)
  source_sfc = get_array(ds, "source_sfc", FT)
  return source_up, source_dn, source_sfc
end

function read_two_stream(ds, FT)
  Rdif = get_array(ds, "Rdif", FT)
  Tdif = get_array(ds, "Tdif", FT)
  Rdir = get_array(ds, "Rdir", FT)
  Tdir = get_array(ds, "Tdir", FT)
  Tnoscat = get_array(ds, "Tnoscat", FT)
  return Rdif, Tdif, Rdir, Tdir, Tnoscat
end

function read_spectral_disc(ds, FT)
    # character(len=*),       intent(in   ) :: fileName
    # class(AbstractOpticalProps), intent(inout) :: spectral_disc

    # integer :: ncid
    # integer :: nband
    # integer,  dimension(:,:), allocatable :: band_lims_gpt
    # real(FT), dimension(:,:), allocatable :: band_lims_wvn

    band_lims_wvn = convert(Array{FT}, ds["band_lims_wvn"][:])
    band_lims_gpt = convert(Array{FT}, ds["band_lims_gpt"][:])
    op = OneScalar(FT,Int)
    init!(op, "spectral_disc", band_lims_wvn, band_lims_gpt)
    return op
end

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

  # class(AbstractOpticalPropsArry), allocatable :: atmos

  # real(FT), dimension(:,:,:), allocatable :: Rdif, Tdif, source_up, source_dn
  # real(FT), dimension(:,  :), allocatable :: source_sfc
  # real(FT), dimension(:,:,:), allocatable :: Rdir, Tdir, Tnoscat
  # real(FT), dimension(:    ), allocatable :: μ_0, tsi, t_sfc
  # real(FT), dimension(:,  :), allocatable :: toa_src
  # real(FT), dimension(:,  :), allocatable :: sfc_alb_dir, sfc_alb_dif,
  #                                            sfc_alb_gpt
  # real(FT)                                :: tsi_scaling

  # real(FT), dimension(:,:,:), allocatable :: flux_up, flux_dn, flux_dn_dir

  # logical :: top_at_1, do_sw
  # integer :: i, j, k, ibnd, igpt
  # type(AbstractOpticalProps) :: spectral_disc
  # # ----------------------------------------------------------------------------------
  ds = Dataset(fileName, "r")
  do_sw = haskey(ds, "solar_zenith_angle")
  top_at_1 = ds.attrib["top_at_1"]==1
  FT = Float64
  spectral_disc = read_spectral_disc(ds, FT)
  source_up, source_dn, source_sfc = read_sources(ds, FT)
  if do_sw
    Rdif, Tdif, Rdir, Tdir, Tnoscat = read_two_stream(ds, FT)
    μ_0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif = read_sw_bc(ds, FT)
    μ_0[:] .= cos.(μ_0[:] * acos(-FT(1))/FT(180))
  else
    # Rdif, Tdif = read_two_stream(fileName)
  #   t_sfc, sfc_alb_dif = read_lw_bc(fileName)
  #   FT = eltype(Rdif)
  #   # Read emissivity; convert to surface albedo
  #   sfc_alb_dif .= FT(1) .- sfc_alb_dif
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
      flux_dn_dir[:,    1,:]  .= toa_src .* spread(μ_0, 2, ngpt)
      for j = 1:nlay
        flux_dn_dir[:,j+1,:] .= Tnoscat[:,j,:] .* flux_dn_dir[:,j,:]
      end
    else
      flux_dn_dir[:,nlay+1,:] .= toa_src .* spread(μ_0, 2, ngpt)
      for j = nlay:-1:1
        flux_dn_dir[:,j,:] .= Tnoscat[:,j,:] .* flux_dn_dir[:,j+1,:]
      end
    end
    #
    # Add direct beam so flux_dn is total
    #
    flux_dn .+= flux_dn_dir
    # write_gpt_fluxes!(fileName, flux_up, flux_dn, flux_dn_dir)
  else
    # write_gpt_fluxes!(fileName, flux_up, flux_dn)
  end

  D = compare(ds, ((flux_dn,"gpt_flux_dn"),
                   (flux_up,"gpt_flux_up")))
  @show D
  close(ds)

end
