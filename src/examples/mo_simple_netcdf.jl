module mo_simple_netcdf

using ..fortran_intrinsics
using NCDatasets

export read_field,
       write_field,
       get_dim_size,
       get_array,
       var_exists

export read_spectral_disc,
       read_two_stream,
       read_sources,
       read_sw_bc,
       read_sw_solar_sources

function read_field(ds, varName, args...)
  @assert haskey(ds, varName)
  field = ds[varName][:]
  return field
end

function write_field(ds, varName, var)
  ndim = length( size(var) )
  if ndim == 1
    defVar(ds,varName, var, ("dim1"))
  elseif ndim == 2
    defVar(ds,varName, var, ("dim1","dim2"))
  elseif ndim == 3
    defVar(ds,varName, var, ("dim1","dim2","dim3"))
  elseif ndim == 4
    defVar(ds,varName, var, ("dim1","dim2","dim3","dim4"))
  else
    error("write_field: variables with more than 4 dimensions not supported at this point" * trim(varName))
  end
end

var_exists(ds, varName) = haskey(ds,varName)

get_array(ds, name, FT) = haskey(ds, name) ? convert(Array{FT}, ds[name][:]) : nothing
get_array(ds, name, FT, s) = haskey(ds, name) ? convert(Array{FT}, ds[name][:]) : zeros(s)
get_dim_size(ds, name) = ds.dim[name]

function read_sw_solar_sources(ds, FT)
  ncol  = get_dim_size(ds, "col")
  toa_src = get_array(ds, "toa_src", FT, (ncol))
  return toa_src
end

function read_sw_bc(ds, FT)
  ncol  = get_dim_size(ds, "col")
  nband = get_dim_size(ds, "band")
  mu0         =  get_array(ds, "mu0", FT, (ncol))
  tsi         =  get_array(ds, "tsi", FT, (ncol))
  sfc_alb_dir =  get_array(ds, "sfc_alb_dir", FT, (nband,ncol))
  sfc_alb_dif =  get_array(ds, "sfc_alb_dif", FT, (nband,ncol))

  tsi_scaling =  get_array(ds, "tsi_scaling", FT)
  return mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif
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
    # class(ty_optical_props), intent(inout) :: spectral_disc

    # integer :: ncid
    # integer :: nband
    # integer,  dimension(:,:), allocatable :: band_lims_gpt
    # real(FT), dimension(:,:), allocatable :: band_lims_wvn

    band_lims_wvn = convert(Array{FT}, ds["band_lims_wvn"][:])
    band_lims_gpt = convert(Array{FT}, ds["band_lims_gpt"][:])
    op = ty_optical_props_1scl(FT,Int)
    init!(op, "spectral_disc", band_lims_wvn, band_lims_gpt)
    return op
end


end # module