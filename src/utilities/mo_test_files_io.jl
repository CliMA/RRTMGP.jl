module mo_test_files_io

using ..fortran_intrinsics
using ..mo_optical_props
using ..mo_source_functions
using ..mo_gas_concentrations
using ..mo_util_reorder
using ..mo_simple_netcdf

export read_atmos,
       is_lw,
       is_sw,
       read_lw_bc,
       read_sw_bc,
       read_lw_rt,
       read_spectral_disc!,
       read_sfc_test_file,
       read_optical_prop_values,
       read_direction,
       read_lw_Planck_sources,
       read_sw_solar_sources,
       read_two_stream,
       read_sources,
       read_lw_Planck_sources!,
       read_gpt_fluxes


#--------------------------------------------------------------------------------------------------------------------
#
# Read profiles for all columns  -- T, p, and gas concentrations
#   Allocation occurs on assignments (says the F2003 standard)
#
function read_atmos(ds)
#                        p_lay, t_lay, p_lev, t_lev,
#                        gas_concs, col_dry)
#    character(len=*),   intent(in   ) :: fileName
#    real(FT), dimension(:,:), allocatable,                 &
#                        intent(inout) :: p_lay, t_lay, p_lev, t_lev, col_dry
#    type(ty_gas_concs), intent(inout) :: gas_concs
#    # -------------------
#    integer :: ncid
#    integer :: ncol, nlay, nlev

  ncol = ds.dim["col"]
  nlay = ds.dim["lay"]
  nlev = ds.dim["lev"]
  if nlev ≠ nlay+1
    error("read_atmos: nlev should be nlay+1")
  end

  #
  # These lines assume that compilers follow the Fortran 2003 standard for
  #   allocating on assignment. This may require explicit compiler support
  #   e.g. -assume realloc_lhs flag for Intel
  #
  p_lay = ds["p_lay"][:]
  t_lay = ds["t_lay"][:]
  p_lev = ds["p_lev"][:]
  t_lev = ds["t_lev"][:]

  var_exists(ds, "vmr_h2o") && set_vmr!(gas_concs,"h2o", read_field(ds, "vmr_h2o"))
  var_exists(ds, "vmr_co2") && set_vmr!(gas_concs,"co2", read_field(ds, "vmr_co2"))
  var_exists(ds, "vmr_o3" ) && set_vmr!(gas_concs,"o3" , read_field(ds, "vmr_o3"))
  var_exists(ds, "vmr_n2o") && set_vmr!(gas_concs,"n2o", read_field(ds, "vmr_n2o"))
  var_exists(ds, "vmr_co" ) && set_vmr!(gas_concs,"co" , read_field(ds, "vmr_co"))
  var_exists(ds, "vmr_ch4") && set_vmr!(gas_concs,"ch4", read_field(ds, "vmr_ch4"))
  var_exists(ds, "vmr_o2" ) && set_vmr!(gas_concs,"o2" , read_field(ds, "vmr_o2"))
  var_exists(ds, "vmr_n2" ) && set_vmr!(gas_concs,"n2" , read_field(ds, "vmr_n2"))
  var_exists(ds, "vmr_ccl4" ) && set_vmr!(gas_concs,"ccl4" , read_field(ds, "vmr_ccl4"))
  var_exists(ds, "vmr_cfc11" ) && set_vmr!(gas_concs,"cfc11" , read_field(ds, "vmr_cfc11"))
  var_exists(ds, "vmr_cfc12" ) && set_vmr!(gas_concs,"cfc12" , read_field(ds, "vmr_cfc12"))
  var_exists(ds, "vmr_cfc22" ) && set_vmr!(gas_concs,"cfc22" , read_field(ds, "vmr_cfc22"))
  var_exists(ds, "vmr_hfc143a" ) && set_vmr!(gas_concs,"hfc143a" , read_field(ds, "vmr_hfc143a"))
  var_exists(ds, "vmr_hfc125" ) && set_vmr!(gas_concs,"hfc125" , read_field(ds, "vmr_hfc125"))
  var_exists(ds, "vmr_hfc23" ) && set_vmr!(gas_concs,"hfc23" , read_field(ds, "vmr_hfc23"))
  var_exists(ds, "vmr_hfc32" ) && set_vmr!(gas_concs,"hfc32" , read_field(ds, "vmr_hfc32"))
  var_exists(ds, "vmr_hfc134a" ) && set_vmr!(gas_concs,"hfc134a" , read_field(ds, "vmr_hfc134a"))
  var_exists(ds, "vmr_cf4" ) && set_vmr!(gas_concs,"cf4" , read_field(ds, "vmr_cf4"))
  var_exists(ds, "vmr_no2" ) && set_vmr!(gas_concs,"no2" , read_field(ds, "vmr_no2"))

  # col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
  var_exists(ds, "col_dry") && (col_dry = read_field(ds, "col_dry"))

  return p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry
end

#
# Write the atmospheric conditions that might be computed on the fly
#
function write_atmos(ds, t_lev, col_dry)
#    character(len=*),         intent(in) :: fileName
#    real(FT), dimension(:,:), intent(in) :: t_lev, col_dry

#    integer :: ncid, ncol, nlev, nlay

  #
  # At present these dimension sizes aren't used
  #   We could certainly check the array sizes against these dimension sizes
  #
  ds = Dataset(fileName,"a")

  ncol  = length( size( ds["col"] ) )
  nlay  = length( size( ds["lay"] ) )
  nlev  = length( size( ds["lev"] ) )

  defVar(ds,"col_dry", col_dry, ("col","lay"))
  defVar(ds,"t_lev",   t_lev,   ("col","lev"))

#    call create_var(ncid, "col_dry", ["col",  "lay"], [ncol, nlay])
#    call create_var(ncid, "t_lev",   ["col",  "lev"], [ncol, nlev])
#    error(write_field(ncid, "col_dry",  col_dry ))
#    error(write_field(ncid, "t_lev",     t_lev ))

end

#
# Does this file contain variables needed to do SW calculations ?
#
is_sw(ds) = haskey(ds,"solar_zenith_angle")
is_lw(ds) = !is_sw(ds)


#
# Read LW boundary conditions for all columns
#
function read_lw_bc(ds)
  t_sfc    =  ds["t_sfc"][:]
  emis_sfc =  ds["emis_sfc"][:]
  return t_sfc, emis_sfc
end

#
# Read LW radiative transfer parameters
#
read_lw_rt(ds) = length( size( ds["angle"] ) )

#--------------------------------------------------------------------------------------------------------------------
#
# Read SW boundary conditions for all columns
#
function read_sw_bc(ds, FT)
ncol  = get_dim_size(ds, "col")
nband = get_dim_size(ds, "band")
mu0         =  get_array(ds, "mu0", FT, (ncol))
tsi         =  get_array(ds, "tsi", FT, (ncol))
sfc_alb_dir =  get_array(ds, "sfc_alb_dir", FT, (nband,ncol))
sfc_alb_dif =  get_array(ds, "sfc_alb_dif", FT, (nband,ncol))

tsi_scaling =  haskey(ds, "tsi_scaling") ? get_array(ds, "tsi_scaling", FT) : nothing
return mu0, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif
end

#
# Read spectral discretization
#

function read_spectral_disc!(ds, spectral_disc, FT, I)
  # character(len=*),       intent(in   ) :: fileName
  # class(ty_optical_props), intent(inout) :: spectral_disc

  # integer :: ncid
  # integer :: nband
  # integer,  dimension(:,:), allocatable :: band_lims_gpt
  # real(FT), dimension(:,:), allocatable :: band_lims_wvn

  band_lims_wvn = convert(Array{FT}, ds["band_lims_wvn"][:])
  band_lims_gpt = convert(Array{I}, ds["band_lims_gpt"][:])
  init!(spectral_disc, "spectral_disc", band_lims_wvn, band_lims_gpt)
end

#
# Read surface SW albedo and LW emissivity spectra from the surface test file
#
function read_sfc_test_file(ds)#, sfc_alb, sfc_emis)
#    character(len=*),           intent(in) :: fileName
#    real(FT), dimension(:,:), allocatable, intent(inout) :: sfc_alb  # Dimensions (nband,nspectra)
#    real(FT), dimension(:,:), allocatable, intent(inout) :: sfc_emis # Dimensions (nband,nspectra)

#    integer :: ncid
#    integer :: nswband, nlwband, nspectra

  if !haskey(ds,"SW_albedo")
    error("read_sfc_test_file: file " * " doesn't contain SW_albedo field.")
  end
  nswband = ds.dim["swband"]
  nspectra = ds.dim["spectra"]
   # allocate on assignment
  sfc_alb = ds["SW_albedo"][:]

  if !haskey(ds,"LW_emissivity") #(.not. var_exists(ncid, 'LW_emissivity')) &
    error("read_sfc_test_file: file " * " doesn't contain LW_emissivity field.")
  end
  nlwband = ds.dim["lwband"]
  sfc_emis = ds["LW_emissivity"][:]

  return sfc_alb, sfc_emis

end

#
# Paired procedures for reading and writing intermediate results used in unit testing
#
function read_optical_prop_values(ds, FT, I)
# character(len=*),                          intent(in ) :: fileName
# class(ty_optical_props_arry), allocatable, intent(out) :: opt_props
# -------------------
#integer :: ncid
#integer :: ncol, nlay, ngpt, nmom, nband
#real(wp), dimension(:,:), allocatable :: band_lims_wvn
#integer,  dimension(:,:), allocatable :: band_lims_gpt
# -------------------

if !var_exists(ds, "tau")
  error("read_optical_prop_values: dataset doesn't contain tau field.")
end

ncol = get_dim_size(ds, "col")
nlay = get_dim_size(ds, "lay")
ngpt = get_dim_size(ds, "gpt")
nband = get_dim_size(ds, "band")

if(var_exists(ds, "p"))
  nmom = get_dim_size(ds, "mom")
  opt_props = ty_optical_props_nstr(FT, I)
elseif (var_exists(ds, "g"))
  opt_props = ty_optical_props_2str(FT, I)
else
  opt_props = ty_optical_props_1scl(FT, I)
end

#
# Spectral discretization
#
if (get_dim_size(ds, "pair") ≠ 2)
  error("read_optical_prop_values: pair dimension not 2 in file " * trim(fileName))
end
band_lims_wvn = Array{I}(read_field(ds, "band_lims_wvn"))
band_lims_gpt = Array{FT}(read_field(ds, "band_lims_gpt"))

name = read_string(ds, "name")
init!(opt_props, name, band_lims_wvn, band_lims_gpt)
if opt_props isa ty_optical_props_1scl      # No scattering
  alloc!(opt_props, ncol, nlay)
elseif opt_props isa ty_optical_props_2str # two-stream calculation
  alloc!(opt_props, ncol, nlay)
  opt_props.ssa = Array{FT}(read_field(ds, "ssa"))
  opt_props.g   = Array{FT}(read_field(ds, "g"))
elseif opt_props isa ty_optical_props_nstr # n-stream calculation
  alloc!(opt_props, nmom, ncol, nlay)
  opt_props.ssa = Array{FT}(read_field(ds, "ssa"))
  opt_props.p   = Array{FT}(read_field(ds, "p"))
else
  error("Bad opt_props in read_optical_prop_values")
end
opt_props.tau = Array{FT}(read_field(ds, "tau"))

return opt_props
end


function read_direction(ds)#, top_at_1)
#    character(len=*),           intent(in ) :: fileName
#    logical,                    intent(out) :: top_at_1
   # -------------------
#    integer :: ncid, status
#    integer :: top

#    status = nf90_get_att(ncid, NF90_GLOBAL, "top_at_1", top)
  return (ds.attrib["top_at_1"] == 1)

end #subroutine read_direction

function read_sources(ds, FT)
  if !haskey(ds,"source_up")
    error("read_sources: file doesn't contain source_up field.")
  end
  source_up =  get_array(ds, "source_up", FT)
  source_dn =  get_array(ds, "source_dn", FT)
  source_sfc = get_array(ds, "source_sfc", FT)
  return source_up, source_dn, source_sfc
end

function read_lw_Planck_sources!(ds, sources::ty_source_func_lw{FT}) where FT
#    character(len=*),        intent(in   ) :: fileName
#    type(ty_source_func_lw), intent(inout) :: sources
#    # -------------------
#    integer :: ncid
#    integer :: ncol, nlay, ngpt, nmom, nband
#    integer,  dimension(:,:), allocatable :: band_lims_gpt
#    real(FT), dimension(:,:), allocatable :: band_lims_wvn

  ncol  = ds.dim["col"]
  nlay  = ds.dim["lay"]
  ngpt  = ds.dim["gpt"]
  nband = ds.dim["band"]
  #
  # Error checking
  #
  if ds.dim["pair"] ≠ 2
    error("read_spectral_disc: pair dimension not 2 in file ")
  end
  if !haskey(ds,"lay_src")
    error("read_lw_Planck_sources: file " * " doesn't contain lay_src field.")
  end

  #
  # Spectral discretization
  #
  band_lims_wvn = get_array(ds, "band_lims_wvn", FT)
  band_lims_gpt = get_array(ds, "band_lims_gpt", FT)

  init!(sources, "Sources" ,band_lims_wvn, band_lims_gpt)
  alloc!(sources,ncol,nlay)

  sources.lay_source[:,:,:]     = ds["lay_src"][:]
  sources.lev_source_inc[:,:,:] = ds["lev_src_inc"][:]
  sources.lev_source_dec[:,:,:] = ds["lev_src_dec"][:]
  sources.sfc_source[:,:]       = ds["sfc_src"][:]
end

function read_sw_solar_sources(ds, FT)
if !haskey(ds,"toa_src")
  error("read_sw_solar_sources: file doesn't contain toa_src field.")
end
ncol  = get_dim_size(ds, "col")
ngpt = get_dim_size(ds, "gpt")
toa_src = get_array(ds, "toa_src", FT, (ncol))
return toa_src
end

function read_two_stream(ds)#, Rdif, Tdif, Rdir, Tdir, Tnoscat)
#    character(len=*),           intent(in ) :: fileName
#    real(FT), dimension(:,:,:), allocatable, &
#                                intent(out) :: Rdif, Tdif
#    real(FT), dimension(:,:,:), allocatable, optional, &
#                                intent(out) :: Rdir, Tdir, Tnoscat

  ncol = ds.dim["col"]
  nlay = ds.dim["lay"]
  ngpt = ds.dim["gpt"]

  Rdif    = ds["Rdif"][:]
  Tdif    = ds["Tdif"][:]

  Rdir = []
  Tdir = []
  Tnoscat = []


  if haskey(ds,"Rdir")
    Rdir = ds["Rdir"][:]
  end
  if haskey(ds,"Tdir")
    Tdir = ds["Tdir"][:]
  end
  if haskey(ds,"Tnoscat")
    Tnoscat = ds["Tnoscat"][:]
  end

  return Rdif, Tdif, Rdir, Tdir, Tnoscat
end

function read_gpt_fluxes(ds)
#    character(len=*),           intent(in ) :: fileName
#    real(FT), dimension(:,:,:), allocatable, &
#                                intent(out) :: gpt_flux_up, gpt_flux_dn
#    real(FT), dimension(:,:,:), allocatable, optional, &
#                                intent(out) :: gpt_flux_dn_dir
#    # -------------------
#    integer :: ncid
#    integer :: ncol, nlev, ngpt
#    # -------------------

  ncol = ds.dim["col"]
  nlev = ds.dim["lev"]
  ngpt = ds.dim["gpt"]

  gpt_flux_up    = ds["gpt_flux_up"][:]
  gpt_flux_dn    = ds["gpt_flux_dn"][:]

  gpt_flux_dn_dir = []
  if haskey(ds,"gpt_flux_dn_dir")
    gpt_flux_dn_dir = ds["gpt_flux_dn_dir"][:]
  end

  return gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir
end

end # module
