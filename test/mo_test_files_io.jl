####
#### Reads profiles needed for RRTMGP and related
#### calculations assuming a certain netCDF file layout
####

"""
    read_atmos(ds, FT, gas_names)

Read profiles for all columns
 - `p_lay` pressure (layers)
 - `t_lay` temperature (layers)
 - `p_lev` pressure (levels)
 - `t_lev` temperature (levels)
 - `col_dry` gas concentrations
"""
function read_atmos(ds, FT, gas_names)

  ncol = ds.dim["col"]
  nlay = ds.dim["lay"]
  nlev = ds.dim["lev"]
  @assert nlev == nlay+1

  p_lay = Array{FT}(ds["p_lay"][:])
  t_lay = Array{FT}(ds["t_lay"][:])
  p_lev = Array{FT}(ds["p_lev"][:])
  t_lev = Array{FT}(ds["t_lev"][:])

  available_gases = ["h2o", "co2", "o3", "n2o", "co",
                     "ch4", "o2", "n2", "ccl4", "cfc11",
                     "cfc12", "cfc22", "hfc143a", "hfc125", "hfc23",
                     "hfc32", "hfc134a", "cf4", "no2"]

  existing_gases = filter(ug->haskey(ds, "vmr_"*ug), available_gases)

  gsc = GasConcSize(ncol, nlay, (ncol, nlay), length(existing_gases))
  gas_concs = ty_gas_concs(FT, gas_names, ncol, nlay, gsc)

  for eg in existing_gases
    set_vmr!(gas_concs, eg, Array{FT}(read_field(ds, "vmr_"*eg)))
  end

  # col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
  col_dry = haskey(ds, "col_dry") ? Array{FT}(read_field(ds, "col_dry")) : nothing

  return p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry
end


#
# Does this file contain variables needed to do SW calculations ?
#
is_sw(ds) = haskey(ds,"solar_zenith_angle")

is_lw(ds) = !is_sw(ds)

#
# Read LW boundary conditions for all columns
#
read_lw_bc(ds) = ds["t_sfc"][:], ds["emis_sfc"][:]

#
# Read LW radiative transfer parameters
#
read_lw_rt(ds) = length( size( ds["angle"] ) )

#
# Read SW boundary conditions for all columns
#
function read_sw_bc(ds)#, sza, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
#    character(len=*),                      intent(in   ) :: fileName
#    real(FT), dimension(:),   allocatable, intent(inout) :: sza, tsi
#    real(FT), dimension(:,:), allocatable, intent(inout) :: sfc_alb_dir, sfc_alb_dif
#    real(FT),                              intent(inout) :: tsi_scaling

  sza         =  ds["solar_zenith_angle"][:]
  tsi         =  ds["total_solar_irradiance"][:]
  sfc_alb_dir =  ds["sfc_alb_direct"]
  sfc_alb_dif =  ds["sfc_alb_diffuse"]

  # read tsi_scaling only if variable is present in the netCDF file
  tsi_scaling = haskey(ds,"tsi_scaling") ? ds["tsi_scaling"][:] : nothing

  return sza,tsi,sfc_alb_dir,sfc_alb_dif,tsi_scaling

end

#
# Read spectral discretization
#
function read_spectral_disc!(ds, spectral_disc)
#    character(len=*),       intent(in   ) :: fileName
#    class(ty_optical_props), intent(inout) :: spectral_disc

#    integer :: ncid
#    integer :: nband
#    integer,  dimension(:,:), allocatable :: band_lims_gpt
#    real(FT), dimension(:,:), allocatable :: band_lims_wvn

  nband = ds.dim["band"]
  @assert ds.dim["pair"] == 2

  band_lims_wvn = ds["band_lims_wvn"][:]
  band_lims_gpt = ds["band_lims_gpt"][:]
  spectral_disc.base = ty_optical_props_base("read_spectral_disc!", band_lims_wvn, band_lims_gpt) # blank string used for now

end

#
# Read surface SW albedo and LW emissivity spectra from the surface test file
#
function read_sfc_test_file(ds)#, sfc_alb, sfc_emis)
#    character(len=*),           intent(in) :: fileName
#    real(FT), dimension(:,:), allocatable, intent(inout) :: sfc_alb  # Dimensions (nband,nspectra)
#    real(FT), dimension(:,:), allocatable, intent(inout) :: sfc_emis # Dimensions (nband,nspectra)

  @assert haskey(ds,"SW_albedo")
  @assert haskey(ds,"LW_emissivity")
  return ds["SW_albedo"][:], ds["LW_emissivity"][:]

end

#  #
#  # Paired procedures for reading and writing intermediate results used in unit testing
#  #
#  #--------------------------------------------------------------------------------------------------------------------
#  #
#  # Optical properties (tau, ssa and g or p if provided)
#  #
function read_optical_prop_values!(ds, opt_props)
#    character(len=*),                          intent(in ) :: fileName
#    class(ty_optical_props_arry), allocatable, intent(out) :: opt_props
#    # -------------------
#    integer :: ncid
#    integer :: ncol, nlay, ngpt, nmom, nband
#    real(FT), dimension(:,:), allocatable :: band_lims_wvn
#    integer,  dimension(:,:), allocatable :: band_lims_gpt
#    # -------------------

  @assert haskey(ds,"tau")

  ncol  = ds.dim["col"]
  nlay  = ds.dim["lay"]
  ngpt  = ds.dim["gpt"]
  nband = ds.dim["band"]

  @assert ds.dim["pair"] == 2

  band_lims_wvn = ds["band_lims_wvn"][:]
  band_lims_gpt = ds["band_lims_gpt"][:]

  opt_props.base = ty_optical_props_base("read_optical_prop_values!",band_lims_wvn, band_lims_gpt) #check for the initializing string

  if(haskey(ds, "p"))
    # n-stream calculation
    nmom = ds.dim["mom"]
    allocate!(opt_props)

    alloc_nstr!(opt_props,nmom, ncol, nlay)
    opt_props.ssa[:,:,:] = ds["ssa"][:] #,       ncol, nlay, ngpt)
    opt_props.p[:,:,:]   = ds["p"][:]   #,   nmom, ncol, nlay, ngpt)

  elseif (haskey(ds, "g"))
    # two-stream calculation
    # allocate(ty_optical_props_2str::opt_props)
    alloc_2str!(opt_props, ncol, nlay)
    opt_props.ssa[:,:,:] = ds["ssa"][:]
    opt_props.g[:,:,:]   = ds["g"][:]
  else
    # No scattering
    allocate(ty_optical_props_1scl::opt_props)
    alloc_1scl!(opt_props,ncol, nlay)
  end
  opt_props.tau[:,:,:] = ds["tau"][:] #, ncol, nlay, ngpt)

  #
  # Spectral discretization
  #
#    band_lims_wvn = read_field(ncid, 'band_lims_wvn', 2, nband)
#    band_lims_gpt = read_field(ncid, 'band_lims_gpt', 2, nband)

#    error(opt_props%init(band_lims_wvn, band_lims_gpt, read_string(ncid, 'name', 32)))
#    select type (opt_props)
#      class is (ty_optical_props_1scl)      # No scattering
#        error(opt_props%alloc_1scl(ncol, nlay))
#      class is (ty_optical_props_2str) # two-stream calculation
#        error(opt_props%alloc_2str(ncol, nlay))
#        opt_props%ssa = read_field(ncid, "ssa", ncol, nlay, ngpt)
#        opt_props%g   = read_field(ncid, "g",   ncol, nlay, ngpt)
#      class is (ty_optical_props_nstr) # n-stream calculation
#        error(opt_props%alloc_nstr(nmom, ncol, nlay))
#        opt_props%ssa = read_field(ncid, "ssa",       ncol, nlay, ngpt)
#        opt_props%p   = read_field(ncid, "p",   nmom, ncol, nlay, ngpt)
#    end select
#    opt_props%tau = read_field(ncid, "tau", ncol, nlay, ngpt)
end

#
# Which direction is up? Stored as a global attribute.
#

read_direction(ds) = (ds.attrib["top_at_1"] == 1)

#
# Sources of upward and downward diffuse radiation, for eah layer and at the surface
#

function read_sources(ds)#, source_up, source_dn, source_sfc)
#    character(len=*),           intent(in ) :: fileName
#    real(FT), dimension(:,:,:), allocatable, &
#                                intent(out) :: source_up, source_dn
#    real(FT), dimension(:,:  ), allocatable, &
#                                intent(out) :: source_sfc

  @assert haskey(ds,"source_up")

  source_up  = ds["source_up"][:]
  source_dn  = ds["source_dn"][:]
  source_sfc = ds["source_sfc"][:]

  return source_up, source_dn, source_sfc

end

#
# Longwave sources at layer centers; edges in two directions; surface
#   Also directionality since this will be needed for solution
#

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
  @assert ds.dim["pair"] == 2
  @assert haskey(ds,"lay_src")

  #
  # Spectral discretization
  #
  band_lims_wvn = ds["band_lims_wvn"][:]
  band_lims_gpt = ds["band_lims_gpt"][:]

  sources.optical_props.base = ty_optical_props_base("ty_source_func_lw", band_lims_wvn, band_lims_gpt)
  sources.optical_props.tau = Array{FT}(undef, ncol,nlay)

  sources.lay_source[:,:,:]     = ds["lay_src"][:]
  sources.lev_source_inc[:,:,:] = ds["lev_src_inc"][:]
  sources.lev_source_dec[:,:,:] = ds["lev_src_dec"][:]
  sources.sfc_source[:,:]       = ds["sfc_src"][:]

#    error(sources%init(band_lims_wvn, band_lims_gpt, read_string(ncid, 'name', 32)))
#    error(sources%alloc(ncol, nlay))
#    sources%lay_source     = read_field(ncid, 'lay_src',     ncol, nlay, ngpt)
#    sources%lev_source_inc = read_field(ncid, 'lev_src_inc', ncol, nlay, ngpt)
#    sources%lev_source_dec = read_field(ncid, 'lev_src_dec', ncol, nlay, ngpt)
#    sources%sfc_source     = read_field(ncid, 'sfc_src',     ncol,       ngpt)

end

#
# Shortwave source at TOA
#   Also directionality since this will be needed for solution
#

read_sw_solar_sources(ds) = ds["toa_src"][:]

#
# Two-stream results: reflection and transmission for diffuse and direct radiation; also extinction
#

function read_two_stream(ds)#, Rdif, Tdif, Rdir, Tdir, Tnoscat)
#    character(len=*),           intent(in ) :: fileName
#    real(FT), dimension(:,:,:), allocatable, &
#                                intent(out) :: Rdif, Tdif
#    real(FT), dimension(:,:,:), allocatable, optional, &
#                                intent(out) :: Rdir, Tdir, Tnoscat

  Rdif    = ds["Rdif"][:]
  Tdif    = ds["Tdif"][:]

  Rdir = haskey(ds,"Rdir") ? ds["Rdir"][:] : nothing
  Tdir = haskey(ds,"Tdir") ? ds["Tdir"][:] : nothing
  Tnoscat = haskey(ds,"Tnoscat") ? ds["Tnoscat"][:] : nothing

  return Rdif, Tdif, Rdir, Tdir, Tnoscat
end

#
# g-point fluxes
#

function read_gpt_fluxes(ds)# gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir)
#    character(len=*),           intent(in ) :: fileName
#    real(FT), dimension(:,:,:), allocatable, &
#                                intent(out) :: gpt_flux_up, gpt_flux_dn
#    real(FT), dimension(:,:,:), allocatable, optional, &
#                                intent(out) :: gpt_flux_dn_dir

  gpt_flux_up    = ds["gpt_flux_up"][:]
  gpt_flux_dn    = ds["gpt_flux_dn"][:]

  gpt_flux_dn_dir = haskey(ds,"gpt_flux_dn_dir") ? ds["gpt_flux_dn_dir"][:] : nothing

  return gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir
end
