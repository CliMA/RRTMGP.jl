"""
    mo_rfmip_io

 This module reads an example file containing atomspheric conditions (temperature, pressure, gas concentrations)
   and surface properties (emissivity, temperature), defined on nlay layers across a set of ncol columns subject to
   nexp perturbations, and returns them in data structures suitable for use in rte and rrtmpg. The input data
   are partitioned into a user-specified number of blocks.
 For the moment only quantities relevant to longwave calculations are provided.

 The example files comes from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
   The protocol for this experiment allows for different specifications of which gases to consider:
 all gases, (CO2, CH4, N2O) + {CFC11eq; CFC12eq + HFC-134eq}. Ozone is always included
 The protocol does not specify the treatmet of gases like CO
"""
module mo_rfmip_io

using ..mo_gas_concentrations
using ..mo_util_string
using ..fortran_intrinsics
using NCDatasets

export read_kdist_gas_names, determine_gas_names, read_size, read_and_block_pt,
          read_and_block_sw_bc, read_and_block_lw_bc, read_and_block_gases_ty,
          unblock_and_write

# Find the size of the problem: columns, layers, perturbations (experiments)
function read_size(ds)
#    character(len=*),          intent(in   ) :: fileName
#    integer,         optional, intent(  out) :: ncol, nlay, nexp

  ncol = Int(ds.dim["site"])
  nlay = Int(ds.dim["layer"])
  nexp = Int(ds.dim["expt"])

  if ds.dim["level"] ≠ nlay+1
    error("read_size: number of levels should be nlay+1")
  end

  return ncol, nlay, nexp
end
#--------------------------------------------------------------------------------------------------------------------
#
# Return layer and level pressures and temperatures as arrays dimensioned (ncol, nlay/+1, nblocks)
#   Input arrays are dimensioned (nlay/+1, ncol, nexp)
#   Output arrays are allocated within this routine
#
function read_and_block_pt(ds, blocksize)
#    character(len=*),           intent(in   ) :: fileName
#    integer,                    intent(in   ) :: blocksize
#    real(FT), dimension(:,:,:), allocatable, & ! [blocksize, nlay/+1, nblocks]
#                                intent(  out) :: p_lay, p_lev, t_lay, t_lev
  # ---------------------------
#    integer :: ncid
#    integer :: b, nblocks
#    real(FT), dimension(:,:  ), allocatable :: temp2d
#    real(FT), dimension(:,:,:), allocatable :: temp3d
  # ---------------------------

  FT = Float64
  ncol_l, nlay_l, nexp_l = read_size(ds)

  if any([ncol_l, nlay_l, nexp_l]  .== 0)
    error("read_and_block_pt: Haven't read problem size yet.")
  end

  if (ncol_l*nexp_l)%blocksize ≠ 0
    error("read_and_block_pt: number of columns doesn't fit evenly into blocks.")
  end

  nblocks = Int((ncol_l*nexp_l)/blocksize)

  p_lay = Array{FT}(undef,blocksize,nlay_l  ,nblocks)
  t_lay = Array{FT}(undef,blocksize,nlay_l  ,nblocks)
  p_lev = Array{FT}(undef,blocksize,nlay_l+1,nblocks)
  t_lev = Array{FT}(undef,blocksize,nlay_l+1,nblocks)


  # Read p, T data; reshape to suit RRTMGP dimensions
  pres_layer = Array{FT}(ds["pres_layer"][:])
  temp3d = reshape( repeat(pres_layer,1,1,nexp_l), nlay_l, blocksize, nblocks )

  for b = 1:nblocks
    p_lay[:,:,b] = transpose(temp3d[:,:,b])
  end

  temp_layer = Array{FT}(ds["temp_layer"][:])
  temp3d = reshape(temp_layer, nlay_l, blocksize, nblocks )

  for b = 1:nblocks
    t_lay[:,:,b] = transpose(temp3d[:,:,b])
  end

  temp3d = []

  pres_level = Array{FT}(ds["pres_level"][:])
  temp3d = reshape( repeat(pres_level,1,1,nexp_l), nlay_l+1, blocksize, nblocks )

  for b = 1:nblocks
    p_lev[:,:,b] = transpose(temp3d[:,:,b])
  end

  temp_level = Array{FT}(ds["temp_level"][:])
  temp3d = reshape(temp_level,nlay_l+1, blocksize, nblocks )

  for b = 1:nblocks
    t_lev[:,:,b] = transpose(temp3d[:,:,b])
  end

  temp3d = []
  return p_lay, p_lev, t_lay, t_lev
end
#--------------------------------------------------------------------------------------------------------------------
#
# Read and reshape shortwave boundary conditions
#
function read_and_block_sw_bc(ds, blocksize)
#    character(len=*),           intent(in   ) :: fileName
#    integer,                    intent(in   ) :: blocksize
#    real(FT), dimension(:,:), allocatable, &
#                                intent(  out) :: surface_albedo, total_solar_irradiance, solar_zenith_angle
#    ! ---------------------------
#    integer :: ncid
#    integer :: nblocks
#    real(FT), dimension(ncol_l, nexp_l) :: temp2D
  # ---------------------------
  FT = Float64
  ncol_l, nlay_l, nexp_l = read_size(ds)

  if any([ncol_l, nlay_l, nexp_l] .== 0)
    error("read_and_block_sw_bc: Haven't read problem size yet.")
  end
  if (ncol_l*nexp_l)%blocksize ≠ 0
    error("read_and_block_sw_bc: number of columns doesn't fit evenly into blocks.")
  end
  nblocks = Int((ncol_l*nexp_l)/blocksize)
  #
  # Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
  #
  surface_albedo = Array{FT}(ds["surface_albedo"][:])
  temp2D = repeat(surface_albedo,1,nexp_l)

  surface_albedo = reshape(temp2D,blocksize,nblocks)

  total_solar_irradiance = Array{FT}(ds["total_solar_irradiance"][:])
  temp2D = repeat(total_solar_irradiance,1,nexp_l)
  total_solar_irradiance = reshape(temp2D, blocksize, nblocks)

  solar_zenith_angle = Array{FT}(ds["solar_zenith_angle"][:])
  temp2D = repeat(solar_zenith_angle,1,nexp_l)
  solar_zenith_angle = reshape(temp2D, blocksize, nblocks)
  return surface_albedo, total_solar_irradiance, solar_zenith_angle
end
#--------------------------------------------------------------------------------------------------------------------
#
# Read and reshape longwave boundary conditions
#
function read_and_block_lw_bc(ds, blocksize)
#    character(len=*),           intent(in   ) :: fileName
#    integer,                    intent(in   ) :: blocksize
#    real(FT), dimension(:,:), allocatable, &
#                                intent(  out) :: surface_emissivity, surface_temperature
#    ! ---------------------------
#    integer :: ncid
#    integer :: nblocks
#    real(FT), dimension(ncol_l, nexp_l) :: temp2D ! Required to make gfortran 8 work, not sure why
  FT = Float64
  # ---------------------------
  ncol_l, nlay_l, nexp_l = read_size(ds)
  if any([ncol_l, nlay_l, nexp_l]  .== 0)
    error("read_and_block_lw_bc: Haven't read problem size yet.")
  end

  if (ncol_l*nexp_l)%blocksize ≠ 0
    error("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
  end

  nblocks = convert(Int,(ncol_l*nexp_l)/blocksize)

  #
  # Allocate on assignment
  #
  temp2D = repeat( ds["surface_emissivity"][:] ,1,nexp_l)
  surface_emissivity  = Array{FT}(reshape(temp2D,blocksize,nblocks))

  surface_temperature = Array{FT}(reshape( ds["surface_temperature"][:], blocksize, nblocks )) # alternate version

  return surface_emissivity, surface_temperature
end
#--------------------------------------------------------------------------------------------------------------------
#
# Create a pair of string arrays - one containing the chemical name of each gas, used by the k-distribution, and
#   one containing the name as contained in the RFMIP input files - depending on the forcing scenario
# Forcing index (1 = all available greenhouse gases;
#                2 = CO2, CH4, N2O, CFC11eq
#                3 = CO2, CH4, N2O, CFC12eq, HFC-134eq
#                All scenarios use 3D values of ozone, water vapor so those aren't listed here
#
function determine_gas_names(ds, forcing_index)
#    character(len=*),                             intent(in   ) :: concentrationFile, kdistFile
#    integer,                                      intent(in   ) :: forcing_index
#    character(len=32), dimension(:), allocatable, intent(inout) :: names_in_kdist, names_in_file

    chem_name = ["co", "ch4", "o2", "n2o", "n2", "co2", "ccl4", "ch4", "ch3br", "ch3cl", "cfc22"]

    conc_name = ["carbon_monoxide", "methane", "oxygen",
        	   "nitrous_oxide", "nitrogen",
      	   "carbon_dioxide", "carbon_tetrachloride",
      	   "methane", "methyl_bromide",
      	   "methyl_chloride", "hcfc22"]
  if forcing_index == 1
    names_in_kdist = read_kdist_gas_names(ds)
    names_in_file = Array{String}(undef,length(names_in_kdist))
#      allocate(names_in_file(size(names_in_kdist)))
    for i = 1:length(names_in_kdist)
      names_in_file[i] = trim(lowercase(names_in_kdist[i]))
      #
      # Use a mapping between chemical formula and name if it exists
      #
      if names_in_file[i] in chem_name
        names_in_file[i] = conc_name[string_loc_in_array(names_in_file[i], chem_name)]
      end
    end
  elseif forcing_index == 2
    num_gases = 6
    #
    # Not part of the RFMIP specification, but oxygen is included because it's a major
    #    gas in some bands in the SW
    #
    names_in_kdist = ["co2", "ch4", "n2o", "o2", "cfc12", "cfc11"]
    names_in_file =  ["carbon_dioxide", "methane", "nitrous_oxide",
                      "oxygen", "cfc12", "cfc11eq"]
  elseif forcing_index == 3
    num_gases = 6
    #
    # Not part of the RFMIP specification, but oxygen is included because it's a major
    #    gas in some bands in the SW
    #
    names_in_kdist = ["co2", "ch4", "n2o", "o2", "cfc12", "hfc134a"]
    names_in_file =  ["carbon_dioxide", "methane", "nitrous_oxide",
                      "oxygen", "cfc12eq", "hfc134aeq"]
  else
    error("determine_gas_names: unknown value of forcing_index")
  end
  return names_in_kdist, names_in_file

end

# Read the names of the gases known to the k-distribution
function read_kdist_gas_names(ds)
  varName = "gas_names"
  ngases = ds.dim["absorber"]
  @assert haskey(ds,varName)
  kdist_gas_names = Array{String}(undef, ngases)
  temp_str = ds[varName][:]
  for i = 1:ngases
    kdist_gas_names[i] = trim(join(temp_str[:,i]))
  end
  return kdist_gas_names
end

# Read and reshape gas concentrations. RRTMGP requires gas concentrations to be supplied via a class
#   (ty_gas_concs). Gas concentrations are set via a call to gas_concs%set_vmr(name, values)
#   where `name` is nominally the chemical formula for the gas in question and `values` may be
#   a scalar, a 1-d profile assumed to apply to all columns, or an array of dimension (ncol, nlay).
# This routine outputs a vector nblocks long of these types so each element of the array can be passed to
#   the rrtmgp gas optics calculation in turn.
#
# This routine exploits RFMIP conventions: only water vapor and ozone vary by column within
#   each experiment.
# Fields in the RFMIP file have a trailing _GM (global mean); some fields use a chemical formula and other
#   a descriptive name, so a map is provided between these.
#

function read_and_block_gases_ty(ds, blocksize, gas_names, names_in_file)
#    character(len=*),           intent(in   ) :: fileName
#    integer,                    intent(in   ) :: blocksize
#    character(len=*),  dimension(:), &
#                                intent(in   ) :: gas_names ! Names used by the k-distribution/gas concentration type
#    character(len=*),  dimension(:), &
#                                intent(in   ) :: names_in_file ! Corresponding names in the RFMIP file
#    type(ty_gas_concs), dimension(:), allocatable, &
#                                intent(  out) :: gas_conc_array

  ncol_l, nlay_l, nexp_l = read_size(ds)

  if any([ncol_l, nlay_l, nexp_l] .== 0)
    error("read_and_block_lw_bc: Haven't read problem size yet.")
  end
  if (ncol_l*nexp_l)%blocksize ≠ 0
    error("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
  end
  nblocks = Int((ncol_l*nexp_l)/blocksize)
  FT = Float64
  gas_concs = ty_gas_concs(FT, gas_names, ncol_l, nlay_l)
  gas_conc_array = Vector([deepcopy(gas_concs) for i in 1:nblocks])

  # Experiment index for each colum
  # ORIGINAL: exp_num = reshape(spread([(b, b = 1, nexp_l)], 1, ncopies = ncol_l), shape = [blocksize, nblocks], order=[1,2])

  exp_num = freshape(spread(convert(Array,collect(1:nexp_l)'), 1, ncol_l), [blocksize, nblocks], order=[1,2])


  #
  # Water vapor and ozone depend on col, lay, exp: look just like other fields
  #
  # gas_conc_temp_3d = reshape(read_field(ncid, "water_vapor", nlay_l, ncol_l, nexp_l), &
  #                           shape = [nlay_l, blocksize, nblocks]) * read_scaling(ncid, "water_vapor")

  water_vapor = Array{FT}(ds["water_vapor"][:])
  gas_conc_temp_3d = reshape(water_vapor, nlay_l, blocksize, nblocks ) .* read_scaling(ds,"water_vapor")

  for b = 1:nblocks
    gas_conc_temp_3d_t = transpose(gas_conc_temp_3d[:,:,b])
    gas_conc_temp_3d_a = convert(Array, gas_conc_temp_3d_t)

    set_vmr!(gas_conc_array[b], "h2o", gas_conc_temp_3d_a, size(gas_conc_temp_3d_a))
  end

  ozone = Array{FT}(ds["ozone"][:])
  gas_conc_temp_3d = reshape(ozone, nlay_l, blocksize, nblocks ) * read_scaling(ds,"ozone")

  for b = 1:nblocks
    set_vmr!(gas_conc_array[b],"o3", convert(Array, transpose(gas_conc_temp_3d[:,:,b])))
  end

  #
  # All other gases are a function of experiment only
  #
  for g = 1:length(gas_names)
    #
    # Skip 3D fields above, also NO2 since RFMIP doesn't have this
    #
    if gas_names[g] in ["h2o", "o3", "no2"]
      continue
    end

    # Read the values as a function of experiment
   # gas_conc_temp_1d = read_field(ncid, trim(names_in_file(g)) // "_GM", nexp_l) * read_scaling(ncid, trim(names_in_file(g)) // "_GM")
    gas_conc_temp_1d = ds[trim(names_in_file[g]) * "_GM"][:] * read_scaling(ds,trim(names_in_file[g]) * "_GM")

    for b = 1:nblocks
      # Does every value in this block belong to the same experiment?
      if all( exp_num[2:end,b] .- exp_num[1,b] .== 0 )
        # Provide a scalar value
        set_vmr!(gas_conc_array[b],gas_names[g], gas_conc_temp_1d[exp_num[1,b]])
      else
        # Create 2D field, blocksize x nlay, with scalar values from each experiment
        set_vmr!(gas_conc_array[b], gas_names[g], spread(gas_conc_temp_1d[exp_num[:,b]], 2, nlay_l))
      end
    end

  end

  # for b = 1:nblocks
  #   for igas = 1:length(gas_conc_array[b].concs)
  #     test_data(gas_conc_array[b].concs[igas].conc, "gas_conc_b"*string(b)*"_i"*string(igas))
  #   end
  # end
  return gas_conc_array
end

function read_scaling(ds, varName)
  FT = Float64
  @assert haskey(ds,varName)
  scaling_str = ds[varName].attrib["units"]
  return parse(FT, scaling_str)
end

# Reshape and reorder values (nominally fluxes) from RTE order (ncol, nlev, nblocks)
#   to RFMIP order (nlev, ncol, nexp), then write them to a user-specified variable
#   in a netCDF file.
function unblock_and_write(ds, varName, values)
#    character(len=*),           intent(in   ) :: fileName, varName
#    real(FT), dimension(:,:,:),  & ! [blocksize, nlay/+1, nblocks]
#                                intent(in   ) :: values
  FT = Float64
  ncol_l, nlay_l, nexp_l = read_size(ds)

  @assert !any([ncol_l, nlay_l, nexp_l] .== 0)
  blocksize = size(values,1)
  nlev      = size(values,2)
  nblocks   = size(values,3)
  @assert nlev == nlay_l+1
  @assert (blocksize*nblocks) == (ncol_l*nexp_l)

  temp2D = Array{FT}(undef,nlev,ncol_l*nexp_l)

  for b = 1:nblocks
    temp2D[1:nlev, ((b-1)*blocksize+1):(b*blocksize)] = transpose(values[1:blocksize,1:nlev,b])
  end

  # Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
  # SK - commenting - assuming it is appending an existing file

  # this portion needs to be double checked
  defVar(ds,varName, reshape(temp2d,nlev,ncol_l,nexp_l), ("nlev","ncol_l","nexp_l"))

end

end
