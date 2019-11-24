"""
    RFMIPIO

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
module RFMIPIO

using ..GasConcentrations
using ..Utilities
using ..FortranIntrinsics
using NCDatasets

export read_kdist_gas_names, determine_gas_names, read_size, read_and_block_pt,
          read_and_block_sw_bc, read_and_block_lw_bc, read_and_block_gases_ty,
          unblock_and_write

"""
    read_size(ds)

Find the size of the problem: columns, layers, perturbations (experiments)
"""
function read_size(ds)
  ncol = Int(ds.dim["site"])
  nlay = Int(ds.dim["layer"])
  nexp = Int(ds.dim["expt"])

  @assert ds.dim["level"] == nlay+1
  return ncol, nlay, nexp
end

"""
    read_and_block_pt(ds, blocksize)

Return layer and level pressures and temperatures as arrays dimensioned (ncol, nlay/+1, nblocks)
   Input arrays are dimensioned (nlay/+1, ncol, nexp)
   Output arrays are allocated within this routine
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    real(FT), dimension(:,:,:), allocatable, & ! [blocksize, nlay/+1, nblocks]
                                intent(  out) :: p_lay, p_lev, t_lay, t_lev
"""
function read_and_block_pt(ds, blocksize)

  FT = Float64
  ncol_l, nlay_l, nexp_l = read_size(ds)

  @assert !any([ncol_l, nlay_l, nexp_l]  .== 0)
  @assert (ncol_l*nexp_l)%blocksize == 0

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

  return p_lay, p_lev, t_lay, t_lev
end

"""
    read_and_block_sw_bc(ds, blocksize)

Read and reshape shortwave boundary conditions
    character(len=*),           intent(in   ) :: fileName
    integer,                    intent(in   ) :: blocksize
    real(FT), dimension(:,:), allocatable, &
                                intent(  out) :: surface_albedo, total_solar_irradiance, solar_zenith_angle
"""
function read_and_block_sw_bc(ds, blocksize)
  FT = Float64
  ncol_l, nlay_l, nexp_l = read_size(ds)

  @assert !any([ncol_l, nlay_l, nexp_l] .== 0)
  @assert (ncol_l*nexp_l)%blocksize == 0
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

"""
    read_and_block_lw_bc(ds, blocksize)

Read and reshape longwave boundary conditions


character(len=*),           intent(in   ) :: fileName
integer,                    intent(in   ) :: blocksize
real(FT), dimension(:,:), allocatable, &
                            intent(  out) :: surface_emissivity, surface_temperature
"""
function read_and_block_lw_bc(ds, blocksize)
  FT = Float64

  ncol_l, nlay_l, nexp_l = read_size(ds)
  @assert !any([ncol_l, nlay_l, nexp_l]  .== 0)

  @assert (ncol_l*nexp_l)%blocksize == 0

  nblocks = convert(Int,(ncol_l*nexp_l)/blocksize)

  # Allocate on assignment
  temp2D = repeat( ds["surface_emissivity"][:] ,1,nexp_l)
  surface_emissivity  = Array{FT}(reshape(temp2D,blocksize,nblocks))

  surface_temperature = Array{FT}(reshape( ds["surface_temperature"][:], blocksize, nblocks )) # alternate version

  return surface_emissivity, surface_temperature
end

"""
    determine_gas_names(ds, forcing_index)

Create a pair of string arrays - one containing the chemical name of each gas, used by the k-distribution, and
 one containing the name as contained in the RFMIP input files - depending on the forcing scenario
Forcing index (1 = all available greenhouse gases;
              2 = CO2, CH4, N2O, CFC11eq
              3 = CO2, CH4, N2O, CFC12eq, HFC-134eq
              All scenarios use 3D values of ozone, water vapor so those aren't listed here

  character(len=*),                             intent(in   ) :: concentrationFile, kdistFile
  integer,                                      intent(in   ) :: forcing_index
  character(len=32), dimension(:), allocatable, intent(inout) :: names_in_kdist, names_in_file
"""
function determine_gas_names(ds, forcing_index)

  chem_name = ["co", "ch4", "o2", "n2o", "n2", "co2", "ccl4", "ch4", "ch3br", "ch3cl", "cfc22"]

  conc_name = ["carbon_monoxide",
               "methane",
               "oxygen",
               "nitrous_oxide",
               "nitrogen",
               "carbon_dioxide",
               "carbon_tetrachloride",
               "methane",
               "methyl_bromide",
               "methyl_chloride",
               "hcfc22"]

  @assert any(forcing_index .== [1,2,3])
  if forcing_index == 1
    names_in_kdist = read_kdist_gas_names(ds)

    # Use a mapping between chemical formula and name if it exists
    names_in_file = map(x->x in chem_name ?
      conc_name[loc_in_array(x,chem_name)] : x, names_in_kdist)

  elseif forcing_index == 2

    # Not part of the RFMIP specification, but oxygen is included because it's a major
    #    gas in some bands in the SW
    names_in_kdist = ["co2", "ch4", "n2o", "o2", "cfc12", "cfc11"]
    names_in_file =  ["carbon_dioxide", "methane", "nitrous_oxide",
                      "oxygen", "cfc12", "cfc11eq"]
  elseif forcing_index == 3

    # Not part of the RFMIP specification, but oxygen is included because it's a major
    #    gas in some bands in the SW
    names_in_kdist = ["co2", "ch4", "n2o", "o2", "cfc12", "hfc134a"]
    names_in_file =  ["carbon_dioxide", "methane", "nitrous_oxide",
                      "oxygen", "cfc12eq", "hfc134aeq"]
  end
  return names_in_kdist, names_in_file

end

# Read the names of the gases known to the k-distribution
read_kdist_gas_names(ds) =
  lowercase.(strip.( String[join(ds["gas_names"][:][:,i]) for i = 1:ds.dim["absorber"]] ))

"""
    read_and_block_gases_ty(ds, blocksize, gas_names, names_in_file)

Read and reshape gas concentrations. RRTMGP requires gas concentrations to be supplied via a class
 (GasConcs). Gas concentrations are set via a call to gas_concs%set_vmr(name, values)
 where `name` is nominally the chemical formula for the gas in question and `values` may be
 a scalar, a 1-d profile assumed to apply to all columns, or an array of dimension (ncol, nlay).
This routine outputs a vector nblocks long of these types so each element of the array can be passed to
 the rrtmgp gas optics calculation in turn.

This routine exploits RFMIP conventions: only water vapor and ozone vary by column within
 each experiment.
Fields in the RFMIP file have a trailing _GM (global mean); some fields use a chemical formula and other
 a descriptive name, so a map is provided between these.

  character(len=*),           intent(in   ) :: fileName
  integer,                    intent(in   ) :: blocksize
  character(len=*),  dimension(:), &
                              intent(in   ) :: gas_names ! Names used by the k-distribution/gas concentration type
  character(len=*),  dimension(:), &
                              intent(in   ) :: names_in_file ! Corresponding names in the RFMIP file
  type(GasConcs), dimension(:), allocatable, &
                              intent(  out) :: gas_conc_array
"""
function read_and_block_gases_ty(ds, blocksize, gas_names, names_in_file)
  ncol_l, nlay_l, nexp_l = read_size(ds)
  @assert !any([ncol_l, nlay_l, nexp_l] .== 0)
  @assert (ncol_l*nexp_l)%blocksize == 0
  nblocks = Int((ncol_l*nexp_l)/blocksize)
  FT = Float64
  gsc = GasConcSize(ncol_l, nlay_l, (blocksize, nlay_l), length(gas_names))
  gas_concs = GasConcs(FT, gas_names, ncol_l, nlay_l, gsc)
  gas_conc_array = Vector([deepcopy(gas_concs) for i in 1:nblocks])

  # Experiment index for each column
  exp_num = freshape(spread(convert(Array,collect(1:nexp_l)'), 1, ncol_l), [blocksize, nblocks], order=[1,2])

  # Water vapor and ozone depend on col, lay, exp: look just like other fields
  water_vapor = Array{FT}(ds["water_vapor"][:])
  gas_conc_temp_3d = reshape(water_vapor, nlay_l, blocksize, nblocks ) .* read_scaling(ds, FT,"water_vapor")

  for b = 1:nblocks
    gas_conc_temp_3d_t = transpose(gas_conc_temp_3d[:,:,b])
    gas_conc_temp_3d_a = convert(Array, gas_conc_temp_3d_t)

    set_vmr!(gas_conc_array[b], "h2o", gas_conc_temp_3d_a)
  end

  ozone = Array{FT}(ds["ozone"][:])
  gas_conc_temp_3d = reshape(ozone, nlay_l, blocksize, nblocks ) * read_scaling(ds, FT,"ozone")

  for b = 1:nblocks
    set_vmr!(gas_conc_array[b],"o3", convert(Array, transpose(gas_conc_temp_3d[:,:,b])))
  end

  # All other gases are a function of experiment only
  for g = 1:length(gas_names)

    # Skip 3D fields above, also NO2 since RFMIP doesn't have this
    if gas_names[g] in ["h2o", "o3", "no2"]
      continue
    end

    # Read the values as a function of experiment
    gas_conc_temp_1d = ds[trim(names_in_file[g]) * "_GM"][:] * read_scaling(ds, FT,trim(names_in_file[g]) * "_GM")

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

  return gas_conc_array
end

read_scaling(ds, FT, varName) = parse(FT, ds[varName].attrib["units"])

end #module
