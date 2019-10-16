# This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
#
# Contacts: Robert Pincus and Eli Mlawer
# email:  rrtmgp@aer.com
#
# Copyright 2015-2018,  Atmospheric and Environmental Research and
# Regents of the University of Colorado.  All right reserved.
#
# Use and duplication is permitted under the terms of the
#    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
# -------------------------------------------------------------------------------------------------
#
# This module reads an example file containing atomspheric conditions (temperature, pressure, gas concentrations)
#   and surface properties (emissivity, temperature), defined on nlay layers across a set of ncol columns subject to
#   nexp perturbations, and returns them in data structures suitable for use in rte and rrtmpg. The input data
#   are partitioned into a user-specified number of blocks.
# For the moment only quantities relevant to longwave calculations are provided.
#
# The example files comes from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
#   The protocol for this experiment allows for different specifications of which gases to consider:
# all gases, (CO2, CH4, N2O) + {CFC11eq; CFC12eq + HFC-134eq}. Ozone is always included
# The protocol does not specify the treatmet of gases like CO
#
# -------------------------------------------------------------------------------------------------

module mo_rfmip_io
#  use mo_rte_kind,   only: FT
#  use mo_gas_concentrations, &
#                        only: ty_gas_concs
#  use mo_util_string,   only: lower_case, string_in_array, string_loc_in_array
#  use mo_simple_netcdf, only: read_field, write_field, get_dim_size
#  use netcdf
#  implicit none
#  private

  using ..mo_gas_concentrations
  using ..mo_util_string
  using ..fortran_intrinsics
  using NCDatasets

  export read_kdist_gas_names, determine_gas_names, read_size, read_and_block_pt,
            read_and_block_sw_bc, read_and_block_lw_bc, read_and_block_gases_ty,
            unblock_and_write

  # ncol_l = 0
  # nlay_l = 0
  # nexp_l = 0 # Local copies
  #--------------------------------------------------------------------------------------------------------------------
  #
  # Find the size of the problem: columns, layers, perturbations (experiments)
  #
  function read_size(ds)
#    character(len=*),          intent(in   ) :: fileName
#    integer,         optional, intent(  out) :: ncol, nlay, nexp
#    ! ---------------------------
#    integer :: ncid
    # ---------------------------
#    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
#      error("read_size: can't find file " // trim(fileName))

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

#    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
#      error("read_and_block_pt: can't find file " // trim(fileName))
    #
    # Read p, T data; reshape to suit RRTMGP dimensions
    #

    temp3d = reshape( repeat(ds["pres_layer"][:],1,1,nexp_l), nlay_l, blocksize, nblocks )

    for b = 1:nblocks
      p_lay[:,:,b] = transpose(temp3d[:,:,b])
    end

    temp3d = reshape( ds["temp_layer"][:], nlay_l, blocksize, nblocks )

    for b = 1:nblocks
      t_lay[:,:,b] = transpose(temp3d[:,:,b])
    end

    temp3d = []

    temp3d = reshape( repeat(ds["pres_level"][:],1,1,nexp_l), nlay_l+1, blocksize, nblocks )

    for b = 1:nblocks
      p_lev[:,:,b] = transpose(temp3d[:,:,b])
    end

    temp3d = reshape( ds["temp_level"][:],nlay_l+1, blocksize, nblocks )

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

    temp2D = repeat(ds["surface_albedo"][:],1,nexp_l)

    surface_albedo = reshape(temp2D,blocksize,nblocks)

    temp2D = repeat(ds["total_solar_irradiance"][:],1,nexp_l)
    total_solar_irradiance = reshape(temp2D, blocksize, nblocks)

    temp2D = repeat(ds["solar_zenith_angle"][:],1,nexp_l)
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

    nblocks = (ncol_l*nexp_l)/blocksize

    #
    # Allocate on assignment
    #
    temp2D = repeat(ds["surface_emissivity"][:],1,nexp_l)
    surface_emissivity  = reshape(temp2D,blocksize,nblocks)

    temp2D = repeat(ds["surface_temperature"][:],1,nexp_l)
    surface_temperature = reshape(temp2D,blocksize,nblocks)

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
#    ! ----------------
#    integer :: num_gases, i
#    character(len=32), dimension(11) :: &
      # chem_name = ["co", "ch4", "o2", "n2o", "n2", "co2", "CCl4", "ch4", "CH3Br", "CH3Cl", "cfc22"]
      chem_name = ["co", "ch4", "o2", "n2o", "n2", "co2", "ccl4", "ch4", "ch3br", "ch3cl", "cfc22"]

      conc_name = ["carbon_monoxide", "methane", "oxygen",
          	   "nitrous_oxide", "nitrogen",
        	   "carbon_dioxide", "carbon_tetrachloride",
        	   "methane", "methyl_bromide",
        	   "methyl_chloride", "hcfc22"]
    # ----------------
    if forcing_index == 1
      names_in_kdist = read_kdist_gas_names(ds)
      names_in_file = Array{String}(undef,length(names_in_kdist))
#      allocate(names_in_file(size(names_in_kdist)))
      for i = 1:length(names_in_kdist)
        names_in_file[i] = trim(lowercase(names_in_kdist[i]))
        #
        # Use a mapping between chemical formula and name if it exists
        #
        if string_in_array(names_in_file[i], chem_name)
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
#--------------------------------------------------------------------------------------------------------------------
#!
#! Read the names of the gases known to the k-distribution
#!
#!
  function read_kdist_gas_names(ds)
#    character(len=*),          intent(in   ) :: fileName
#    character(len=32), dimension(:), allocatable, &
#                               intent(  out) :: kdist_gas_names
#    ! ---------------------------
#    integer :: ncid, varid
#    character(len=9), parameter :: varName = "gas_names"
#    ! ---------------------------
    varName = "gas_names"


#    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
#      error("read_kdist_gas_names: can't open file " // trim(fileName))

#    allocate(kdist_gas_names(get_dim_size(ncid, 'absorber')))

    ngases = ds.dim["absorber"]

    kdist_gas_names = Array{String}(undef, ngases)

    if !haskey(ds,varName)
      error("read_kdist_gas_names: can't find variable " * trim(varName))
    end

    temp_str = ds[varName][:]

    for i = 1:ngases
      kdist_gas_names[i] = trim(join(temp_str[:,i]))
    end
    return kdist_gas_names


#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_kdist_gas_names: can't find variable " // trim(varName))
#    if(nf90_get_var(ncid, varid, kdist_gas_names)  /= NF90_NOERR) &
#      error("read_kdist_gas_names: can't read variable " // trim(varName))

#    ncid = nf90_close(ncid)
  end
  #!--------------------------------------------------------------------------------------------------------------------
  #!
  #! Read and reshape gas concentrations. RRTMGP requires gas concentrations to be supplied via a class
  #!   (ty_gas_concs). Gas concentrations are set via a call to gas_concs%set_vmr(name, values)
  #!   where `name` is nominally the chemical formula for the gas in question and `values` may be
  #!   a scalar, a 1-d profile assumed to apply to all columns, or an array of dimension (ncol, nlay).
  #! This routine outputs a vector nblocks long of these types so each element of the array can be passed to
  #!   the rrtmgp gas optics calculation in turn.
  #!
  #! This routine exploits RFMIP conventions: only water vapor and ozone vary by column within
  #!   each experiment.
  #! Fields in the RFMIP file have a trailing _GM (global mean); some fields use a chemical formula and other
  #!   a descriptive name, so a map is provided between these.
  #!
  function read_and_block_gases_ty(ds, blocksize, gas_names, names_in_file)
#    character(len=*),           intent(in   ) :: fileName
#    integer,                    intent(in   ) :: blocksize
#    character(len=*),  dimension(:), &
#                                intent(in   ) :: gas_names ! Names used by the k-distribution/gas concentration type
#    character(len=*),  dimension(:), &
#                                intent(in   ) :: names_in_file ! Corresponding names in the RFMIP file
#    type(ty_gas_concs), dimension(:), allocatable, &
#                                intent(  out) :: gas_conc_array

#    ! ---------------------------
#    integer :: ncid
#    integer :: nblocks
#    integer :: b, g
#    integer,  dimension(:,:),   allocatable :: exp_num
#    real(FT), dimension(:),     allocatable :: gas_conc_temp_1d
#    real(FT), dimension(:,:,:), allocatable :: gas_conc_temp_3d
#    ! ---------------------------
    ncol_l, nlay_l, nexp_l = read_size(ds)

    if any([ncol_l, nlay_l, nexp_l] .== 0)
      error("read_and_block_lw_bc: Haven't read problem size yet.")
    end
    if (ncol_l*nexp_l)%blocksize ≠ 0
      error("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
    end
    nblocks = Int((ncol_l*nexp_l)/blocksize)
    gas_concs = ty_gas_concs(Float64, ncol_l, nlay_l)
    gas_conc_array = Vector([deepcopy(gas_concs) for i in 1:nblocks])
    # gas_conc_array = Vector()
#    allocate(gas_conc_array(nblocks))
    # Experiment index for each colum

    # ORIGINAL: exp_num = reshape(spread([(b, b = 1, nexp_l)], 1, ncopies = ncol_l), shape = [blocksize, nblocks], order=[1,2])

    exp_num = freshape(spread(convert(Array,collect(1:nexp_l)'), 1, ncol_l), [blocksize, nblocks], order=[1,2])
    # test_data(exp_num, "exp_num")


#    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
#      error("read_and_block_gases_ty: can't find file " // trim(fileName))

    #
    # Water vapor and ozone depend on col, lay, exp: look just like other fields
    #
    # gas_conc_temp_3d = reshape(read_field(ncid, "water_vapor", nlay_l, ncol_l, nexp_l), &
    #                           shape = [nlay_l, blocksize, nblocks]) * read_scaling(ncid, "water_vapor")

    gas_conc_temp_3d = reshape( ds["water_vapor"][:], nlay_l, blocksize, nblocks ) .* read_scaling(ds,"water_vapor")
    # test_data(gas_conc_temp_3d, "water_vapor_gas_conc_temp_3d")

    for b = 1:nblocks
      gas_conc_temp_3d_t = transpose(gas_conc_temp_3d[:,:,b])
      gas_conc_temp_3d_a = convert(Array, gas_conc_temp_3d_t)

      set_vmr!(gas_conc_array[b], "h2o", gas_conc_temp_3d_a, size(gas_conc_temp_3d_a))
    end


    gas_conc_temp_3d = reshape( ds["ozone"][:], nlay_l, blocksize, nblocks ) * read_scaling(ds,"ozone")
    # test_data(gas_conc_temp_3d, "ozone_gas_conc_temp_3d")

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
      if string_in_array(gas_names[g], ["h2o", "o3", "no2"])
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
#          error(gas_conc_array(b)%set_vmr(gas_names(g), &
#          spread(gas_conc_temp_1d(exp_num(:,b)), 2, ncopies = nlay_l)))
          # gas_conc_temp_1d_rep = repeat(gas_conc_temp_1d, 1, nlay_l)
          # set_vmr!(gas_conc_array[b],gas_names[g], gas_conc_temp_1d_rep )

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
  #--------------------------------------------------------------------------------------------------------------------
#  function read_scaling(ncid, varName)
  function read_scaling(ds, varName) # modified to pass DS (dataset instead of ncid)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    real(FT)                     :: read_scaling

#    integer           :: varid
#    character(len=16) :: charUnits

    FT = Float64

    if haskey(ds,varName)
      scaling_str = ds[varName].attrib["units"]
    else
      error("read_scaling: can't find variable " * trim(varName))
    end

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_scaling: can't find variable " // trim(varName))
#    if(nf90_get_att(ncid, varid, "units", charUnits)  /= NF90_NOERR) &
#      error("read_scaling: can't read attribute 'units' from variable " // trim(varName))
#    read(charUnits, *) read_scaling
    return parse(FT, scaling_str)

  end
#--------------------------------------------------------------------------------------------------------------------
  #
  # Reshape and reorder values (nominally fluxes) from RTE order (ncol, nlev, nblocks)
  #   to RFMIP order (nlev, ncol, nexp), then write them to a user-specified variable
  #   in a netCDF file.
  #
  function unblock_and_write(ds, varName, values)
#    character(len=*),           intent(in   ) :: fileName, varName
#    real(FT), dimension(:,:,:),  & ! [blocksize, nlay/+1, nblocks]
#                                intent(in   ) :: values
#    ! ---------------------------
#    integer :: ncid
#    integer :: b, blocksize, nlev, nblocks
#    real(FT), dimension(:,:), allocatable :: temp2d
#    ! ---------------------------
    FT = Float64
    ncol_l, nlay_l, nexp_l = read_size(ds)

    if any([ncol_l, nlay_l, nexp_l] .== 0)
      error("unblock_and_write: Haven't read problem size yet.")
    end
    blocksize = size(values,1)
    nlev      = size(values,2)
    nblocks   = size(values,3)
    if nlev ≠ nlay_l+1
      error("unblock_and_write: array values has the wrong number of levels")
    end
    if (blocksize*nblocks) ≠ (ncol_l*nexp_l)
      error("unblock_and_write: array values has the wrong number of blocks/size")
    end

    temp2D = Array{FT}(undef,nlev,ncol_l*nexp_l)

    for b = 1:nblocks
      temp2D[1:nlev, ((b-1)*blocksize+1):(b*blocksize)] = transpose(values[1:blocksize,1:nlev,b])
    end
    #
    # Check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
    #
# SK - commenting - assuming it is appending an existing file
#    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
#      error("unblock_and_write: can't find file " // trim(fileName))
#    error(write_field(ncid, varName,  &
#                                 reshape(temp2d, shape = [nlev, ncol_l, nexp_l])))
#
#    ncid = nf90_close(ncid)
    # this portion needs to be double checked

    defVar(ds,varName, reshape(temp2d,nlev,ncol_l,nexp_l), ("nlev","ncol_l","nexp_l"))

  end
#-------------------------------------------
end
