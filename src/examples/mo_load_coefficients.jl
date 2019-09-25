module mo_load_coefficients

  using ..mo_gas_concentrations
  using ..mo_gas_optics_rrtmgp
  using ..fortran_intrinsics
  using NCDatasets 

  export load_and_init, stop_on_err

#-------------------------
  function stop_on_err(msg)
#    use iso_fortran_env, only : error_unit
#    character(len=*), intent(in) :: msg

    if msg != ""
      error(msg)
    end
  end


#--------------------------
  function load_and_init(kdist::ty_gas_optics_rrtmgp, filename::String, available_gases::ty_gas_concs )

  # Reading the properties from the NetCDF file

  ds = Dataset(filename,"r")

  gas_names   			  =      ds["gas_names"][:]
  key_species 			  =      ds["key_species"][:]
  band_lims   			  =      ds["bnd_limits_wavenumber"][:]
  band2gpt          		  =  Int(ds["bnd_limits_gpt")][:]
  press_ref         		  = 	 ds["press_ref"][:]
  temp_ref          		  = 	 ds["temp_ref"][:]
  temp_ref_p        		  = 	 ds["absorption_coefficient_ref_P"][:]
  temp_ref_t        		  = 	 ds["absorption_coefficient_ref_T"][:]
  press_ref_trop    		  =      ds["press_ref_trop"][:]
  kminor_lower      		  =      ds["kminor_lower"][:]
  kminor_upper       		  =      ds["kminor_upper"][:]
  gas_minor 	    		  = 	 ds["gas_minor"][:]
  identifier_minor  		  = 	 ds["identifier_minor"][:]
  minor_gases_lower 		  = 	 ds["minor_gases_lower"][:]
  minor_gases_upper 		  = 	 ds["minor_gases_upper"][:]
  minor_limits_gpt_lower 	  =  Int(ds["minor_limits_gpt_lower")][:]
  minor_limits_gpt_upper  	  =  Int(ds["minor_limits_gpt_upper")][:]
  minor_scales_with_density_lower = Bool(ds["minor_scales_with_density_lower")][:]
  minor_scales_with_density_upper = Bool(ds["minor_scales_with_density_upper")][:]
  scale_by_complement_lower       = Bool(ds["scale_by_complement_lower")][:]
  scale_by_complement_upper       = Bool(ds["scale_by_complement_upper")][:]
  scaling_gas_lower               =      ds["scaling_gas_lower"][:]
  scaling_gas_upper               =      ds["scaling_gas_upper"][:]
  kminor_start_lower              =      ds["kminor_start_lower"][:]
  kminor_start_upper              =      ds["kminor_start_upper"][:]
  vmr_ref                         =      ds["vmr_ref"][:]
  kmajor            	          =      ds["kmajor"][:]

  if haskey(ds,"rayl_lower")
    rayl_lower         	          =      ds["rayl_lower"][:]
    rayl_upper         	          =      ds["rayl_upper"][:]
  end

 #   ! Initialize the gas optics class with data. The calls look slightly different depending
 #   !   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
 #   ! gas_optics%load() returns a string; a non-empty string indicates an error.

  if haskey(ds,"totplnk")
    totplnk         	          =      ds["totplnk"][:]
    planck_frac        	          =      ds["plank_fraction"][:]

    stop_on_err(load!(kdist, 
	  available_gases, 
	  gas_names, 
	  key_species, 
	  band2gpt, 
	  band_lims,   
          press_ref, 
	  press_ref_trop, 
	  temp_ref, 
	  temp_ref_p, 
	  temp_ref_t,     
          vmr_ref, 
	  kmajor, 
	  kminor_lower, 
	  kminor_upper, 
	  gas_minor,
	  identifier_minor,
          minor_gases_lower, 
	  minor_gases_upper, 
	  minor_limits_gpt_lower,
          minor_limits_gpt_upper, 
	  minor_scales_with_density_lower, 
          minor_scales_with_density_upper, 
          scaling_gas_lower, 
	  scaling_gas_upper, 
          scale_by_complement_lower, 
          scale_by_complement_upper, 
          kminor_start_lower, 
          kminor_start_upper, 
          totplnk, planck_frac,       
          rayl_lower, 
	  rayl_upper))

  else
	
    solar_src                     =      ds["solar_source"][:]

    stop_on_err(load!(kdist, 
	  available_gases, 
          gas_names,   
          key_species, 
          band2gpt,    
          band_lims,  
          press_ref, 
          press_ref_trop, 
          temp_ref,    
          temp_ref_p, 
	  temp_ref_t,     
          vmr_ref, 
	  kmajor,            
          kminor_lower, 
	  kminor_upper, 
          gas_minor,
	  identifier_minor,
          minor_gases_lower, 
	  minor_gases_upper, 
          minor_limits_gpt_lower, 
          minor_limits_gpt_upper, 
          minor_scales_with_density_lower, 
          minor_scales_with_density_upper, 
          scaling_gas_lower, 
	  scaling_gas_upper, 
          scale_by_complement_lower, 
          scale_by_complement_upper, 
          kminor_start_lower, 
          kminor_start_upper, 
          solar_src, 
          rayl_lower, 
	  rayl_upper))


  end

  #-------------------------------

  end
#--------------------------


end
