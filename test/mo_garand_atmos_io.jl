#### NetCDF I/O routines, shared with other RTE+RRTMGP examples

"""
    read_atmos!()

Read profiles for all columns
 - `p_lay` pressure (layers)
 - `t_lay` temperature (layers)
 - `p_lev` pressure (levels)
 - `t_lev` temperature (levels)
 - `col_dry` gas concentrations
"""
function read_atmos!(ds, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry)
  ngas = 8
  gas_names = ["h2o", "co2", "o3 ", "n2o", "co ", "ch4", "o2 ", "n2 "]

  ncol = ds.dim["col"]
  nlay = ds.dim["lay"]
  nlev = ds.dim["lev"]
  nlev ≠ nlay+1 && error("read_atmos: nlev should be nlay+1")

  p_lay = read_field(ds, "p_lay", ncol, nlay)
  t_lay = read_field(ds, "t_lay", ncol, nlay)
  p_lev = read_field(ds, "p_lev", ncol, nlev)
  t_lev = read_field(ds, "t_lev", ncol, nlev)

  init!(gas_concs, gas_names)
  for igas = 1:ngas
    vmr_name = "vmr_" * strip(gas_names(igas))
    !haskey(ds, strip(vmr_name)) && error("read_atmos: can't read concentration of "*strip(gas_names(igas)))
    set_vmr!(gas_concs, strip(gas_names(igas)), read_field(ds, strip(vmr_name), ncol, nlay))
  end

  # col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
  haskey(ds, "col_dry") && (col_dry = read_field(ds, "col_dry", ncol, nlay))
end
