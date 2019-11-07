"""
    mo_garand_atmos_io

NetCDF I/O routines, shared with other RTE+RRTMGP examples
"""
module mo_garand_atmos_io

using ..mo_gas_concentrations
using ..mo_util_reorder
using ..mo_optical_props

export read_atmos!

"""
    read_atmos!()

Read profiles for all columns  -- T, p, and gas concentrations
  Allocation occurs on assignments (says the F2003 standard)
"""
function read_atmos!(ds, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry)
  # character(len=*),   intent(in   ) :: fileName
  # real(wp), dimension(:,:), allocatable,                 &
  #                     intent(inout) :: p_lay, t_lay, p_lev, t_lev, col_dry
  # type(ty_gas_concs), intent(inout) :: gas_concs
  # # -------------------
  # integer :: igas
  # integer :: ds, ncol, nlay, nlev
  # integer, parameter :: ngas = 8
  # character(len=3), dimension(ngas) &
  #                    :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'co ', 'ch4', 'o2 ', 'n2 ']
  # character(len=7) :: vmr_name
  # # -------------------
  ngas = 8
  gas_names = ["h2o", "co2", "o3 ", "n2o", "co ", "ch4", "o2 ", "n2 "]

  ncol = get_dim_size(ds, "col")
  nlay = get_dim_size(ds, "lay")
  nlev = get_dim_size(ds, "lev")
  nlev ≠ nlay+1 && error("read_atmos: nlev should be nlay+1")

  #
  # These lines assume that compilers follow the Fortran 2003 standard for
  #   allocating on assignment. This may require explicit compiler support
  #   e.g. -assume realloc_lhs flag for Intel
  #
  p_lay = read_field(ds, "p_lay", ncol, nlay)
  t_lay = read_field(ds, "t_lay", ncol, nlay)
  p_lev = read_field(ds, "p_lev", ncol, nlev)
  t_lev = read_field(ds, "t_lev", ncol, nlev)

  init!(gas_concs, gas_names)
  for igas = 1:ngas
    vmr_name = "vmr_" * strip(gas_names(igas))
    !var_exists(ds, strip(vmr_name)) && error("read_atmos: can't read concentration of "*strip(gas_names(igas)))
    set_vmr!(gas_concs, strip(gas_names(igas)), read_field(ds, strip(vmr_name), ncol, nlay))
  end

  # col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
  var_exists(ds, "col_dry") && (col_dry = read_field(ds, "col_dry", ncol, nlay))
end

end # module