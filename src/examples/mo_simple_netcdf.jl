module mo_simple_netcdf

using ..fortran_intrinsics
using NCDatasets

export read_field, read_string, write_field, var_exists, get_dim_size, get_array

#  public :: dim_exists, get_dim_size, create_dim, &
#            var_exists, get_var_size, create_var, &
#            read_field, read_string, read_char_vec, read_logical_vec, write_field

function read_string(ds, varName)

  if !haskey(ds,varName)
    return varName
  end

  field = ds[varName][:]
  return field
end

function read_field(ds, varName)

  if !haskey(ds,varName)
    error("read_field: can't find variable " * trim(varName))
  end

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

var_exists(ds, varName) = haskey(ds,varName) ? true : false

get_array(ds, name, FT) = haskey(ds, name) ? convert(Array{FT}, ds[name][:]) : nothing
get_array(ds, name, FT, s) = haskey(ds, name) ? convert(Array{FT}, ds[name][:]) : zeros(s)
get_dim_size(ds, name) = ds.dim[name]

#contains
#  !--------------------------------------------------------------------------------------------------------------------
#  function read_scalar(ncid, varName)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    real(FT)                     :: read_scalar

#    integer :: varid

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_field: can't find variable " // trim(varName))
#    if(nf90_get_var(ncid, varid, read_scalar)  /= NF90_NOERR) &
#      error("read_field: can't read variable " // trim(varName))

#  end function read_scalar
#  !--------------------------------------------------------------------------------------------------------------------
#  function read_1d_field(ncid, varName, nx)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    integer,          intent(in) :: nx
#    real(FT), dimension(nx)      :: read_1d_field

#    integer :: varid

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_field: can't find variable " // trim(varName))
#    if(nf90_get_var(ncid, varid, read_1d_field)  /= NF90_NOERR) &
#      error("read_field: can't read variable " // trim(varName))

#  end function read_1d_field
#  !--------------------------------------------------------------------------------------------------------------------
#  function read_2d_field(ncid, varName, nx, ny)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    integer,          intent(in) :: nx, ny
#    real(FT), dimension(nx, ny)  :: read_2d_field

#    integer :: varid

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_field: can't find variable " // trim(varName))
#    if(nf90_get_var(ncid, varid, read_2d_field)  /= NF90_NOERR) &
#      error("read_field: can't read variable " // trim(varName))

#  end function read_2d_field
#  !--------------------------------------------------------------------------------------------------------------------
#  function read_3d_field(ncid, varName, nx, ny, nz)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    integer,          intent(in) :: nx, ny, nz
#    real(FT), dimension(nx, ny, nz)  :: read_3d_field

#    integer :: varid

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_field: can't find variable " // trim(varName))
#    if(nf90_get_var(ncid, varid, read_3d_field)  /= NF90_NOERR) &
#      error("read_field: can't read variable " // trim(varName))

#  end function read_3d_field
#  !--------------------------------------------------------------------------------------------------------------------
#  function read_4d_field(ncid, varName, nw, nx, ny, nz)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    integer,          intent(in) :: nw, nx, ny, nz
#    real(FT), dimension(nw, nx, ny, nz)  :: read_4d_field

#    integer :: varid

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_field: can't find variable " // trim(varName))
#    if(nf90_get_var(ncid, varid, read_4d_field)  /= NF90_NOERR) &
#      error("read_field: can't read variable " // trim(varName))
#
#  end function read_4d_field
#  !--------------------------------------------------------------------------------------------------------------------
#  function read_string(ncid, varName, nc)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    integer,          intent(in) :: nc
#    character(len=nc)            :: read_string

#    integer :: varid

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) then
#      read_string = ""
#      return
#    end if
#    if(nf90_get_var(ncid, varid, read_string)  /= NF90_NOERR) &
#      error("read_field: can't read variable " // trim(varName))
#  end function read_string

#  function read_logical_vec(ncid, varName, nx)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    integer,          intent(in) :: nx
#    integer,      dimension(nx) :: read_logical_tmp
#    logical(wl),  dimension(nx) :: read_logical_vec

#    integer :: varid
#    integer :: ix

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_logical_vec: can't find variable " // trim(varName))
#    if(nf90_get_var(ncid, varid, read_logical_tmp)  /= NF90_NOERR) &
#      error("read_logical_vec: can't read variable " // trim(varName))
#    do ix = 1, nx
#      if (read_logical_tmp(ix) .eq. 0) then
#        read_logical_vec(ix) = .false.
#      else
#        read_logical_vec(ix) = .true.
#      endif
#    enddo

#  end function read_logical_vec
#  !--------------------------------------------------------------------------------------------------------------------
#  function read_char_vec(ncid, varName, nx)
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    integer,          intent(in) :: nx
#    character(len=32), dimension(nx) :: read_char_vec

#    integer :: varid

#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("read_char_vec: can't find variable " // trim(varName))
#    if(nf90_get_var(ncid, varid, read_char_vec)  /= NF90_NOERR) &
#      error("read_char_vec: can't read variable " // trim(varName))

#  end function read_char_vec
#  !--------------------------------------------------------------------------------------------------------------------
#  function dim_exists(ncid, dimName)
#    !
#    ! Does this dimension exist (have a valid dim_id) in the open netCDF file?
#    !
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: dimName
#    logical                      :: dim_exists

#    integer :: dimid
#    dim_exists = nf90_inq_dimid(ncid, trim(dimName), dimid) == NF90_NOERR
#  end function dim_exists
#  !--------------------------------------------------------------------------------------------------------------------

#  !--------------------------------------------------------------------------------------------------------------------
#  subroutine create_dim(ncid, dimName, dimLength)
#    !
#    ! Check to see if a dimiable with this name exists in the file
#    !   If so, check against current size
#    !   If not, create with specified dimensions
#    !
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: dimName
#    integer,          intent(in) :: dimLength

#    integer                 :: i, dimid

#    if(dim_exists(ncid, dimName)) then
#      if (dimLength /= get_dim_size(ncid, trim(dimName))) &
#          error("dim " // trim(dimName) // " is present but incorrectly sized.")
#    else
#      if(nf90_redef(ncid) /= NF90_NOERR) &
#        error("create_dim: can't put file into redefine mode")
#      if(nf90_def_dim(ncid, dimName, dimLength, dimid) /= NF90_NOERR) &
#        error("create_dim: can't define dimension " // trim(dimName))
#      if(nf90_enddef(ncid) /= NF90_NOERR) &
#        error("create_dim: can't end redefinition??")
#    end if
#  end subroutine create_dim
#  !--------------------------------------------------------------------------------------------------------------------
#  subroutine create_var(ncid, varName, dimNames, dimLengths, dataType)
#    !
#    ! Check to see if a variable with this name exists in the file
#    !   If so, check against current size
#    !   If not, create with specified dimensions
#    ! datatype: NF90_DOUBLE, NF90_FLOAT, NF90_INT, etc.
#    !
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    character(len=*), intent(in) :: dimNames(:)
#    integer,          intent(in) :: dimLengths(:)
#    integer, optional, intent(in) :: dataType

#    integer :: i, varid, xtype
#    integer :: dimIds(size(dimNames))

#    if(var_exists(ncid, varName)) then
#      do i = 1, size(dimNames)
#        if (dimLengths(i) /= get_dim_size(ncid, trim(dimNames(i)))) &
#          error("Variable " // trim(varName) // " is present but incorrectly sized.")
#      end do
#    else
#      do i = 1, size(dimNames)
#        if(nf90_inq_dimid(ncid, trim(dimNames(i)), dimIds(i)) /= NF90_NOERR) &
#          error("create_var: Can't get id for dimension " // trim(dimnames(i)))
#      end do
#      if(nf90_redef(ncid) /= NF90_NOERR) &
#        error("create_var: can't put file into redefine mode")
#      xtype = NF90_DOUBLE
#      if(present(dataType)) xtype = dataType
#      if(nf90_def_var(ncid, varName, xtype, dimIds, varid) /= NF90_NOERR) &
#        error("create_var: can't define variable " // trim(varName))
#      if(nf90_enddef(ncid) /= NF90_NOERR) &
#        error("create_dim: can't end redefinition??")
#    end if
#  end subroutine create_var

#  !--------------------------------------------------------------------------------------------------------------------
#  function get_var_size(ncid, varName, n)
#    !
#    ! Returns the extents of a netcdf variable on disk
#    !
#    integer,          intent(in) :: ncid
#    character(len=*), intent(in) :: varName
#    integer,          intent(in) :: n
#    integer                      :: get_var_size(n)

#    integer :: i
#    integer :: varid, ndims, dimids(n)

#    get_var_size(n) = -1
#    if(nf90_inq_varid(ncid, trim(varName), varid) /= NF90_NOERR) &
#      error("get_var_size: can't find variable " // trim(varName))
#    if(nf90_inquire_variable(ncid, varid, ndims = ndims) /= NF90_NOERR) &
#      error("get_var_size: can't get information for variable " // trim(varName))
#    if(ndims /= n) &
#      error("get_var_size:  variable " // trim(varName) // " has the wrong number of dimensions" )
#    if(nf90_inquire_variable(ncid, varid, dimids = dimids) /= NF90_NOERR) &
#      error("get_var_size: can't read dimension ids for variable " // trim(varName))
#    do i = 1, n
#      if(nf90_inquire_dimension(ncid, dimids(i), len = get_var_size(i)) /= NF90_NOERR) &
#        error("get_var_size: can't get dimension lengths for variable " // trim(varName))
#    end do

#  end function get_var_size
#  !--------------------------------------------------------------------------------------------------------------------
end
