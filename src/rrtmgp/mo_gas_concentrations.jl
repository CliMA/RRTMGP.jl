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
# Encapsulates a collection of volume mixing ratios (concentrations) of gases.
#   Each concentration is associated with a name, normally the chemical formula.
#
# Values may be provided as scalars, 1-dimensional profiles (nlay), or 2-D fields (ncol,nlay).
#   (nlay and ncol are determined from the input arrays; self-consistency is enforced)
#   example:
#   error_msg = gas_concs%set_vmr('h2o', values(:,:))
#   error_msg = gas_concs%set_vmr('o3' , values(:)  )
#   error_msg = gas_concs%set_vmr('co2', value      )
#
# Values can be requested as profiles (valid only if there are no 2D fields present in the object)
#   or as 2D fields. Values for all columns are returned although the entire collection
#   can be subsetted in the column dimension
#
# Subsets can be extracted in the column dimension
#
# Functions return strings. Non-empty strings indicate an error.
#
# -------------------------------------------------------------------------------------------------

module mo_gas_concentrations
  # use mo_rte_kind,    only: wp
  # use mo_util_string, only: lower_case
  # implicit none
  const GAS_NOT_IN_LIST = -1

  struct conc_field{T}
    conc::Array{T,2}
  end
  # type, private :: conc_field
  #   real(wp), dimension(:,:), allocatable :: conc
  # end type conc_field

  export ty_gas_concs
  struct ty_gas_concs{T, I}
    gas_name::Vector{String}
    concs::Vector{conc_field{T}}
    ncol::I
    nlay::I
  end
#   type, public :: ty_gas_concs
#     #
#     # Data
#     #
#     character(len=32), dimension(:), allocatable :: gas_name
#     type(conc_field),  dimension(:), allocatable :: concs
#     integer :: ncol = 0, nlay = 0
#     contains
#       #
#       # Procedures
#       #
#       procedure, private :: increase_list_size
#       procedure, private :: find_gas
#       procedure, private :: set_vmr_scalar
#       procedure, private :: set_vmr_1d
#       procedure, private :: set_vmr_2d
#       procedure, private :: get_vmr_1d
#       procedure, private :: get_vmr_2d
#       procedure, private :: get_subset_range
#       #
#       # public interface
#       #
#       procedure, public :: reset
#       generic,   public :: set_vmr => set_vmr_scalar,
#                                       set_vmr_1d,
#                                       set_vmr_2d
#       generic,   public :: get_vmr => get_vmr_1d,
#                                       get_vmr_2d
#       generic,   public :: get_subset => get_subset_range
#       procedure, public :: get_num_gases
#       procedure, public :: get_gas_names
#   end type ty_gas_concs
# contains
  # -------------------------------------------------------------------------------------
  #
  # Set concentrations --- scalar, 1D, 2D
  #
  # -------------------------------------------------------------------------------------
  function set_vmr!(this::ty_gas_concs{DT}, gas, w::DT) where DT # result(error_msg)
    # class(ty_gas_concs), intent(inout) :: this
    # character(len=*),    intent(in   ) :: gas
    # real(wp),            intent(in   ) :: w
    # character(len=128)                 :: error_msg
    # # ---------
    # integer :: igas
    # # ---------
    if w < DT(0) || w > DT(1)
      error("ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1")
    end

    igas = find_gas(this, gas)
    if igas == GAS_NOT_IN_LIST
      increase_list_size!(this)
      igas = size(this.gas_name)
    end
    #
    # Deallocate anything existing -- could be more efficient to test if it's already the correct size
    #

    this.concs[igas].conc = Array(undef, 1,1)
    this.concs[igas].conc[:,:] .= w
    this.gas_name[igas] = trim(gas)
  end
  # -------------------------------------------------------------------------------------
  function set_vmr!(this::ty_gas_concs{DT}, gas, w::Vector{DT}) where DT
    # class(ty_gas_concs), intent(inout) :: this
    # character(len=*),    intent(in   ) :: gas
    # real(wp), dimension(:),
        #                      intent(in   ) :: w
    # character(len=128)                 :: error_msg
    # # ---------
    # integer :: igas
    # # ---------

    if any(w < DT(0) || w > DT(1))
      error("ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1")
    end
    if this.nlay > 0
      if (size(w) ≠ this.nlay)
        error("ty_gas_concs%set_vmr: different dimension (nlay)")
      end
    else
      this.nlay = size(w)
    end

    igas = this.find_gas(gas)
    if igas == GAS_NOT_IN_LIST
      increase_list_size!(this)
      igas = size(this.gas_name)
    end
    #
    # Deallocate anything existing -- could be more efficient to test if it's already the correct size
    #

    this.concs[igas].conc = Array(undef,1,this.nlay)
    this.concs[igas].conc[1,:] .= w
    this.gas_name[igas] = trim(gas)
  end
  # --------------------
  function set_vmr!(this::ty_gas_concs, gas, w::Array{DT}) where DT
    # class(ty_gas_concs), intent(inout) :: this
    # character(len=*),    intent(in   ) :: gas
    # real(wp), dimension(:,:),
        #                      intent(in   ) :: w
    # character(len=128)                 :: error_msg
    # # ---------
    # integer :: igas
    # # ---------

    if any(w < DT(0) || w > DT(1))
      error("ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1")
    end
    if this.ncol > 0 && size(w, 1) ≠ this.ncol
      error("ty_gas_concs%set_vmr: different dimension (ncol)")
    else
      this.ncol = size(w, 1)
    end
    if this.nlay > 0 && size(w, 2) ≠ this.nlay
      error("ty_gas_concs%set_vmr: different dimension (nlay)")
    else
      this.nlay = size(w, 2)
    end

    igas = find_gas(this, gas)
    if igas == GAS_NOT_IN_LIST
      increase_list_size!(this)
      igas = size(this.gas_name)
    end
    #
    # Deallocate anything existing -- could be more efficient to test if it's already the correct size
    #
    this.concs[igas].conc = Array(undef, this.ncol,this.nlay)
    this.concs[igas].conc[:,:] .= w
    this.gas_name[igas] = trim(gas)
  end
  # -------------------------------------------------------------------------------------
  #
  # Return volume mixing ratio as 1D or 2D array
  #
  # -------------------------------------------------------------------------------------
  #
  # 1D array ( lay depdendence only)
  #
  function get_vmr_1d(this::ty_gas_concs{DT}, gas, array) where DT #result(error_msg)
    # class(ty_gas_concs) :: this
    # character(len=*),         intent(in ) :: gas
    # real(wp), dimension(:),   intent(out) :: array
    # character(len=128) :: error_msg
    # # ---------------------
    # integer :: igas
    # # ---------------------

    igas = find_gas(this, gas)
    if igas == GAS_NOT_IN_LIST
      error("ty_gas_concs%get_vmr; gas " * trim(gas) * " not found")
      array[:] .= DT(0)
    elseif size(this.concs(igas).conc, 1) > 1 # Are we requesting a single profile when many are present?
      error("ty_gas_concs%get_vmr; gas " * trim(gas) * " requesting single profile but many are available")
      array[:] .= DT(0)
    end
    if this.nlay > 0 && this.nlay ≠ size(array)
      error("ty_gas_concs%get_vmr; gas " * trim(gas) * " array is wrong size (nlay)")
      array[:] .= DT(0)
    end

    if size(this.concs[igas].conc, 2) > 1
      array[:] .= this.concs[igas].conc[1,:]
    else
      array[:] .= this.concs[igas].conc[1,1]
    end

  end
  # -------------------------------------------------------------------------------------
  #
  # 2D array (col, lay)
  #
  function get_vmr_2d(this, gas, array)
    # class(ty_gas_concs) :: this
    # character(len=*),         intent(in ) :: gas
    # real(wp), dimension(:,:), intent(out) :: array
    # character(len=128)                    :: error_msg
    # # ---------------------
    # integer :: igas
    # # ---------------------

    igas = find_gas(this, gas)
    if igas == GAS_NOT_IN_LIST
      error("ty_gas_concs%get_vmr; gas " * trim(gas) * " not found")
      array[:,:] .= DT(0)
    end
    #
    # Is the requested array the correct size?
    #
    if this.ncol > 0 && this.ncol ≠ size(array,1)
      error("ty_gas_concs%get_vmr; gas " * trim(gas) * " array is wrong size (ncol)")
      array[:,:] .= DT(0)
    end
    if this.nlay > 0 && this.nlay ≠ size(array,2)
      error("ty_gas_concs%get_vmr; gas " * trim(gas) * " array is wrong size (nlay)")
      array[:,:] .= DT(0)
    end

    if size(this.concs[igas].conc, 1) > 1      # Concentration stored as 2D
      array[:,:] .= this.concs[igas].conc[:,:]
    elseif size(this.concs(igas).conc, 2) > 1 # Concentration stored as 1D
      array[:,:] .= spread(this.concs(igas).conc(1,:), dim=1, ncopies=max(this.ncol, size(array,1)))
    else                                                   # Concentration stored as scalar
      array[:,:] .= this.concs[igas].conc[1,1]
    end

  end
  # -------------------------------------------------------------------------------------
  #
  # Extract a subset of n columns starting with column 'start'
  #
  # -------------------------------------------------------------------------------------
  function get_subset_range(this, start, n, subset)
    # class(ty_gas_concs),      intent(in   ) :: this
    # integer,                  intent(in   ) :: start, n
    # class(ty_gas_concs),      intent(inout) :: subset
    # character(len=128)                      :: error_msg
    # # ---------------------
    # integer :: i
    # # ---------------------
    if (n <= 0)
      error("gas_concs%get_vmr: Asking for 0 or fewer columns ")
    end
    if (start < 1 )
      error("gas_concs%get_vmr: Asking for columns outside range")
    end
    if (this.ncol > 0 && start > this.ncol || start+n-1 > this.ncol )
      error("gas_concs%get_vmr: Asking for columns outside range")
    end

    reset!(subset)
    # These two arrays should be the same length
    subset.gas_name = Array(undef, size(this.gas_name))
    subset.concs = Array(undef, size(this.concs))
    subset.nlay = this.nlay
    subset.ncol = fmerge(n, 0, this.ncol > 0)
    subset.gas_name[:] .= this.gas_name[:]

    for i = 1:size(this.gas_name)
      #
      # Preserve scalar/1D/2D representation in subset,
      #   but need to ensure at least extent 1 in col dimension (ncol = 0 means no gas exploits this dimension)
      #
      allocate(subset.concs[i].conc(min(max(subset.ncol,1), size(this.concs[i].conc, 1)),
                                          min(    subset.nlay,    size(this.concs[i].conc, 2))))
      if size(this.concs[i].conc, 1) > 1      # Concentration stored as 2D
        subset.concs[i].conc[:,:] .= this.concs[i].conc[start:(start+n-1),:]
      else
        subset.concs[i].conc[:,:] .= this.concs[i].conc[:,:]
      end
    end

  end
  # -------------------------------------------------------------------------------------
  #
  # Deallocate memory
  #
  # -------------------------------------------------------------------------------------
  function reset!(this)
    # class(ty_gas_concs), intent(inout) :: this
    # # -----------------
    # integer :: i
    # # -----------------
    this.nlay = 0
    this.ncol = 0
    if allocated(this%concs)
      for i = 1:size(this.concs)
      end
    end
  end
  # -------------------------------------------------------------------------------------
  #
  # Inquiry functions
  #
  # -------------------------------------------------------------------------------------
  function get_num_gases(this)
    # class(ty_gas_concs), intent(in) :: this
    # integer :: get_num_gases

    return size(this.gas_name)
  end
  # -------------------------------------------------------------------------------------
  function get_gas_names(this)
    # class(ty_gas_concs), intent(in) :: this
    # character(len=32), dimension(this%get_num_gases()) :: get_gas_names

    return this.gas_name[:]
  end
  # -------------------------------------------------------------------------------------
  #
  # Private procedures
  #
  # -------------------------------------------------------------------------------------
  #
  # This routine is called when adding a new concentration if the
  #   the gas isn't in the list already
  #
  function increase_list_size(this::ty_gas_concs{DT}) where DT
    # class(ty_gas_concs), intent(inout) :: this
    # # -----------------
    # character(len=32), dimension(:), allocatable :: new_names
    # type(conc_field),  dimension(:), allocatable :: new_concs
    # # -----------------

    if allocated(this.gas_name)
      new_names = Array{DT}(undef, size(this.gas_name)+1)
      new_concs = Array{DT}(undef, size(this.gas_name)+1)
      new_names[1:size(this.gas_name)] = this.gas_name[:]
      new_concs[1:size(this.gas_name)] = this.concs[:]
      move_alloc!(new_names, this.gas_name)
      move_alloc!(new_concs, this.concs)
    else
      this.gas_name = Array{DT}(undef, 1)
      this.concs = Array{DT}(undef, 1)
    end
  end
  # -------------------------------------------------------------------------------------
  #
  # find gas in list; GAS_NOT_IN_LIST if not found
  #
  function find_gas(this::ty_gas_concs, gas)
    # character(len=*),   intent(in) :: gas
    # class(ty_gas_concs), intent(in) :: this
    # integer                        :: find_gas
    # # -----------------
    # integer :: igas
    # # -----------------
    find_gas = GAS_NOT_IN_LIST
    if !allocated(this.gas_name)
      return
    end
    for igas = 1:size(this.gas_name)
      if lowercase(trim(this.gas_name[igas])) == lowercase(trim(gas))
        find_gas = igas
      end
    end
  end

end # module
