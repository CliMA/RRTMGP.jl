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
  using ..fortran_intrinsics
  # use mo_rte_kind,    only: FT
  # use mo_util_string, only: lower_case
  # implicit none
  const GAS_NOT_IN_LIST = -1

  mutable struct conc_field{FT}
    conc::Array{FT,2}
  end
  # type, private :: conc_field
  #   real(FT), dimension(:,:), allocatable :: conc
  # end type conc_field
  export set_vmr!, get_vmr

  export ty_gas_concs
  mutable struct ty_gas_concs{FT}
    gas_name::Vector{String}
    # concs::Vector{conc_field{FT}}
    concs
    ncol#::Int
    nlay#::Int
  end
  function ty_gas_concs(::Type{FT}, ncol, nlay) where FT
    gas_name = Vector([])
    # conc = Array{FT,2}(undef, ncol, nlay)
    # concs = Vector([conc_field(conc)])
    concs = Vector([])
    # concs = nothing
    ncol, nlay = nothing, nothing
    return ty_gas_concs{FT}(gas_name, concs, ncol, nlay)
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
  function set_vmr!(this::ty_gas_concs{FT}, gas, w, sw) where FT
    if length(sw)==1 && length(w)==1
      set_vmr!(this, gas, w[1])
    elseif length(sw)==1 && length(w)>1
      set_vmr!(this, gas, w)
    else
      set_vmr!(this, gas, w)
    end
  end

  function set_vmr!(this::ty_gas_concs{FT}, gas, w::FT) where FT # result(error_msg)
    # class(ty_gas_concs), intent(inout) :: this
    # character(len=*),    intent(in   ) :: gas
    # real(FT),            intent(in   ) :: w

    if w < FT(0) || w > FT(1)
      error("ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1")
    end

    igas = find_gas(this, gas)

    conc = Array{FT}(undef, 1,1)
    fill!(conc, w)
    gas_name = trim(gas)

    # if igas == GAS_NOT_IN_LIST || igas>length(this.concs)
    if igas == GAS_NOT_IN_LIST
      push!(this.concs, conc_field(conc))
      push!(this.gas_name, gas_name)
      igas = length(this.concs)
      # L = lowercase.(strip.(this.gas_name)) .== lowercase(trim(gas))
      # @show L
      # @show gas
      # @show igas
    end
    this.concs[igas].conc = conc
    this.gas_name[igas] = gas_name
    #
    # Deallocate anything existing -- could be more efficient to test if it's already the correct size
    #

  end

function set_vmr!(this::ty_gas_concs{FT}, gas, w::Vector{FT}) where FT
  # class(ty_gas_concs), intent(inout) :: this
  # character(len=*),    intent(in   ) :: gas
  # real(FT), dimension(:),
      #                      intent(in   ) :: w

  if any(w .< FT(0)) || any(w .> FT(1))
    error("ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1")
  end
  if this.nlay ≠ nothing && length(w) ≠ this.nlay
    error("ty_gas_concs%set_vmr: different dimension (nlay)")
  else
    this.nlay = length(w)
  end

  igas = find_gas(this, gas)
  conc = Array{FT}(undef, 1, this.nlay)
  conc .= reshape(w, 1, this.nlay)
  gas_name = trim(gas)

  if igas == GAS_NOT_IN_LIST || igas == length(this.concs)+1
    push!(this.concs, conc_field(conc))
    push!(this.gas_name, gas_name)
    igas = length(this.concs)
  end
  this.concs[igas].conc = conc
  this.gas_name[igas] = gas_name
end

function set_vmr!(this::ty_gas_concs, gas, w::Array{FT}) where FT
  # class(ty_gas_concs), intent(inout) :: this
  # character(len=*),    intent(in   ) :: gas
  # real(FT), dimension(:,:),
  #                      intent(in   ) :: w

  if any(w .< FT(0)) || any(w .> FT(1))
    error("set_vmr: concentrations should be >= 0, <= 1")
  end
  if this.ncol ≠ nothing && size(w, 1) ≠ this.ncol
    @show size(w)
    @show size(w, 1)
    @show this.ncol
    error("set_vmr: different dimension (ncol)")
  else
    this.ncol = size(w, 1)
  end
  if this.nlay ≠ nothing && size(w, 2) ≠ this.nlay
    error("set_vmr: different dimension (nlay)")
  else
    this.nlay = size(w, 2)
  end

  conc = w
  gas_name = trim(gas)
  igas = find_gas(this, gas)
  # @show "Array", gas, igas, igas == GAS_NOT_IN_LIST

  # if igas == GAS_NOT_IN_LIST || igas>length(this.concs)
  if igas == GAS_NOT_IN_LIST
    push!(this.concs, conc_field(conc))
    push!(this.gas_name, gas_name)
    igas = length(this.concs)
  end
  this.concs[igas].conc = conc
  this.gas_name[igas] = gas_name
end

  #
  # Return volume mixing ratio as 1D or 2D array
  #
  # -------------------------------------------------------------------------------------
  #
  # 1D array ( lay depdendence only)
  #
  function get_vmr(this::ty_gas_concs{FT}, gas) where FT #result(error_msg)
    # class(ty_gas_concs) :: this
    # character(len=*),         intent(in ) :: gas
    # real(FT), dimension(:),   intent(out) :: array
    # character(len=128) :: error_msg
    # # ---------------------
    # integer :: igas
    # # ---------------------

    igas = find_gas(this, gas)
    if igas == GAS_NOT_IN_LIST
      error("get_vmr: gas " * trim(gas) * " not found")
    # elseif size(this.concs[igas].conc, 1) > 1 # Are we requesting a single profile when many are present?
    #   error("get_vmr: gas " * trim(gas) * " requesting single profile but many are available")
    end
    conc = this.concs[igas].conc

    # if size(this.concs[igas].conc, 2) > 1
    #   array = this.concs[igas].conc[1,:]
    # else
    #   array = this.concs[igas].conc[1,1]
    # end

    if size(this.concs[igas].conc, 2) > 1
      array = this.concs[igas].conc[1,:]
    else
      array = this.concs[igas].conc[1,1]
    end


    # if size(conc, 1) > 1     # Concentration stored as 2D
    #   array .= conc[1,:]
    # elseif size(conc, 2) > 1 # Concentration stored as 1D
    #   array .= spread(conc[1,:], dim=1, ncopies=max(this.ncol, size(array,1)))
    #   # array .= conc[1,:]
    # else                     # Concentration stored as scalar
    #   fill!(array, conc[1,1])
    # end

    # if this.ncol ≠ nothing && this.ncol ≠ size(array,1)
    #   @show size(conc, 1) > 1
    #   @show size(conc, 2) > 1
    #   @show gas, igas
    #   @show size(conc)
    #   @show size(array)
    #   @show this.ncol
    #   @show this.nlay
    #   error("get_vmr: gas " * trim(gas) * " array is wrong size (ncol)")
    # end
    # if this.nlay ≠ nothing && this.nlay ≠ size(array,2)
    #   @show size(conc, 1) > 1
    #   @show size(conc, 2) > 1
    #   @show gas, igas
    #   @show size(conc)
    #   @show size(array)
    #   @show this.ncol
    #   @show this.nlay
    #   error("get_vmr: gas " * trim(gas) * " array is wrong size (nlay)")
    # end

    return array

  end
  # -------------------------------------------------------------------------------------
  #
  # 2D array (col, lay)
  #
  function get_vmr(this::ty_gas_concs, gas, array::Array{FT,2}) where FT
    # class(ty_gas_concs) :: this
    # character(len=*),         intent(in ) :: gas
    # real(FT), dimension(:,:), intent(out) :: array
    # character(len=128)                    :: error_msg
    # # ---------------------
    # integer :: igas
    # # ---------------------

    igas = find_gas(this, gas)
    if igas == GAS_NOT_IN_LIST
      error("get_vmr: gas " * trim(gas) * " not found")
    end
    #
    # Is the requested array the correct size?
    #
    if this.ncol ≠ nothing && this.ncol ≠ size(array,1)
      error("get_vmr: gas " * trim(gas) * " array is wrong size (ncol)")
    end
    if this.nlay ≠ nothing && this.nlay ≠ size(array,2)
      error("get_vmr: gas " * trim(gas) * " array is wrong size (nlay)")
    end

    if size(this.concs[igas].conc, 1) > 1      # Concentration stored as 2D
      array = this.concs[igas].conc[:,:]
    elseif size(this.concs(igas).conc, 2) > 1 # Concentration stored as 1D
      array = spread(this.concs(igas).conc(1,:), dim=1, ncopies=max(this.ncol, size(array,1)))
    else                                                   # Concentration stored as scalar
      array = this.concs[igas].conc[1,1]
    end
    return array

  end
  # -------------------------------------------------------------------------------------
  #
  # Extract a subset of n columns starting with column 'start'
  #
  # -------------------------------------------------------------------------------------
  function get_subset_range(this::ty_gas_concs, start, n, subset)
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
    if (this.ncol ≠ nothing && start > this.ncol || start+n-1 > this.ncol )
      error("gas_concs%get_vmr: Asking for columns outside range")
    end

    reset!(subset)
    # These two arrays should be the same length
    subset.gas_name = Array(undef, size(this.gas_name))
    subset.concs = Array(undef, size(this.concs))
    subset.nlay = this.nlay
    subset.ncol = fmerge(n, 0, this.ncol ≠ nothing)
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
  # find gas in list; GAS_NOT_IN_LIST if not found
  #
  function find_gas(this::ty_gas_concs, gas)
    L = lowercase.(strip.(this.gas_name)) .== lowercase(trim(gas))
    !any(L) && return GAS_NOT_IN_LIST
    return argmax(L)
  end

  # function find_gas(this::ty_gas_concs, gas)
  #   # character(len=*),   intent(in) :: gas
  #   # class(ty_gas_concs), intent(in) :: this
  #   # integer                        :: find_gas
  #   # # -----------------
  #   # integer :: igas
  #   # # -----------------
  #   gas_names = this.gas_name
  #   return
  #   if !allocated(this.gas_name)
  #     return GAS_NOT_IN_LIST
  #   end
  #   println("********************************")
  #   @show gas
  #   @show this.gas_name
  #   for igas = 1:length(this.gas_name)
  #     @show igas, this.gas_name[igas], gas
  #     if lowercase(trim(this.gas_name[igas])) == lowercase(trim(gas))
  #       return igas
  #     end
  #   end
  #   println("********************************")
  # end

end # module
