"""
    ty_optical_props_nstr{FT,I} <: ty_optical_props_arry{FT,I}

Holds extinction optical depth `tau`, `ssa`, and phase function moments `p` with leading dimension nmom.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ty_optical_props_nstr{FT,I} <: ty_optical_props_arry{FT,I}
  base
  tau#::Array{FT,3}
  ssa#::Array{FT,3}
  p#::Array{FT,4}
end


"""
  alloc!()

  class(ty_optical_props_nstr)    :: this
  integer,             intent(in) :: nmom # number of moments
  integer,             intent(in) :: ncol, nlay
  character(len=128)              :: err_message
"""
function alloc!(this::ty_optical_props_nstr, nmom, ncol, nlay)
  if any([ncol, nlay] .<= 0)
    error("optical_props%alloc: must provide positive extents for ncol, nlay")
  else
    allocated(this.tau) && deallocate!(this.tau)
    this.tau = Array(undef, ncol,nlay,get_ngpt(this))
  end
  allocated(this.ssa) && deallocate!(this.ssa)
  this.ssa = Array(undef, ncol,nlay,get_ngpt(this))
  allocated(this.p) && deallocate!(this.p)
  this.p = Array(undef, nmom,ncol,nlay,get_ngpt(this))
end


"""
    copy_and_alloc_nstr!(...)

    class(ty_optical_props_nstr)             :: this
    integer,                      intent(in) :: nmom, ncol, nlay
    class(ty_optical_props     ), intent(in) :: spectral_desc
    character(len=*), optional,   intent(in) :: name
    character(len=128)                       :: err_message
"""
function copy_and_alloc!(this::ty_optical_props_nstr, nmom, ncol, nlay, spectral_desc::ty_optical_props, name=nothing)
  is_initialized(this) && finalize!(this)
  init!(this, name, get_band_lims_wavenumber(spectral_desc),
                    get_band_lims_gpoint(spectral_desc))
  alloc!(this, nmom, ncol, nlay)
end

"""
    delta_scale!(...)

    class(ty_optical_props_nstr), intent(inout) :: this
    real(FT), dimension(:,:,:), optional,
                                 intent(in   ) :: for_
    character(128)                             :: err_message
"""
delta_scale!(this::ty_optical_props_nstr, for_) = "delta_scale_nstr: Not yet implemented"


"""

  class(ty_optical_props_nstr), intent(in) :: this
  character(len=128)                       :: err_message

  integer :: varSizes(3)
"""
function validate!(this::ty_optical_props_nstr{FT}) where FT

  #
  # Array allocation status, sizing
  #
  if !all([allocated(this.tau), allocated(this.ssa), allocated(this.p)])
    error("validate: arrays not allocated/initialized")
  end
  varSizes =   [size(this.tau, 1), size(this.tau, 2), size(this.tau, 3)]
  if      !all([size(this.ssa, 1), size(this.ssa, 2), size(this.ssa, 3)] == varSizes) ||
          !all([size(this.p,   2), size(this.p,   3), size(this.p,   4)] == varSizes)
    error("validate: arrays not sized consistently")
  end
  #
  # Valid values
  #
  any_vals_less_than(this.tau,  FT(0)) && error("validate: tau values out of range")
  any_vals_outside(this.ssa,  FT(0), FT(1)) && error("validate: ssa values out of range")
  any_vals_outside(this.p[1,:,:,:], FT(-1), FT(1)) && error("validate: p(1,:,:,:)  = g values out of range")
end

# ------------------------------------------------------------------------------------------
#
#  Routines for array classes: subsetting of optical properties arrays along x (col) direction
#
# Allocate class, then arrays; copy. Could probably be more efficient if
#   classes used pointers internally.
#
# This set takes start position and number as scalars
#
# ------------------------------------------------------------------------------------------

"""
  class(ty_optical_props_1scl), intent(inout) :: full
  integer,                      intent(in   ) :: start, n
  class(ty_optical_props_arry), intent(inout) :: subset
  character(128)                              :: err_message

  integer :: ncol, nlay, ngpt, nmom
"""
function subset_range!(full::ty_optical_props_1scl{FT}, start, n, subset) where FT
  if !is_initialized(full)
    error("optical_props%subset: Asking for a subset of uninitialized data")
  end
  ncol = get_ncol(full)
  nlay = get_nlay(full)
  ngpt = get_ngpt(full)
  if start < 1 || start + n-1 > get_ncol(full)
     error("optical_props%subset: Asking for columns outside range")
  end

  is_initialized(subset) && finalize!(subset)
  init!(subset, full)
  # Seems like the deallocation statements should be needed under Fortran 2003
  #   but Intel compiler doesn't run without them

  allocated(subset.tau) && deallocate!(subset.tau)
  if subset isa ty_optical_props_1scl
      alloc!(subset, n, nlay)
  elseif subset isa ty_optical_props_2str
      allocated(subset.ssa) && deallocate!(subset.ssa)
      allocated(subset.g  ) && deallocate!(subset.g  )
      alloc!(subset, n, nlay)
      subset.ssa[1:n,:,:] .= FT(0)
      subset.g[1:n,:,:] .= FT(0)
  elseif subset isa ty_optical_props_nstr
      allocated(subset.ssa) && deallocate!(subset.ssa)
      if allocated(subset.p)
        nmom = get_nmom(subset)
        allocated(subset.p  ) && deallocate!(subset.p  )
      else
        nmom = 1
      end
      alloc!(subset, nmom, n, nlay)
      subset.ssa[1:n,:,:] .= FT(0)
      subset.p[:,1:n,:,:] .= FT(0)
  else
    error("Uncaught case in subset_range!(full::ty_optical_props_1scl{FT}")
  end
  extract_subset!(ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)

end
# ------------------------------------------------------------------------------------------

"""
  subset_range!(...)

  class(ty_optical_props_2str), intent(inout) :: full
  integer,                      intent(in   ) :: start, n
  class(ty_optical_props_arry), intent(inout) :: subset
  character(128)                              :: err_message

  integer :: ncol, nlay, ngpt, nmom
"""
function subset_range!(full::ty_optical_props_2str{FT}, start, n, subset) where FT

  if !is_initialized(full)
    error("optical_props%subset: Asking for a subset of uninitialized data")
  end
  ncol = get_ncol(full)
  nlay = get_nlay(full)
  ngpt = get_ngpt(full)
  if start < 1 || start + n-1 > get_ncol(full)
     error("optical_props%subset: Asking for columns outside range")
  end

  is_initialized(subset) && finalize!(subset)
  init!(subset, full)

  if subset isa ty_optical_props_1scl # TODO: check logic
    alloc!(subset, n, nlay)
    extract_subset!(ncol, nlay, ngpt, full.tau, full.ssa, start, start+n-1, subset.tau)
  elseif subset isa ty_optical_props_2str
    allocated(subset.ssa) && deallocate!(subset.ssa)
    allocated(subset.g  ) && deallocate!(subset.g  )
    alloc!(subset, n, nlay)
    extract_subset!(ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)
    extract_subset!(ncol, nlay, ngpt, full.ssa, start, start+n-1, subset.ssa)
    extract_subset!(ncol, nlay, ngpt, full.g  , start, start+n-1, subset.g  )
  elseif subset isa ty_optical_props_nstr
    allocated(subset.ssa) && deallocate!(subset.ssa)
    if allocated(subset.p)
      nmom = get_nmom(subset)
      allocated(subset.p  ) && deallocate!(subset.p  )
    else
      nmom = 1
    end
    alloc!(subset, nmom, n, nlay)
    extract_subset!(ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)
    extract_subset!(ncol, nlay, ngpt, full.ssa, start, start+n-1, subset.ssa)
    subset.p[1,1:n,:,:] .= full.g[start:start+n-1,:,:]
    subset.p[2:end,:, :,:] .= FT(0) # TODO Verify this line
  else
    error("Uncaught case in subset_range!(full::ty_optical_props_2str{FT}")
  end
end

"""
    subset_range!(full::ty_optical_props_nstr, start, n, subset)

    class(ty_optical_props_nstr), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props_arry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom
"""
function subset_range!(full::ty_optical_props_nstr, start, n, subset)

  if !is_initialized(full)
    error("optical_props%subset: Asking for a subset of uninitialized data")
  end
  ncol = get_ncol(full)
  nlay = get_nlay(full)
  ngpt = get_ngpt(full)
  if(start < 1 || start + n-1 > get_ncol(full))
     error("optical_props%subset: Asking for columns outside range")
  end

  is_initialized(subset) && finalize!(subset)
  init!(subset, full)

  allocated(subset.tau) && deallocate!(subset.tau)
  if subset isa ty_optical_props_1scl # TODO: check logic
    alloc!(subset, n, nlay)
    extract_subset!(ncol, nlay, ngpt, full.tau, full.ssa, start, start+n-1, subset.tau)
  elseif subset isa ty_optical_props_2str
    allocated(subset.ssa) && deallocate!(subset.ssa)
    allocated(subset.g  ) && deallocate!(subset.g  )
    alloc!(subset, n, nlay)
    extract_subset!(ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)
    extract_subset!(ncol, nlay, ngpt, full.ssa, start, start+n-1, subset.ssa)
    subset.g[1:n,:,:] .= full.p[1,start:start+n-1,:,:]
  elseif subset isa ty_optical_props_nstr
    allocated(subset.ssa) && deallocate!(subset.ssa)
    allocated(subset.p  ) && deallocate!(subset.p  )
    alloc!(subset, nmom, n, nlay)
    extract_subset!(      ncol, nlay, ngpt, full.tau, start, start+n-1, subset.tau)
    extract_subset!(      ncol, nlay, ngpt, full.ssa, start, start+n-1, subset.ssa)
    extract_subset!(nmom, ncol, nlay, ngpt, full.p  , start, start+n-1, subset.p  )
  else
    error("Uncaught case in subset_range!(full::ty_optical_props_nstr, start, n, subset)")
  end
end


"""
    class(ty_optical_props_nstr), intent(in   ) :: this
    integer                                     :: get_nmom
"""
get_nmom(this::ty_optical_props_nstr) = allocated(this.p) ? size(this.p, 1) : 0


"""
    combine_and_reorder_nstr!(...)

Combine absoprtion and Rayleigh optical depths for total tau, ssa, p
using Rayleigh scattering phase function

integer, intent(in) :: ncol, nlay, ngpt, nmom
real(FT), dimension(ngpt,nlay,ncol), intent(in ) :: tau_abs, tau_rayleigh
real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: tau, ssa
real(FT), dimension(ncol,nlay,ngpt,nmom),
                                     intent(inout) :: p
# -----------------------
integer :: icol, ilay, igpt, imom
real(FT) :: t
# -----------------------
"""
function combine_and_reorder_nstr!(ncol, nlay, ngpt, nmom, tau_abs, tau_rayleigh, tau, ssa, p)
  FT = eltype(tau_abs)

  for icol in 1:ncol
    for ilay in 1:nlay
      for igpt in 1:ngpt
        t = tau_abs[igpt,ilay,icol] + tau_rayleigh[igpt,ilay,icol]
        tau[icol,ilay,igpt] = t
        if (t > FT(2) * realmin(FT))
          ssa[icol,ilay,igpt] = tau_rayleigh[igpt,ilay,icol] / t
        else
          ssa[icol,ilay,igpt] = FT(0)
        end
        for imom = 1:nmom
          p[imom,icol,ilay,igpt] = FT(0)
        end
        if (nmom >= 2)
          p[2,icol,ilay,igpt] = FT(0.1)
        end
      end
    end
  end
end
