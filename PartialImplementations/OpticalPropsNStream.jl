"""
    NStream{FT,I} <: AbstractOpticalPropsArry{FT,I}

Holds extinction optical depth `tau`, `ssa`, and phase function moments `p` with leading dimension nmom.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct NStream{FT, I} <: AbstractOpticalPropsArry{FT, I}
    base::Any
    tau::Any#::Array{FT,3}
    ssa::Any#::Array{FT,3}
    p::Any#::Array{FT,4}
end

"""
    delta_scale!(...)

    class(NStream), intent(inout) :: this
    real(FT), dimension(:,:,:), optional,
                                 intent(in   ) :: for_
    character(128)                             :: err_message
"""
delta_scale!(this::NStream, for_) = "delta_scale_nstr: Not yet implemented"


"""

  class(NStream), intent(in) :: this
  character(len=128)                       :: err_message

  integer :: varSizes(3)
"""
function validate!(this::NStream{FT}) where {FT}

    #
    # Array allocation status, sizing
    #
    varSizes = [size(this.tau, 1), size(this.tau, 2), size(this.tau, 3)]
    if !all([size(this.ssa, 1), size(this.ssa, 2), size(this.ssa, 3)] == varSizes) ||
       !all([size(this.p, 2), size(this.p, 3), size(this.p, 4)] == varSizes)
        error("validate: arrays not sized consistently")
    end
    #
    # Valid values
    #
    any_vals_less_than(this.tau, FT(0)) && error("validate: tau values out of range")
    any_vals_outside(this.ssa, FT(0), FT(1)) && error("validate: ssa values out of range")
    any_vals_outside(this.p[1, :, :, :], FT(-1), FT(1)) && error("validate: p(1,:,:,:)  = g values out of range")
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
  class(OneScalar), intent(inout) :: full
  integer,                      intent(in   ) :: start, n
  class(AbstractOpticalPropsArry), intent(inout) :: subset
  character(128)                              :: err_message

  integer :: ncol, nlay, ngpt, nmom
"""
function subset_range!(full::OneScalar{FT}, start, n, subset) where {FT}
    ncol = get_ncol(full)
    nlay = get_nlay(full)
    ngpt = get_ngpt(full)
    if start < 1 || start + n - 1 > get_ncol(full)
        error("optical_props%subset: Asking for columns outside range")
    end

    init!(subset, full)
    # Seems like the deallocation statements should be needed under Fortran 2003
    #   but Intel compiler doesn't run without them

    if subset isa OneScalar
        alloc!(subset, n, nlay)
    elseif subset isa TwoStream
        alloc!(subset, n, nlay)
        subset.ssa[1:n, :, :] .= FT(0)
        subset.g[1:n, :, :] .= FT(0)
    elseif subset isa NStream
        if allocated(subset.p)
            nmom = get_nmom(subset)
        else
            nmom = 1
        end
        alloc!(subset, nmom, n, nlay)
        subset.ssa[1:n, :, :] .= FT(0)
        subset.p[:, 1:n, :, :] .= FT(0)
    else
        error("Uncaught case in subset_range!(full::OneScalar{FT}")
    end
    extract_subset!(ncol, nlay, ngpt, full.tau, start, start + n - 1, subset.tau)

end
# ------------------------------------------------------------------------------------------

"""
  subset_range!(...)

  class(TwoStream), intent(inout) :: full
  integer,                      intent(in   ) :: start, n
  class(AbstractOpticalPropsArry), intent(inout) :: subset
  character(128)                              :: err_message

  integer :: ncol, nlay, ngpt, nmom
"""
function subset_range!(full::TwoStream{FT}, start, n, subset) where {FT}

    ncol = get_ncol(full)
    nlay = get_nlay(full)
    ngpt = get_ngpt(full)
    if start < 1 || start + n - 1 > get_ncol(full)
        error("optical_props%subset: Asking for columns outside range")
    end

    init!(subset, full)

    if subset isa OneScalar # TODO: check logic
        alloc!(subset, n, nlay)
        extract_subset!(ncol, nlay, ngpt, full.tau, full.ssa, start, start + n - 1, subset.tau)
    elseif subset isa TwoStream
        alloc!(subset, n, nlay)
        extract_subset!(ncol, nlay, ngpt, full.tau, start, start + n - 1, subset.tau)
        extract_subset!(ncol, nlay, ngpt, full.ssa, start, start + n - 1, subset.ssa)
        extract_subset!(ncol, nlay, ngpt, full.g, start, start + n - 1, subset.g)
    elseif subset isa NStream
        if allocated(subset.p)
            nmom = get_nmom(subset)
        else
            nmom = 1
        end
        alloc!(subset, nmom, n, nlay)
        extract_subset!(ncol, nlay, ngpt, full.tau, start, start + n - 1, subset.tau)
        extract_subset!(ncol, nlay, ngpt, full.ssa, start, start + n - 1, subset.ssa)
        subset.p[1, 1:n, :, :] .= full.g[start:(start + n - 1), :, :]
        subset.p[2:end, :, :, :] .= FT(0) # TODO Verify this line
    else
        error("Uncaught case in subset_range!(full::TwoStream{FT}")
    end
end

"""
    subset_range!(full::NStream, start, n, subset)

    class(NStream), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(AbstractOpticalPropsArry), intent(inout) :: subset
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom
"""
function subset_range!(full::NStream, start, n, subset)

    ncol = get_ncol(full)
    nlay = get_nlay(full)
    ngpt = get_ngpt(full)
    if (start < 1 || start + n - 1 > get_ncol(full))
        error("optical_props%subset: Asking for columns outside range")
    end

    init!(subset, full)

    if subset isa OneScalar # TODO: check logic
        alloc!(subset, n, nlay)
        extract_subset!(ncol, nlay, ngpt, full.tau, full.ssa, start, start + n - 1, subset.tau)
    elseif subset isa TwoStream
        alloc!(subset, n, nlay)
        extract_subset!(ncol, nlay, ngpt, full.tau, start, start + n - 1, subset.tau)
        extract_subset!(ncol, nlay, ngpt, full.ssa, start, start + n - 1, subset.ssa)
        subset.g[1:n, :, :] .= full.p[1, start:(start + n - 1), :, :]
    elseif subset isa NStream
        alloc!(subset, nmom, n, nlay)
        extract_subset!(ncol, nlay, ngpt, full.tau, start, start + n - 1, subset.tau)
        extract_subset!(ncol, nlay, ngpt, full.ssa, start, start + n - 1, subset.ssa)
        extract_subset!(nmom, ncol, nlay, ngpt, full.p, start, start + n - 1, subset.p)
    else
        error("Uncaught case in subset_range!(full::NStream, start, n, subset)")
    end
end


"""
    class(NStream), intent(in   ) :: this
    integer                                     :: get_nmom
"""
get_nmom(this::NStream) = size(this.p, 1)


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
                t = tau_abs[igpt, ilay, icol] + tau_rayleigh[igpt, ilay, icol]
                tau[icol, ilay, igpt] = t
                if (t > FT(2) * realmin(FT))
                    ssa[icol, ilay, igpt] = tau_rayleigh[igpt, ilay, icol] / t
                else
                    ssa[icol, ilay, igpt] = FT(0)
                end
                for imom in 1:nmom
                    p[imom, icol, ilay, igpt] = FT(0)
                end
                if (nmom >= 2)
                    p[2, icol, ilay, igpt] = FT(0.1)
                end
            end
        end
    end
end
