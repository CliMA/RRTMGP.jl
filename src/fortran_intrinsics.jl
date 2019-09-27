module fortran_intrinsics

export fmerge,
       fminloc,
       fmaxloc,
       spread,
       associated,
       is_initialized,
       allocated,
       present

"""
    fmerge(t_source, f_source, mask)
Based on: http://gcc.gnu.org/onlinedocs/gfortran/MERGE.html
TODO: FIX IMPLEMENTATION
"""
fmerge(t_source, f_source, mask) = mask ? t_source : f_source

"""
    fminloc(a, dim, mask=nothing)
Based on: https://gcc.gnu.org/onlinedocs/gfortran/MINLOC.html
TODO: FIX IMPLEMENTATION
"""
fminloc(a; dim, mask=nothing) = mask==nothing ? argmin(a, dims=dim) : argmin(a[mask], dims=dim)

"""
    fmaxloc(a, dim=nothing, mask=nothing)
Based on: https://gcc.gnu.org/onlinedocs/gfortran/MAXLOC.html
TODO: FIX IMPLEMENTATION
"""
fmaxloc(a; dim, mask=nothing) = mask==nothing ? argmax(a, dims=dim) : argmax(a[mask], dims=dim)

"""
    fspread(a, dim=nothing, mask=nothing)
Based on: https://gcc.gnu.org/onlinedocs/gcc-4.4.0/gfortran/SPREAD.html
TODO: FIX IMPLEMENTATION
"""

function spread(source::T, dim::Int, ncopies::Int) where T<:Number
  @assert dim <= N+1
  Vector{T}([source for x in 1:ncopies])
end
function spread(source::Array{T,N}, dim::Int, ncopies::Int) where {T,N}
  @assert dim <= N+1
  counts = [1 for i in 1:dim]
  counts[dim] = ncopies
  return repeat(source, counts...)
end

associated(a::Array) = length(a)>0
allocated(a::Array) = a ≠ nothing
deallocate!(a) = (a = nothing)
is_initialized(a::Array) = allocated(a) && length(a)>0
present(arg) = arg ≠ nothing

"""
Based on: https://gnu.huihoo.org/gcc/gcc-4.4.5/gfortran/SPACING.html
TODO: Verify
"""
spacing(x) = nextfloat(x)-x

"""
    move_alloc!(from, to)
Based on: https://gcc.gnu.org/onlinedocs/gfortran/MOVE_005fALLOC.html
TODO: FIX IMPLEMENTATION
"""
function move_alloc!(from, to)
  to = deepcopy(from)
  from = nothing
end

# https://gcc.gnu.org/onlinedocs/gfortran/PACK.html
# pack()


trim(s) = strip(s)
len_trim(s) = length(strip(s))
maxval(v) = max(v)

end
