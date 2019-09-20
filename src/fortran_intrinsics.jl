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
spread(source, dim, ncopies) = ...

associated(a::Array) = length(a)>0
is_initialized(a::Array) = length(a)>0
allocated(a::Array) = is_initialized(a)
present(arg) = arg ≠ nothing

end