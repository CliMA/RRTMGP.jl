module FortranIntrinsics
using DelimitedFiles
using Test

export fmerge,
       fminloc_wrapper,
       fmaxloc_wrapper,
       spread,
       spread_new,
       associated,
       allocated,
       pack,
       spacing,
       fint,
       freshape,
       present

"""
    squeeze(A)

Same as Python's squeeze, only using temporarily
TODO: Remove dependence on this function
"""
squeeze(A) = dropdims(A, dims=tuple([i for i in 1:length(size(A)) if size(A)[i]==1]...))

"""
    fint

Based on https://gcc.gnu.org/onlinedocs/gfortran/INT.html
"""
fint(x::Real) = abs(x) < 1 ? 0 : floor(x)

"""
    fmerge(t_source, f_source, mask)
Based on: http://gcc.gnu.org/onlinedocs/gfortran/MERGE.html
TODO: FIX IMPLEMENTATION
"""
@inline fmerge(t_source, f_source, mask) = mask ? t_source : f_source

"""
    fminloc(a, dim, mask=nothing)
Based on: https://gcc.gnu.org/onlinedocs/gfortran/MINLOC.html
"""
function fminloc(a; dim, mask=nothing)
  if mask==nothing
    return argmin(a, dims=dim)
  else
    tmp = deepcopy(a)
    tmp[.!mask] .= Inf
    return squeeze(argmin(tmp, dims=dim))
  end
end

function fminloc_wrapper(a; dim, mask=nothing)
  r = fminloc(a; dim=dim, mask=mask)
  res = [r[i][2] for i in 1:length(r)] # get index from CartesianIndex
  return res
end

"""
    fmaxloc(a, dim=nothing, mask=nothing)
Based on: https://gcc.gnu.org/onlinedocs/gfortran/MAXLOC.html
"""
function fmaxloc(a; dim, mask=nothing)
  if mask==nothing
    return argmax(a, dims=dim)
  else
    tmp = deepcopy(a)
    tmp[.!mask] .= -Inf
    return squeeze(argmax(tmp, dims=dim))
  end
end

function fmaxloc_wrapper(a; dim, mask=nothing)
  r = fmaxloc(a; dim=dim, mask=mask)
  res = [r[i][2] for i in 1:length(r)]  # get index from CartesianIndex
  return res
end

"""
    spread(a, dim=nothing, mask=nothing)
Based on: https://gcc.gnu.org/onlinedocs/gcc-4.4.0/gfortran/SPREAD.html
"""
function spread(source::Array{FT,N}, dim::Int, ncopies::Int) where {FT,N}
  @assert 1 <= dim <= N+1
  counts = [1 for i in 1:dim+1]
  counts[dim] = ncopies
  return squeeze(repeat(source, counts...))
end
function spread_new(source::Vector{FT}, dim::Int, ncopies::Int) where {FT}
  @assert 1 <= dim <= 2
  counts = [1 for i in 1:dim+1]
  counts[dim] = ncopies
  s = collect(size(source))
  mat = reshape(source, 1, length(source))
  return repeat(mat, counts...)
end

"""
    freshape
Based on: https://gcc.gnu.org/onlinedocs/gcc-4.3.6/gfortran/RESHAPE.html
"""
function freshape(source, shape; pad=nothing, order = nothing)
  if pad ≠ nothing && order ≠ nothing
    return reshape(source, shape...)
  elseif order ≠ nothing
    return permutedims(reshape(source, shape...), order)
  elseif pad ≠ nothing
    return reshape(source, shape...)
  else
    return reshape(source, shape...)
  end

end

"""
    spacing(x)

Based on: https://gnu.huihoo.org/gcc/gcc-4.4.5/gfortran/SPACING.html
TODO: Verify
"""
spacing(x) = nextfloat(x)-x


"""
    pack

Based on https://gcc.gnu.org/onlinedocs/gfortran/PACK.html
"""
function pack(arr::Array, mask::Array, v::Union{Nothing, Vector}=nothing)
  if v==nothing
    return reshape([x for (m,x) in zip(mask,arr) if m], count(mask))
  else
    @assert length(v) >= count(mask)
    return reshape([x for (m,x) in zip(mask,arr) if m], length(v))
  end
end

end
