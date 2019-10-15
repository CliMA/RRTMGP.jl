module fortran_intrinsics
using DelimitedFiles
using OffsetArrays
using Test

export fmerge,
       fminloc,
       fmaxloc,
       fminloc_wrapper,
       fmaxloc_wrapper,
       spread,
       associated,
       is_initialized,
       allocated,
       trim,
       pack,
       move_alloc!,
       maxval,
       spacing,
       fint,
       linear_show,
       freshape,
       test_data,
       present

"""
    squeeze(A)

Same as Python's squeeze, only using temporarily
TODO: Remove dependence on this function
"""
squeeze(A) = dropdims(A, dims=tuple([i for i in 1:length(size(A)) if size(A)[i]==1]...))

# https://gcc.gnu.org/onlinedocs/gfortran/INT.html
fint(x::Integer) = x
fint(x::Real) = abs(x) < 1 ? 0 : floor(x)

test_data(x::OffsetArray, name) = test_data(reshape([x[:]...], size(x)...), name)

function test_data(x, name)
  x_correct, s = readdlm(joinpath("..","TestData",name)*".dat", '\t', Float64, header=true, '\n')
  s = parse.(Int, filter(i-> i≠"", split(strip(s[1]), " ")))
  x_correct = reshape(x_correct, s...)
  @assert all(size(x_correct) .== size(x))
  FT = Float32
  x_max = max(x_correct...)
  d = x .- x_correct
  if x_max>1
    d = d ./x_max
  end
  ad = abs.(d) .< eps(Float32)
  # L = sum(abs.(d)) < eps(Float32)
  L = all(ad)
  if !L
    println("\n\n\n\n\n")
    println(" ***************************** Failed comparison: "*name)
    @show max(x_correct...)
    @show sum(abs.(d))
    @show max(abs.(d)...)
    @show x_correct[1:10]
    @show x[1:10]
    @show d[1:10]
    # @show count(ad)
    # @show count(ad)/length(ad)
    @show sum(abs.(x))
    @show sum(abs.(x_correct))
    @show max(abs.(x .- x_correct)...)
    @show sum(abs.(x .- x_correct))
  else
    if !occursin("gas_conc",name)
      println(" ----------------------------- Successful comparison: "*name)
    end
  end
  @test L
end

function linear_show(x, name)
  println("-------------------------------- Variable: ", name)
  @show x[:]
  # for i in eachindex(x)
  #   @show i, x[i]
  # end
  println("--------------------------------")
end

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
TODO: FIX IMPLEMENTATION
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
    fspread(a, dim=nothing, mask=nothing)
Based on: https://gcc.gnu.org/onlinedocs/gcc-4.4.0/gfortran/SPREAD.html
TODO: FIX IMPLEMENTATION
"""

function spread(source::FT, dim::Int, ncopies::Int) where FT<:Number
  @assert dim <= N+1
  Vector{FT}([source for x in 1:ncopies])
end

function spread(source::Array{FT,N}, dim::Int, ncopies::Int) where {FT,N}
  @assert 1 <= dim <= N+1
  counts = [1 for i in 1:dim+1]
  counts[dim] = ncopies
  return squeeze(repeat(source, counts...))
end

# function spread(source::Array{FT,N}, dim::Int, ncopies::Int) where {FT,N}
#   @assert 1 <= dim <= N+1
#   # if dim==2
#   #   x = reshape(source, 1, length(source))
#   #   return repeat(x, 1, ncopies)
#   # else
#     counts = [1 for i in 1:dim+1]
#     counts[dim] = ncopies
#     # @show counts
#     return squeeze(repeat(source, counts...))
#   # end
# end
# function spread(source::FT, dim::Int, ncopies::Int) where FT<:Number
#   @assert dim <= N+1
#   Vector{FT}([source for x in 1:ncopies])
# end
# function spread(source::Array{FT,N}, dim::Int, ncopies::Int) where {FT,N}
#   @assert dim <= N+1
#   counts = [1 for i in 1:dim]
#   counts[dim] = ncopies
#   return repeat(source, counts...)
# end

function freshape(source, shape; pad=nothing, order = nothing)
  if pad ≠ nothing && order ≠ nothing
    return reshape(source, shape...)
  elseif order ≠ nothing
    return permutedims(reshape(source, shape...), order)
    # return reshape(source, shape...), order)
  elseif pad ≠ nothing
    return reshape(source, shape...)
  else
    return reshape(source, shape...)
  end

end

associated(a::Nothing) = false
associated(a::Array) = length(a)>0
allocated(a::Array) = a ≠ nothing
allocated(a::Nothing) = false
deallocate!(a) = (a = nothing)
is_initialized(a::Array) = allocated(a) && length(a)>0
is_initialized(a::Nothing) = false
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
# function move_alloc!(from, to)
#   to = deepcopy(from)
#   from = nothing
# end

# https://gcc.gnu.org/onlinedocs/gfortran/PACK.html
function pack(arr::Array, mask::Array, v::Union{Nothing, Vector}=nothing)
  if v==nothing
    return reshape([x for (m,x) in zip(mask,arr) if m], count(mask))
  else
    @assert length(v) >= count(mask)
    return reshape([x for (m,x) in zip(mask,arr) if m], length(v))
  end
end
function pack(arr::Array, mask::Bool, v::Union{Nothing, Vector}=nothing)
  if v==nothing
    return reshape([x for (m,x) in zip(mask,arr) if m], count(mask))
  else
    @assert length(v) >= length(arr)
    return reshape([x for (m,x) in zip(mask,arr) if m], length(v))
  end
end


trim(s::AbstractString) = strip(s)
trim(s::Char) = s
len_trim(s) = length(strip(s))
maxval(v) = max(v...)
# function maxval(v, dim=nothing, mask=nothing)
#   if dim==nothing
#     maximum(v, dims=dim)
#   else
#     maximum(v, dims=dim)
#   end
# end
# function maxval(v, mask=nothing)
#   if dim==nothing
#     max(v...)
#   else
#     max(v...)
#   end
# end

end
