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
       freshape,
       test_data,
       present

test_data(x::OffsetArray, args...) = test_data(reshape([x[:]...], size(x)...), args...)

function test_data(x, name, tol = nothing)
  if eltype(x) <: AbstractFloat
    FT = Float64
    if eltype(x) ≠ Float64
      @show eltype(x)
      error("Data not Float64: $(name)")
    end
    tol==nothing && (tol = eps(FT))
  elseif eltype(x) <: Integer
    FT = Int64
    tol = eps(Float64)
  else
    @show eltype(x)
    error("Data not Float64: $(name)")
  end
  x_correct, s = readdlm(joinpath("/Users","charliekawczynski","TestData",name)*".dat", '\t', FT, header=true, '\n')
  s = parse.(Int, filter(i-> i≠"", split(strip(s[1]), " ")))
  x_correct = reshape(x_correct, s...)
  @assert all(size(x_correct) .== size(x))
  x_max = max(x_correct...)
  ad = abs.(x .- x_correct)
  cond = ad .< tol
  L = all(cond)
  if !L
    println("\n\n\n\n\n")
    println(" ***************************** Failed comparison: "*name)
    @show size(ad)
    @show max(x_correct...)
    @show sum(ad)
    @show max(ad...)
    @show x_correct[1:10]
    @show x[1:10]
    @show ad[1:10]
    # @show ad
    @show count(cond)
    @show length(cond)
    @show count(cond)/length(cond)
    @show sum(abs.(x))
    @show sum(abs.(x_correct))
    @show max(abs.(x .- x_correct)...)
    @show sum(abs.(x .- x_correct))
    error("Failed comparison in test_data")
  else
    if !occursin("gas_conc",name)
      println(" ----------------------------- Successful comparison: "*name)
    end
  end
  @test L
end


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
fint(x::Integer) = x
fint(x::Real) = abs(x) < 1 ? 0 : floor(x)

"""
    fmerge(t_source, f_source, mask)
Based on: http://gcc.gnu.org/onlinedocs/gfortran/MERGE.html
TODO: FIX IMPLEMENTATION
"""
fmerge(t_source, f_source, mask) = mask ? t_source : f_source

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
function pack(arr::Array, mask::Bool, v::Union{Nothing, Vector}=nothing)
  if v==nothing
    return reshape([x for (m,x) in zip(mask,arr) if m], count(mask))
  else
    @assert length(v) >= length(arr)
    return reshape([x for (m,x) in zip(mask,arr) if m], length(v))
  end
end


associated(a::Nothing) = false
associated(a::AbstractArray) = length(a)>0
allocated(a::Array) = a ≠ nothing
allocated(a::Nothing) = false
deallocate!(a) = (a = nothing)
is_initialized(a::Array) = allocated(a) && length(a)>0
is_initialized(a::Nothing) = false
present(arg) = arg ≠ nothing
trim(s::AbstractString) = strip(s)
trim(s::Char) = s
len_trim(s) = length(strip(s))
maxval(v) = max(v...)

end
