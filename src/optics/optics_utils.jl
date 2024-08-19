
@inline loc_lower(xi, Δx, n, x) = @inbounds max(min(unsafe_trunc(Int, (xi - x[1]) / Δx) + 1, n - 1), 1)

"""
    interp1d(xi::FT, x, y) where {FT<:AbstractFloat}

perform 1D linear interpolation.
"""
@inline function interp1d(xi::FT, x, y) where {FT <: AbstractFloat}
    @inbounds Δx = x[2] - x[1]
    n = length(x)
    loc = loc_lower(xi, Δx, n, x)
    @inbounds factor = (xi - x[loc]) / Δx
    return @inbounds (y[loc] * (FT(1) - factor) + y[loc + 1] * factor)
end

@inline function interp1d_loc_factor(xi::FT, x::AbstractArray{FT, 1}) where {FT <: AbstractFloat}
    @inbounds Δx = x[2] - x[1]
    n = length(x)
    loc = loc_lower(xi, Δx, n, x)
    @inbounds factor = (xi - x[loc]) / Δx
    return loc, factor
end

"""
    interp2d(
        fη1::FT,
        fη2::FT,
        ftemp::FT,
        coeff::FTA2D,
        jη1::Int,
        jη2::Int,
        jtemp::Int,
    ) where {FT,FTA2D}

Perform 2D linear interpolation.

`fminor[1, 1] = (FT(1) - fη1) * (1-ftemp)`

`fminor[2, 1] = fη1 * (1-ftemp)`

`fminor[1, 2] = (FT(1) - fη2) * ftemp`

`fminor[2, 2] = fη2 * ftemp`
"""
@inline function interp2d(fη1::FT, fη2::FT, ftemp::FT, coeff::FTA2D, jη1::Int, jη2::Int, jtemp::Int) where {FT, FTA2D}
    return @inbounds (FT(1) - fη1) * (1 - ftemp) * coeff[jη1, jtemp] +
                     fη1 * (1 - ftemp) * coeff[jη1 + 1, jtemp] +
                     (FT(1) - fη2) * ftemp * coeff[jη2, jtemp + 1] +
                     fη2 * ftemp * coeff[jη2 + 1, jtemp + 1]
end

"""
    interp3d(
        jη1::Int,
        jη2::Int,
        fη1::FT,
        fη2::FT,
        jtemp::Int,
        ftemp::FT,
        jpresst::Int,
        fpress::FT,
        coeff::FTA3D,
        s1::FT = FT(1),
        s2::FT = FT(1),
    ) where {FT<:AbstractFloat,FTA4D<:AbstractArray{FT,4}}

Perform 3D linear interpolation.

`fmajor[1, 1, itemp] = (FT(1) - fpress) * fminor[1, itemp]`

`fmajor[2, 1, itemp] = (FT(1) - fpress) * fminor[2, itemp]`

`fmajor[1, 2, itemp] = fpress * fminor[1, itemp]`

`fmajor[2, 2, itemp] = fpress * fminor[2, itemp]`

where,

`fminor[1, itemp=1] = (FT(1) - fη1) * (1-ftemp)`

`fminor[2, itemp=1] = fη1 * (1-ftemp)`

`fminor[1, itemp=2] = (FT(1) - fη2) * ftemp`

`fminor[2, itemp=2] = fη2 * ftemp`

"""
@inline function interp3d(
    jη1::Int,
    jη2::Int,
    fη1::FT,
    fη2::FT,
    jtemp::Int,
    ftemp::FT,
    jpresst::Int,
    fpress::FT,
    coeff::FTA3D,
    s1::FT = FT(1),
    s2::FT = FT(1),
) where {FT, FTA3D}
    omftemp = FT(1) - ftemp
    omfpress = FT(1) - fpress
    omfη1 = FT(1) - fη1
    omfη2 = FT(1) - fη2
    return @inbounds s1 * (
        omfpress * (omftemp * (omfη1 * coeff[jη1, jpresst - 1, jtemp] + fη1 * coeff[jη1 + 1, jpresst - 1, jtemp])) +
        fpress * (omftemp * (omfη1 * coeff[jη1, jpresst, jtemp] + fη1 * coeff[jη1 + 1, jpresst, jtemp]))
    ) +
                     s2 * (
        omfpress *
        (ftemp * (omfη2 * coeff[jη2, jpresst - 1, jtemp + 1] + fη2 * coeff[jη2 + 1, jpresst - 1, jtemp + 1])) +
        fpress * (ftemp * (omfη2 * coeff[jη2, jpresst, jtemp + 1] + fη2 * coeff[jη2 + 1, jpresst, jtemp + 1]))
    )
end
"""
    increment_2stream!(τ1::FT, ssa1::FT, g1::FT, τ2::FT, ssa2::FT, g2::FT) where {FT}

Increment TwoStream optical properties `τ1`, `ssa1` and `g1` 
with `τ2`, `ssa2` and `g2`. Here `τ` is the optical thickness,
`ssa` is the single-scattering albedo, and `g` is the symmetry parameter.
"""
function increment_2stream(τ1::FT, ssa1::FT, g1::FT, τ2::FT, ssa2::FT, g2::FT) where {FT}
    τ = τ1 + τ2
    ssa = τ1 * ssa1 + τ2 * ssa2
    ssag = (τ1 * ssa1 * g1 + τ2 * ssa2 * g2) / max(eps(FT), ssa)
    ssa /= max(eps(FT), τ)
    return τ, ssa, ssag
end
"""
    delta_scale(τ, ssa, g)

delta-scale two stream optical properties.
"""
function delta_scale(τ, ssa, g)
    FT = typeof(τ)
    f = g * g
    wf = ssa * f
    τ_s = (FT(1) - wf) * τ
    ssa_s = (ssa - wf) / max(eps(FT), FT(1) - wf)
    g_s = (g - f) / max(eps(FT), FT(1) - f)
    return (τ_s, ssa_s, g_s)
end
