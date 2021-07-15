
@inline loc_lower(xi, Δx, n, x) =
    @inbounds max(min(unsafe_trunc(Int, (xi - x[1]) / Δx) + 1, n - 1), 1)

"""
    interp1d(xi::FT, x, y, col) where {FT<:AbstractFloat}

perform 1D linear interpolation.
"""
@inline function interp1d(xi::FT, x, y, col) where {FT<:AbstractFloat}
    @inbounds Δx = x[2] - x[1]
    n = length(x)
    loc = loc_lower(xi, Δx, n, x)
    @inbounds factor = (xi - x[loc]) / Δx
    return @inbounds (y[loc, col] * (FT(1) - factor) + y[loc+1, col] * factor)
end

"""
    interp2d(
        fη1::FT,
        fη2::FT,
        ftemp::FT,
        coeff::FTA3D,
        igpt::I,
        jη1::I,
        jη2::I,
        jtemp::I,
    ) where {FT<:AbstractFloat,FTA3D<:AbstractArray{FT,3},I<:Int}

Perform 2D linear interpolation.

`fminor[1, 1] = (FT(1) - fη1) * (1-ftemp)`

`fminor[2, 1] = fη1 * (1-ftemp)`

`fminor[1, 2] = (FT(1) - fη2) * ftemp`

`fminor[2, 2] = fη2 * ftemp`
"""
@inline function interp2d(
    fη1::FT,
    fη2::FT,
    ftemp::FT,
    coeff::FTA3D,
    igpt::I,
    jη1::I,
    jη2::I,
    jtemp::I,
) where {FT<:AbstractFloat,FTA3D<:AbstractArray{FT,3},I<:Int}
    return @inbounds (FT(1) - fη1) * (1 - ftemp) * coeff[igpt, jη1, jtemp] +
                     fη1 * (1 - ftemp) * coeff[igpt, jη1+1, jtemp] +
                     (FT(1) - fη2) * ftemp * coeff[igpt, jη2, jtemp+1] +
                     fη2 * ftemp * coeff[igpt, jη2+1, jtemp+1]
end

"""
    interp3d(
        jη1::I,
        jη2::I,
        fη1::FT,
        fη2::FT,
        jtemp::I,
        ftemp::FT,
        jpresst::I,
        fpress::FT,
        coeff::FTA4D,
        igpt::I,
        s1::FT = FT(1),
        s2::FT = FT(1),
    ) where {FT<:AbstractFloat,FTA4D<:AbstractArray{FT,4},I<:Int}

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
    jη1::I,
    jη2::I,
    fη1::FT,
    fη2::FT,
    jtemp::I,
    ftemp::FT,
    jpresst::I,
    fpress::FT,
    coeff::FTA4D,
    igpt::I,
    s1::FT = FT(1),
    s2::FT = FT(1),
) where {FT<:AbstractFloat,FTA4D<:AbstractArray{FT,4},I<:Int}
    omftemp = FT(1) - ftemp
    omfpress = FT(1) - fpress
    omfη1 = FT(1) - fη1
    omfη2 = FT(1) - fη2
    return @inbounds s1 * (
        omfpress * (
            omftemp * (
                omfη1 * coeff[igpt, jη1, jpresst-1, jtemp] +
                fη1 * coeff[igpt, jη1+1, jpresst-1, jtemp]
            )
        ) +
        fpress * (
            omftemp * (
                omfη1 * coeff[igpt, jη1, jpresst, jtemp] +
                fη1 * coeff[igpt, jη1+1, jpresst, jtemp]
            )
        )
    ) +
                     s2 * (
        omfpress * (
            ftemp * (
                omfη2 * coeff[igpt, jη2, jpresst-1, jtemp+1] +
                fη2 * coeff[igpt, jη2+1, jpresst-1, jtemp+1]
            )
        ) +
        fpress * (
            ftemp * (
                omfη2 * coeff[igpt, jη2, jpresst, jtemp+1] +
                fη2 * coeff[igpt, jη2+1, jpresst, jtemp+1]
            )
        )
    )
end
"""
    increment!(
        op::TwoStream{FT},
        τ2,
        ssa2,
        g2,
        glay,
        gcol,
        igpt,
    ) where {FT<:AbstractFloat}

Increment TwoStream optical properties `op` for layer `glay` 
and column `gcol` with optical thickness `τ2`, single scattering
albedo `ssa2` and symmetry parameter `g2`.
"""
function increment!(
    op::TwoStream{FT},
    τ2,
    ssa2,
    g2,
    glay,
    gcol,
    igpt,
) where {FT<:AbstractFloat}
    τ1, ssa1, g1 = op.τ[glay, gcol], op.ssa[glay, gcol], op.g[glay, gcol]
    τ = τ1 + τ2
    ssa = τ1 * ssa1 + τ2 * ssa2
    ssag = (τ1 * ssa1 * g1 + τ2 * ssa2 * g2) / max(eps(FT), ssa)
    ssa /= max(eps(FT), τ)

    op.τ[glay, gcol] = τ
    op.ssa[glay, gcol] = ssa
    op.g[glay, gcol] = ssag
    return nothing
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
