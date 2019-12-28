module AngularDiscretizations

export AbstractAngularDiscretization
export GaussQuadrature

abstract type AbstractAngularDiscretization{FT<:AbstractFloat, I<:Int} end

"""
    GaussQuadrature{FT, I}

Weights and angle secants for first order (k=1) Gaussian quadrature.
Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
after Abramowitz & Stegun 1972, page 921

"""
struct GaussQuadrature{FT, I} <: AbstractAngularDiscretization{FT, I}
  "number of quadrature angles"
  n_gauss_angles::I
  "quadrature secants / secant of propagation angle"
  gauss_Ds::Array{FT,2}
  "quadrature weights"
  gauss_wts::Array{FT,2}
  function GaussQuadrature(::Type{FT}, n_gauss_angles::I) where {FT<:AbstractFloat,I<:Int}
    max_gauss_pts = I(4)
    # Diffusivity angle, not Gaussian angle
    gauss_Ds  = reshape(FT[1.66000000, 0.00000000, 0.00000000, 0.00000000,
                           1.18350343, 2.81649655, 0.00000000, 0.00000000,
                           1.09719858, 1.69338507, 4.70941630, 0.00000000,
                           1.06056257, 1.38282560, 2.40148179, 7.15513024],
                           max_gauss_pts, max_gauss_pts)

    gauss_wts = reshape(FT[0.5000000000, 0.0000000000, 0.0000000000, 0.0000000000,
                           0.3180413817, 0.1819586183, 0.0000000000, 0.0000000000,
                           0.2009319137, 0.2292411064, 0.0698269799, 0.0000000000,
                           0.1355069134, 0.2034645680, 0.1298475476, 0.0311809710],
                           max_gauss_pts, max_gauss_pts)
    @assert !(n_gauss_angles > max_gauss_pts)
    @assert !(n_gauss_angles < 1)
    return new{FT,I}(n_gauss_angles, gauss_Ds, gauss_wts)
  end
end

end #module