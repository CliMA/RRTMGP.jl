module AngularDiscretizations

export AbstractAngularDiscretization
export GaussQuadrature

abstract type AbstractAngularDiscretization{FT<:AbstractFloat,I<:Int} end

"""
    GaussQuadrature{FT, I}

Weights and angle secants for first order (k=1) Gaussian quadrature.
Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
after Abramowitz & Stegun 1972, page 921

"""
struct GaussQuadrature{FT,I} <: AbstractAngularDiscretization{FT,I}
    "number of quadrature angles"
    n_gauss_angles::I
    "quadrature secants / secant of propagation angle"
    gauss_Ds::Array{FT,2}
    "quadrature weights"
    gauss_wts::Array{FT,2}
    function GaussQuadrature(
        ::Type{FT},
        n_gauss_angles::I,
    ) where {FT<:AbstractFloat,I<:Int}
        max_gauss_pts = I(4)
    # Diffusivity angle, not Gaussian angle
        gauss_Ds = reshape(
            FT[
               1.66000000,
               0.00000000,
               0.00000000,
               0.00000000,
               1.18350343,
               2.81649655,
               0.00000000,
               0.00000000,
               1.09719858,
               1.69338507,
               4.70941630,
               0.00000000,
               1.06056257,
               1.38282560,
               2.40148179,
               7.15513024,
            ],
            max_gauss_pts,
            max_gauss_pts,
        )

        table2 = (
            FT[
               0.5000000000 60.0000000 1.0000000000 0.6666666667 48.1896852 0.5000000000 # 1 points
               0.2113248654 77.7999960 0.5000000000 0.3550510257 69.2034283 0.1819586183 # 2 points
               0.7886751346 37.9381274 0.5000000000 0.8449489743 32.3335307 0.3180413817
               0.1127016654 83.5289027 0.2777777778 0.2123405382 77.7404509 0.0698269799 # 3 points
               0.5000000000 60.0000000 0.4444444444 0.5905331356 53.8051497 0.2292411064
               0.8872983346 27.4643042 0.2777777778 0.9114120405 24.2987798 0.2009319137
               0.0694318442 86.0186452 0.1739274226 0.1397598643 81.9660491 0.0311809710 # 4 points
               0.3300094782 70.7306492 0.3260725774 0.4164095676 65.3918850 0.1298475476
               0.6699905218 47.9336668 0.3260725774 0.7231569864 43.6842540 0.2034645680
               0.9305681558 21.4764462 0.1739274226 0.9428958039 19.4563020 0.1355069134
            ]
        )
        moment_k0 = table2[:, 1:3]
        moment_k1 = table2[:, 4:6]
        gauss_wts = reshape(
            FT[
               moment_k1[1, end],
               zeros(FT, 3)...,
               reverse(moment_k1[2:3, end])...,
               zeros(FT, 2)...,
               reverse(moment_k1[4:6, end])...,
               zeros(FT, 1)...,
               reverse(moment_k1[7:10, end])...,
            ],
            max_gauss_pts,
            max_gauss_pts,
        )

        @assert !(n_gauss_angles > max_gauss_pts)
        @assert !(n_gauss_angles < 1)
        return new{FT,I}(n_gauss_angles, gauss_Ds, gauss_wts)
    end
end

end #module
