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
    gauss_Ds::Array{FT,1}
    "quadrature weights"
    gauss_wts::Array{FT,1}
    function GaussQuadrature(
        ::Type{FT},
        n_gauss_angles::I,
    ) where {FT<:AbstractFloat,I<:Int}
        max_gauss_pts = I(4)
        @assert 1 ≤ n_gauss_angles ≤ max_gauss_pts

        if n_gauss_angles == 1
            Ds = [1.660000000000] # Diffusivity angle, not Gaussian angle
            wts = [0.500000000000]
        elseif n_gauss_angles == 2
            Ds = [1.183503430000 2.816496550000]
            wts = [0.318041381700 0.181958618300]
        elseif n_gauss_angles == 3
            Ds = [1.097198580000 1.693385070000 4.709416300000]
            wts = [0.200931913700 0.229241106400 0.069826979900]
        else
            Ds = [1.060562570000 1.382825600000 2.401481790000 7.155130240000]
            wts = [0.135506913400 0.203464568000 0.129847547600 0.031180971000]
        end

        return new{FT,I}(n_gauss_angles, Ds, wts)
    end
end

end #module
