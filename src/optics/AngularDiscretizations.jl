module AngularDiscretizations

using DocStringExtensions
using Adapt

export AngularDiscretization

"""
    AngularDiscretization{FT,FTA1D}

Weights and angle secants for first order (k=1) Gaussian quadrature.
Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
after Abramowitz & Stegun 1972, page 921

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AngularDiscretization{FT <: AbstractFloat, FTA1D <: AbstractArray{FT, 1}}
    "number of quadrature angles"
    n_gauss_angles::Int
    "quadrature secants / secant of propagation angle"
    gauss_Ds::FTA1D
    "quadrature weights"
    gauss_wts::FTA1D
end
Adapt.@adapt_structure AngularDiscretization

function AngularDiscretization(::Type{FT}, ::Type{DA}, n_gauss_angles::Int) where {FT <: AbstractFloat, DA}
    max_gauss_pts = 4
    @assert 1 ≤ n_gauss_angles ≤ max_gauss_pts
    if n_gauss_angles == 1
        gauss_Ds = DA{FT}([1.660000000000])
        gauss_wts = DA{FT}([0.500000000000])
    elseif n_gauss_angles == 2
        gauss_Ds = DA{FT}(FT[1.183503430000 2.816496550000])
        gauss_wts = DA{FT}([0.318041381700 0.181958618300])
    elseif n_gauss_angles == 3
        gauss_Ds = DA{FT}([1.097198580000 1.693385070000 4.709416300000])
        gauss_wts = DA{FT}([0.200931913700 0.229241106400 0.069826979900])
    else
        gauss_Ds = DA{FT}([1.060562570000 1.382825600000 2.401481790000 7.155130240000])
        gauss_wts = DA{FT}([0.135506913400 0.203464568000 0.129847547600 0.031180971000])
    end

    return AngularDiscretization{eltype(gauss_Ds), typeof(gauss_Ds)}(
        n_gauss_angles,
        gauss_Ds, # Diffusivity angle, not Gaussian angle
        gauss_wts,
    )
end

end #module
