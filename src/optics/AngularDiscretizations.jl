module AngularDiscretizations

using DocStringExtensions
using Adapt

export AngularDiscretization

"""
    AngularDiscretization{FT,FTA1D}

Weights and angle secants for "Gauss-Jacobi-5" quadrature.
Values from Table 1, R. J. Hogan 2023, doi:10.1002/qj.4598

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
        gauss_Ds = DA{FT}(FT(1) ./ [0.6096748751])
        gauss_wts = DA{FT}([1])
    elseif n_gauss_angles == 2
        gauss_Ds = DA{FT}(FT(1) ./ [0.2509907356, 0.7908473988])
        gauss_wts = DA{FT}([0.2300253764, 0.7699746236])
    elseif n_gauss_angles == 3
        gauss_Ds = DA{FT}(FT(1) ./ [0.1024922169, 0.4417960320, 0.8633751621])
        gauss_wts = DA{FT}([0.0437820218, 0.3875796738, 0.5686383044])
    else
        gauss_Ds = DA{FT}(FT(1) ./ [0.0454586727, 0.2322334416, 0.5740198775, 0.9030775973])
        gauss_wts = DA{FT}([0.0092068785, 0.1285704278, 0.4323381850, 0.4298845087])
    end

    return AngularDiscretization{eltype(gauss_Ds), typeof(gauss_Ds)}(
        n_gauss_angles,
        gauss_Ds, # Diffusivity angle, not Gaussian angle
        gauss_wts,
    )
end

end #module
