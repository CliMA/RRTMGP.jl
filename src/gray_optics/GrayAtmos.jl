module GrayAtmos
using ..Device: array_type
using ..GrayAtmosphericStates
using DocStringExtensions
import GaussQuadrature
using ..GraySources
using ..GrayFluxes
using ..GrayAngularDiscretizations
using ..GrayOptics
using ..GrayBCs

using CLIMAParameters: Stefan
using CLIMAParameters.Planet: grav, R_d

export gas_optics_gray_atmos!,
    compute_gas_τs_gray_atmos!, tlay_to_tlev!, GrayRRTMGP

function gas_optics_gray_atmos!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I,B},
    optical_props::AbstractGrayOpticalProps{FT,FTA2D},
    source::AbstractGraySource{FT},
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    B<:Bool,
}

    compute_gas_τs_gray_atmos!(as, optical_props)
    compute_Planck_source!(source, as)
    if typeof(optical_props) == GrayTwoStream{FT,FTA2D}
        optical_props.ssa .= FT(0)
        optical_props.g .= FT(0)
    end
    return nothing
end

function compute_gas_τs_gray_atmos!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I,B},
    optical_props::AbstractGrayOpticalProps{FT,FTA2D},
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
    B<:Bool,
}
    ncol = as.ncol
    nlay = as.nlay
    top_at_1 = as.top_at_1
    p_lay = as.p_lay
    p_lev = as.p_lev
    α = as.α

    for icol = 1:ncol
        p0 = (top_at_1 ? as.p_lev[icol, end] : as.p_lev[icol, 1])
        d0 = as.d0[icol]
        for ilay = 1:nlay
            optical_props.τ[icol, ilay] = abs(
                (α * d0 * (p_lay[icol, ilay] ./ p0) .^ α ./ p_lay[icol, ilay]) * (p_lev[icol, ilay+1] - p_lev[icol, ilay]),
            )
        end
    end
    return nothing
end

function compute_Planck_source!(
    source::AbstractGraySource{FT},
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I,B},
) where {
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    B<:Bool,
}
    ncol = as.ncol
    nlay = as.nlay
    sfc_source = source.sfc_source
    lay_source = source.lay_source
    lev_source_inc = source.lev_source_inc
    lev_source_dec = source.lev_source_dec
    t_lay = as.t_lay
    t_lev = as.t_lev

    for icol = 1:ncol
        sfc_source[icol, 1] = Stefan() * as.t_sfc[icol]^4 / FT(π) # computing sfc_source
        for ilay = 1:nlay
            lay_source[icol, ilay, 1] =
                Stefan() * t_lay[icol, ilay]^FT(4) / FT(π) # computing lay_source

            lev_source_inc[icol, ilay, 1] =
                Stefan() * t_lev[icol, ilay+1]^4.0 / FT(π)
            lev_source_dec[icol, ilay, 1] =
                Stefan() * t_lev[icol, ilay]^4.0 / FT(π)
        end
    end
    return nothing
end

function tlay_to_tlev!(
    p_lay::FTA2D,
    p_lev::FTA2D,
    t_lay::FTA2D,
    t_lev::FTA2D,
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    ncol, nlay, nlev = size(t_lay, 1), size(t_lay, 2), size(t_lev, 2)
    #    for ilay in nlay:-1:2
    #        t_lev[1, ilay] = 2 * t_lay[1, ilay] - t_lev[1, ilay+1]
    #    end
    #    t_lev[1, 1] = 2 * t_lay[1, 1] - t_lev[1, 2]

    for icol = 1:ncol
        for j = 2:nlay-1
            t_lev[icol, j] =
                FT(1 / 3) * t_lay[icol, j-1] + FT(5 / 6) * t_lay[icol, j] -
                FT(1 / 6) * t_lay[icol, j+1]
        end
        t_lev[icol, 1] = FT(2) * t_lay[icol, 1] - t_lev[icol, 2]
        t_lev[icol, nlay] =
            FT(1 / 3) * t_lay[icol, nlay] + FT(5 / 6) * t_lay[icol, nlay-1] -
            FT(1 / 6) * t_lay[icol, nlay-2]
        t_lev[icol, nlay+1] = FT(2) * t_lay[icol, nlay] - t_lev[icol, nlay]
    end

end

struct GrayRRTMGP{
    FT<:AbstractFloat,
    I<:Int,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
    B<:Bool,
    OP<:AbstractGrayOpticalProps{FT,FTA2D},
    GS<:AbstractGraySource{FT},
    BC<:AbstractGrayBCs{FT,FTA1D},
}
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I,B}
    op::OP
    src::GS
    bcs::BC
    flux::GrayFlux{FT,FTA2D}
    angle_disc::AngularDiscretization{FT,FTA1D,I}

    function GrayRRTMGP(
        nlay::Int,
        lat::FTA1D,
        p0::FT,
        pe::FT,
        optical_props_constructor,
        sfc_emis::Array{FT,1},
        flux::GrayFlux{FT,FTA2D},
        n_gauss_angles::Int,
        param_set,
        top_at_1::Bool,
        stype,
        ::Type{DA},
    ) where {
        FT<:AbstractFloat,
        FTA1D<:AbstractArray{FT,1},
        FTA2D<:AbstractArray{FT,2},
        DA,
    }
        ncol, ngpt, nbnd = length(lat), 1, 1
        if stype == "lw"
            as = setup_gray_as_pr_grid(
                nlay,
                lat,
                p0,
                pe,
                param_set,
                top_at_1,
                DA,
            )
            optical_props = optical_props_constructor(FT, ncol, nlay, DA)
            sf = source_func_longwave_gray_atmos(
                FT,
                ncol,
                nlay,
                ngpt,
                optical_props_constructor,
                DA,
            )
            bc = GrayLwBCs(DA, sfc_emis)
            ang_disc = AngularDiscretization(FT, n_gauss_angles, DA)
            return new{
                FT,
                Int,
                DA{FT,1},
                DA{FT,2},
                DA{FT,3},
                Bool,
                typeof(optical_props),
                typeof(sf),
                typeof(bc),
            }(
                as,
                optical_props,
                sf,
                bc,
                flux,
                ang_disc,
            )
        else
            error("Shortwave solver not yet supported")
        end
    end

    function GrayRRTMGP(
        lat::FTA1D,
        p0::FT,
        zend::FT,
        ncls::Int,
        poly_order::Int,
        optical_props_constructor,
        sfc_emis::Array{FT,1},
        flux::GrayFlux{FT,FTA2D},
        n_gauss_angles::Int,
        param_set,
        top_at_1::Bool,
        stype,
        ::Type{DA},
    ) where {
        FT<:AbstractFloat,
        FTA1D<:AbstractArray{FT,1},
        FTA2D<:AbstractArray{FT,2},
        DA,
    }
        nlev = Int(ncls + 1 + (poly_order - 1) * ncls)
        nlay = Int(nlev - 1)
        ncol, ngpt, nbnd = length(lat), 1, 1
        if stype == "lw"
            as = setup_gray_as_alt_grid(
                lat,
                p0,
                zend,
                ncls,
                poly_order,
                param_set,
                top_at_1,
                DA,
            )
            optical_props = optical_props_constructor(FT, ncol, nlay, DA)
            sf = source_func_longwave_gray_atmos(
                FT,
                ncol,
                nlay,
                ngpt,
                optical_props_constructor,
                DA,
            )
            bc = GrayLwBCs(DA, sfc_emis)
            ang_disc = AngularDiscretization(FT, n_gauss_angles, DA)
            return new{
                FT,
                Int,
                DA{FT,1},
                DA{FT,2},
                DA{FT,3},
                Bool,
                typeof(optical_props),
                typeof(sf),
                typeof(bc),
            }(
                as,
                optical_props,
                sf,
                bc,
                flux,
                ang_disc,
            )
        else
            error("Shortwave solver not yet supported")
        end
    end

end

end
