module GraySources

using DocStringExtensions
using ..Device: array_type
using ..GrayOptics: GrayOneScalar, GrayTwoStream, AbstractGrayOpticalProps

export AbstractGraySource,
    GraySourceLWNoScat, GraySourceLW2Str, source_func_longwave_gray_atmos

abstract type AbstractGraySource{FT} end

"""
    GraySourceLWNoScat{FTA2D, FTA3D} <: AbstractGraySource{FTA1D, FTA2D, FTA3D}

Longwave sources: computed at layer center, layer edges, 
and at the surface for no scattering calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct GraySourceLWNoScat{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
} <: AbstractGraySource{FT}
    "Planck source at layer average temperature `[W/m2]` `(ncol, nlay, ngpt)`"
    lay_source::FTA3D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(ncol, nlay+1, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA3D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(ncol, nlay+1, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA3D
    "Surface source"
    sfc_source::FTA2D
    "temporary storage array used in no scattering GPU calculations"
    source_up::FTA2D
    "temporary storage array used in no scattering GPU calculations"
    source_dn::FTA2D
    "transmissivity, a temporary storage array used in no scattering GPU calculations"
    trans::FTA2D
end

"""
    GraySourceLW2Str{FTA1D, FTA2D, FTA3D} <: AbstractGraySource{FTA1D, FTA2D, FTA3D}

Longwave sources: computed at layer center, layer edges, 
and at the surface for 2-stream calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct GraySourceLW2Str{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
} <: AbstractGraySource{FT}
    "Planck source at layer average temperature `[W/m2]` `(ncol, nlay, ngpt)`"
    lay_source::FTA3D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(ncol, nlay+1, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA3D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(ncol, nlay+1, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA3D
    "Surface source"
    sfc_source::FTA2D
    "level source `[W/m2]` `(ncol, nlay+1, ngpt)`, used in 2 stream calculations"
    lev_source::FTA2D
    "temporary storage array, used in 2 stream calculations"
    src_up::FTA1D
    "temporary storage array, used in 2 stream calculations"
    src_dn::FTA1D
    "temporary storage array, used in 2 stream calculations"
    Rdif::FTA1D
    "temporary storage array, used in 2 stream calculations"
    Tdif::FTA1D
    "temporary storage array, used in 2 stream calculations"
    albedo::FTA1D
    "temporary storage array, used in 2 stream calculations"
    src::FTA1D
end

function source_func_longwave_gray_atmos(
    ::Type{FT},
    ncol::I,
    nlay::I,
    ngpt::I,
    ::Type{OPC},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,OPC<:AbstractGrayOpticalProps,DA}
    println("DA = $DA")
    lay_source = DA{FT}(undef, ncol, nlay, ngpt)
    lev_source_inc = DA{FT}(undef, ncol, nlay, ngpt)
    lev_source_dec = DA{FT}(undef, ncol, nlay, ngpt)
    sfc_source = DA{FT}(undef, ncol, ngpt)

    if OPC == GrayOneScalar
        source_up = DA{FT}(undef, ncol, nlay)
        source_dn = DA{FT}(undef, ncol, nlay)
        trans = DA{FT}(undef, ncol, nlay)
        return GraySourceLWNoScat{
            eltype(lay_source),
            DA{FT,1},
            DA{FT,2},
            typeof(lay_source),
        }(
            lay_source,
            lev_source_inc,
            lev_source_dec,
            sfc_source,
            source_up,
            source_dn,
            trans,
        )
    else
        lev_source = DA{FT}(undef, ncol, nlay + 1)
        src_up = DA{FT}(undef, nlay)
        src_dn = DA{FT}(undef, nlay)
        Rdif = DA{FT}(undef, nlay)
        Tdif = DA{FT}(undef, nlay)
        albedo = DA{FT}(undef, nlay + 1)
        src = DA{FT}(undef, nlay + 1)
        return GraySourceLW2Str{FT,DA{FT,1},DA{FT,2},DA{FT,3}}(
            lay_source,
            lev_source_inc,
            lev_source_dec,
            sfc_source,
            lev_source,
            src_up,
            src_dn,
            Rdif,
            Tdif,
            albedo,
            src,
        )
    end
end

end
