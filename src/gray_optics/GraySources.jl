module GraySources

using DocStringExtensions
using ..Device: array_type
using ..GrayOptics: GrayOneScalar, GrayTwoStream, AbstractGrayOpticalProps

export AbstractGraySource,
    GraySourceLWNoScat,
    GraySourceLW2Str,
    source_func_longwave_gray_atmos,
    GraySourceSW2Str,
    source_func_shortwave_gray_atmos

abstract type AbstractGraySource{FT} end

"""
    GraySourceLWNoScat{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
} <: AbstractGraySource{FT}

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
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol, ngpt)`"
    lay_source::FTA3D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA3D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
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
    GraySourceLW2Str{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
} <: AbstractGraySource{FT}

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
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol, ngpt)`"
    lay_source::FTA3D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA3D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA3D
    "Surface source"
    sfc_source::FTA2D
    "level source `[W/m2]` `(nlay+1, ncol, ngpt)`, used in 2 stream calculations"
    lev_source::FTA2D
    "temporary storage array, used in 2 stream calculations"
    src_up::FTA2D
    "temporary storage array, used in 2 stream calculations"
    src_dn::FTA2D
    "temporary storage array, used in 2 stream calculations"
    Rdif::FTA2D
    "temporary storage array, used in 2 stream calculations"
    Tdif::FTA2D
    "temporary storage array, used in 2 stream calculations"
    albedo::FTA2D
    "temporary storage array, used in 2 stream calculations"
    src::FTA2D
end

function source_func_longwave_gray_atmos(
    ::Type{FT},
    ncol::I,
    nlay::I,
    ngpt::I,
    ::Type{OPC},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,OPC<:AbstractGrayOpticalProps,DA}
    lay_source = DA{FT}(undef, nlay, ncol, ngpt)
    lev_source_inc = DA{FT}(undef, nlay, ncol, ngpt)
    lev_source_dec = DA{FT}(undef, nlay, ncol, ngpt)
    sfc_source = DA{FT}(undef, ncol, ngpt)

    if OPC == GrayOneScalar
        source_up = DA{FT}(undef, nlay, ncol)
        source_dn = DA{FT}(undef, nlay, ncol)
        trans = DA{FT}(undef, nlay, ncol)
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
        lev_source = DA{FT}(undef, nlay + 1, ncol)
        src_up = DA{FT}(undef, nlay, ncol)
        src_dn = DA{FT}(undef, nlay, ncol)
        Rdif = DA{FT}(undef, nlay, ncol)
        Tdif = DA{FT}(undef, nlay, ncol)
        albedo = DA{FT}(undef, nlay + 1, ncol)
        src = DA{FT}(undef, nlay + 1, ncol)
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

"""
    GraySourceSW2Str{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    FTA3D<:AbstractArray{FT,3},
} <: AbstractGraySource{FT}


# Fields

$(DocStringExtensions.FIELDS)
"""
struct GraySourceSW2Str{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
} <: AbstractGraySource{FT}
    "Surface source"
    src_sfc::FTA1D
    "Surface albedo"
    sfc_albedo::FTA1D
    "temporary storage array, used in 2 stream calculations"
    src_up::FTA2D
    "temporary storage array, used in 2 stream calculations"
    src_dn::FTA2D
    "temporary storage array, used in 2 stream calculations"
    Rdif::FTA2D
    "temporary storage array, used in 2 stream calculations"
    Tdif::FTA2D
    "temporary storage array, used in 2 stream calculations"
    Rdir::FTA2D
    "temporary storage array, used in 2 stream calculations"
    Tdir::FTA2D
    "temporary storage array, used in 2 stream calculations"
    Tnoscat::FTA2D
end

function source_func_shortwave_gray_atmos(
    ::Type{FT},
    ncol::I,
    nlay::I,
    ::Type{OPC},
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,OPC<:AbstractGrayOpticalProps,DA}

    if OPC == GrayOneScalar
        return nothing
    else
        ngpt = 1
        src_sfc = DA{FT}(undef, ncol)
        sfc_albedo = DA{FT}(undef, ncol)
        src_up = DA{FT}(undef, nlay, ncol)
        src_dn = DA{FT}(undef, nlay, ncol)
        Rdif = DA{FT}(undef, nlay, ncol)
        Tdif = DA{FT}(undef, nlay, ncol)
        Rdir = DA{FT}(undef, nlay, ncol)
        Tdir = DA{FT}(undef, nlay, ncol)
        Tnoscat = DA{FT}(undef, nlay, ncol)
        return GraySourceSW2Str{FT,typeof(src_sfc),typeof(src_up)}(
            src_sfc,
            sfc_albedo,
            src_up,
            src_dn,
            Rdif,
            Tdif,
            Rdir,
            Tdir,
            Tnoscat,
        )
    end
end

end
