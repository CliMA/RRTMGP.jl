module GraySources

using DocStringExtensions
using ..Device: array_type, array_device
using ..GrayOptics: GrayOneScalar, GrayTwoStream, AbstractGrayOpticalProps

using Adapt

export AbstractGraySource,
    GraySourceLWNoScat, GraySourceLW2Str, source_func_longwave_gray_atmos

abstract type AbstractGraySource{FT} end

"""
    GraySourceLWNoScat{
    FT<:AbstractFloat,
    FTA2D<:AbstractArray{FT,2},
} <: AbstractGraySource{FT}

Longwave sources: computed at layer center, layer edges, 
and at the surface for no scattering calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct GraySourceLWNoScat{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractGraySource{FT}
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol, ngpt)`"
    lay_source::FTA2D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA2D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA2D
    "Surface source"
    sfc_source::FTA2D
    "temporary storage array used in no scattering GPU calculations"
    source_up::FTA2D
    "temporary storage array used in no scattering GPU calculations"
    source_dn::FTA2D
end

function Adapt.adapt_structure(to, x::GraySourceLWNoScat)
    FT = eltype(x.lay_source)
    FTA2D = typeof(adapt(to, x.sfc_source))
    GraySourceLWNoScat{FT,FTA2D}(
        adapt(to, x.lay_source),
        adapt(to, x.lev_source_inc),
        adapt(to, x.lev_source_dec),
        adapt(to, x.sfc_source),
        adapt(to, x.source_up),
        adapt(to, x.source_dn),
    )
end

"""
    GraySourceLW2Str{
    FT<:AbstractFloat,
    FTA2D<:AbstractArray{FT,2},
} <: AbstractGraySource{FT}

Longwave sources: computed at layer center, layer edges, 
and at the surface for 2-stream calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct GraySourceLW2Str{FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}} <:
       AbstractGraySource{FT}
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol, ngpt)`"
    lay_source::FTA2D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA2D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA2D
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

function Adapt.adapt_structure(to, x::GraySourceLW2Str)
    FT = eltype(x.lay_source)
    FTA2D = typeof(adapt(to, x.sfc_source))
    GraySourceLW2Str{FT,FTA2D}(
        adapt(to, x.lay_source),
        adapt(to, x.lev_source_inc),
        adapt(to, x.lev_source_dec),
        adapt(to, x.sfc_source),
        adapt(to, x.lev_source),
        adapt(to, x.src_up),
        adapt(to, x.src_dn),
        adapt(to, x.Rdif),
        adapt(to, x.Tdif),
        adapt(to, x.albedo),
        adapt(to, x.src),
    )
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
    lay_source = DA{FT}(undef, nlay, ncol)
    lev_source_inc = DA{FT}(undef, nlay, ncol)
    lev_source_dec = DA{FT}(undef, nlay, ncol)
    sfc_source = DA{FT}(undef, ncol, ngpt)

    if OPC == GrayOneScalar
        source_up = DA{FT}(undef, nlay, ncol)
        source_dn = DA{FT}(undef, nlay, ncol)
        return GraySourceLWNoScat{eltype(lay_source),DA{FT,2}}(
            lay_source,
            lev_source_inc,
            lev_source_dec,
            sfc_source,
            source_up,
            source_dn,
        )
    else
        lev_source = DA{FT}(undef, nlay + 1, ncol)
        src_up = DA{FT}(undef, nlay, ncol)
        src_dn = DA{FT}(undef, nlay, ncol)
        Rdif = DA{FT}(undef, nlay, ncol)
        Tdif = DA{FT}(undef, nlay, ncol)
        albedo = DA{FT}(undef, nlay + 1, ncol)
        src = DA{FT}(undef, nlay + 1, ncol)
        return GraySourceLW2Str{FT,DA{FT,2}}(
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
