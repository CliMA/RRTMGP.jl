module Sources

using DocStringExtensions
using ..Device: array_type, array_device
#using ..Optics: OneScalar, TwoStream, AbstractOpticalProps

using Adapt

export AbstractSource,
       SourceLWNoScat, 
       SourceLW2Str, 
       SourceSW2Str,
       source_func_longwave,
       source_func_shortwave

abstract type AbstractSource{FT} end

"""
    SourceLWNoScat{
    FT<:AbstractFloat,
    FTA2D<:AbstractArray{FT,2},
} <: AbstractSource{FT}

Longwave sources: computed at layer center, layer edges, 
and at the surface for no scattering calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceLWNoScat{FT<:AbstractFloat,
                      FTA1D<:AbstractArray{FT,1},
                      FTA2D<:AbstractArray{FT,2}} <: AbstractSource{FT}
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol, ngpt)`"
    lay_source::FTA2D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA2D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA2D
    "Surface source"
    sfc_source::FTA1D
    "temporary storage array used in no scattering GPU calculations"
    source_up::FTA2D
    "temporary storage array used in no scattering GPU calculations"
    source_dn::FTA2D
end
Adapt.@adapt_structure SourceLWNoScat

"""
    SourceLW2Str{
    FT<:AbstractFloat,
    FTA2D<:AbstractArray{FT,2},
} <: AbstractSource{FT}

Longwave sources: computed at layer center, layer edges, 
and at the surface for 2-stream calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceLW2Str{FT<:AbstractFloat,
                    FTA1D<:AbstractArray{FT,1},
                    FTA2D<:AbstractArray{FT,2},
} <: AbstractSource{FT}
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol, ngpt)`"
    lay_source::FTA2D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA2D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay+1, ncol, ngpt)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA2D
    "Surface source"
    sfc_source::FTA1D
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
Adapt.@adapt_structure SourceLW2Str

function source_func_longwave(
    ::Type{FT},
    ncol::I,
    nlay::I,
    ngpt::I,
    OPC::Symbol,
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}
    println("DA = $DA")
    lay_source = DA{FT}(undef, nlay, ncol)
    lev_source_inc = DA{FT}(undef, nlay, ncol)
    lev_source_dec = DA{FT}(undef, nlay, ncol)
    sfc_source = DA{FT}(undef, ncol)

    if OPC == :OneScalar
        source_up = DA{FT}(undef, nlay, ncol)
        source_dn = DA{FT}(undef, nlay, ncol)
        return SourceLWNoScat{eltype(lay_source),typeof(sfc_source),DA{FT,2}}(
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
        return SourceLW2Str{FT,typeof(sfc_source),DA{FT,2}}(
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

struct SourceSW2Str{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
} <: AbstractSource{FT}
    sfc_source::FTA1D
    albedo::FTA2D
    src_up::FTA2D
    src_dn::FTA2D
    Rdif::FTA2D
    Tdif::FTA2D
    Rdir::FTA2D
    Tdir::FTA2D
    Tnoscat::FTA2D
    "temporary storage array, used in 2 stream calculations"
    src::FTA2D    
end
Adapt.@adapt_structure SourceSW2Str

function source_func_shortwave(
    ::Type{FT},
    ncol::I,
    nlay::I,
    opc::Symbol,
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}

    if opc == :OneScalar
        return nothing
    else
        sfc_source = DA{FT}(undef, ncol)
        albedo = DA{FT}(undef, nlay + 1, ncol)
        src_up = DA{FT}(undef, nlay, ncol)
        src_dn = DA{FT}(undef, nlay, ncol)
        Rdif = DA{FT}(undef, nlay, ncol)
        Tdif = DA{FT}(undef, nlay, ncol)
        Rdir = DA{FT}(undef, nlay, ncol)
        Tdir = DA{FT}(undef, nlay, ncol)
        Tnoscat = DA{FT}(undef, nlay, ncol)
        src = DA{FT}(undef, nlay + 1, ncol)
        return SourceSW2Str{FT,typeof(sfc_source),typeof(src_up)}(
            sfc_source,
            albedo,
            src_up,
            src_dn,
            Rdif,
            Tdif,
            Rdir,
            Tdir,
            Tnoscat,
            src,
        )
    end
end


end
