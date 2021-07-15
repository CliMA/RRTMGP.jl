module Sources

using Adapt
using DocStringExtensions

export AbstractSourceLW,
    SourceLWNoScat,
    SourceLW2Str,
    SourceSW2Str,
    source_func_longwave,
    source_func_shortwave

"""
    AbstractSourceLW{FT,FTA1D,FTA2D}

Abstract longwave source for no-scattering and two stream longwave solvers.
"""
abstract type AbstractSourceLW{FT,FTA1D,FTA2D} end

"""
    SourceLWNoScat{FT,FTA1D,FTA2D} <: AbstractSourceLW{FT,FTA1D,FTA2D}

Longwave sources: computed at layer center, layer edges, 
and at the surface for no scattering calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceLWNoScat{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
} <: AbstractSourceLW{FT,FTA1D,FTA2D}
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol)`"
    lay_source::FTA2D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay, ncol)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA2D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay, ncol)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA2D
    "Surface source `(ncol)`"
    sfc_source::FTA1D
    "temporary storage array used in no scattering GPU calculations `(nlay, ncol)`"
    src_up::FTA2D
    "temporary storage array used in no scattering GPU calculations `(nlay, ncol)`"
    src_dn::FTA2D
end
Adapt.@adapt_structure SourceLWNoScat

function SourceLWNoScat(
    ::Type{FT},
    ::Type{DA},
    nlay::I,
    ncol::I,
) where {FT<:AbstractFloat,I<:Int,DA}
    FTA1D = DA{FT,1}
    FTA2D = DA{FT,2}
    return SourceLWNoScat{FT,FTA1D,FTA2D}(
        FTA2D(undef, nlay, ncol), # lay_source
        FTA2D(undef, nlay, ncol), # lev_source_inc
        FTA2D(undef, nlay, ncol), # lev_source_dec
        FTA1D(undef, ncol),       # sfc_source
        FTA2D(undef, nlay, ncol), # src_up
        FTA2D(undef, nlay, ncol), # src_dn
    )
end

"""
    SourceLW2Str{FT,FTA1D,FTA2D} <: AbstractSourceLW{FT,FTA1D,FTA2D}

Longwave sources: computed at layer center, layer edges, 
and at the surface for 2-stream calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceLW2Str{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
} <: AbstractSourceLW{FT,FTA1D,FTA2D}
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol)`"
    lay_source::FTA2D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay, ncol)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA2D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay, ncol)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA2D
    "Surface source `(ncol)`"
    sfc_source::FTA1D
    "level source `[W/m2]` `(nlay+1, ncol)`, used in 2 stream calculations"
    lev_source::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    src_up::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    src_dn::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    Rdif::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    Tdif::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    albedo::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    src::FTA2D
end
Adapt.@adapt_structure SourceLW2Str

function SourceLW2Str(
    ::Type{FT},
    ::Type{DA},
    nlay::I,
    ncol::I,
) where {FT<:AbstractFloat,I<:Int,DA}
    FTA1D = DA{FT,1}
    FTA2D = DA{FT,2}
    return SourceLW2Str{FT,FTA1D,FTA2D}(
        FTA2D(undef, nlay, ncol), # lay_source
        FTA2D(undef, nlay, ncol), # lev_source_inc
        FTA2D(undef, nlay, ncol), # lev_source_dec
        FTA1D(undef, ncol), # sfc_source
        FTA2D(undef, nlay + 1, ncol), # lev_source
        FTA2D(undef, nlay, ncol), # src_up
        FTA2D(undef, nlay, ncol), # src_dn
        FTA2D(undef, nlay, ncol), # Rdif
        FTA2D(undef, nlay, ncol), # Tdif
        FTA2D(undef, nlay + 1, ncol), # albedo
        FTA2D(undef, nlay + 1, ncol), # src
    )
end

"""
    source_func_longwave(
        ::Type{FT},
        ncol::I,
        nlay::I,
        OPC::Symbol,
        ::Type{DA},
    ) where {FT<:AbstractFloat,I<:Int,DA}

Initializes the longwave source for one scalar and two stream simulations.
"""
function source_func_longwave(
    ::Type{FT},
    ncol::I,
    nlay::I,
    OPC::Symbol,
    ::Type{DA},
) where {FT<:AbstractFloat,I<:Int,DA}
    lay_source = DA{FT,2}(undef, nlay, ncol)
    lev_source_inc = DA{FT,2}(undef, nlay, ncol)
    lev_source_dec = DA{FT,2}(undef, nlay, ncol)
    sfc_source = DA{FT,1}(undef, ncol)
    src_up = DA{FT,2}(undef, nlay, ncol)
    src_dn = DA{FT,2}(undef, nlay, ncol)

    if OPC === :OneScalar
        return SourceLWNoScat{eltype(lay_source),typeof(sfc_source),DA{FT,2}}(
            lay_source,
            lev_source_inc,
            lev_source_dec,
            sfc_source,
            src_up,
            src_dn,
        )
    else
        lev_source = DA{FT,2}(undef, nlay + 1, ncol)
        Rdif = DA{FT,2}(undef, nlay, ncol)
        Tdif = DA{FT,2}(undef, nlay, ncol)
        albedo = DA{FT,2}(undef, nlay + 1, ncol)
        src = DA{FT,2}(undef, nlay + 1, ncol)
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

"""
    SourceSW2Str{FT,FTA1D,FTA2D}

Shortwave sources: computed at layer center, layer edges, 
and at the surface for 2-stream calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceSW2Str{
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
}
    "surface source `(ncol)`"
    sfc_source::FTA1D
    "albedo `(nlay + 1, ncol)`"
    albedo::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    src_up::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    src_dn::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    Rdif::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    Tdif::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    Rdir::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    Tdir::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay, ncol)`"
    Tnoscat::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    src::FTA2D
end
Adapt.@adapt_structure SourceSW2Str

function SourceSW2Str(
    ::Type{FT},
    ::Type{DA},
    nlay::I,
    ncol::I,
) where {FT<:AbstractFloat,I<:Int,DA}
    FTA1D = DA{FT,1}
    FTA2D = DA{FT,2}

    return SourceSW2Str{FT,FTA1D,FTA2D}(
        FTA1D(undef, ncol),
        FTA2D(undef, nlay + 1, ncol),
        FTA2D(undef, nlay, ncol),
        FTA2D(undef, nlay, ncol),
        FTA2D(undef, nlay, ncol),
        FTA2D(undef, nlay, ncol),
        FTA2D(undef, nlay, ncol),
        FTA2D(undef, nlay, ncol),
        FTA2D(undef, nlay, ncol),
        FTA2D(undef, nlay + 1, ncol),
    )
end

"""
    source_func_shortwave(
        ::Type{FT},
        ncol::I,
        nlay::I,
        opc::Symbol,
        ::Type{DA},
    ) where {FT<:AbstractFloat,I<:Int,DA}

Initializes the shortwave source for one scalar and two stream simulations.
"""
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
