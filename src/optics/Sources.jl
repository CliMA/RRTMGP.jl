module Sources

using Adapt
using DocStringExtensions

import ..Parameters as RP

export AbstractSourceLW, SourceLWNoScat, SourceLW2Str, SourceSW2Str, source_func_longwave, source_func_shortwave

"""
    AbstractSourceLW{FT,FTA1D,FTA2D}

Abstract longwave source for no-scattering and two stream longwave solvers.
"""
abstract type AbstractSourceLW{FT, FTA1D, FTA2D} end

"""
    SourceLWNoScat{FT,FTA1D,FTA2D} <: AbstractSourceLW{FT,FTA1D,FTA2D}

Longwave sources: computed at layer center, layer edges, 
and at the surface for no scattering calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceLWNoScat{
    FT <: AbstractFloat,
    FTA1D <: AbstractArray{FT, 1},
    FTA2D <: AbstractArray{FT, 2},
    PS <: RP.ARP,
} <: AbstractSourceLW{FT, FTA1D, FTA2D}
    "Parameter set"
    param_set::PS
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol)`"
    lay_source::FTA2D
    "Planck source at layer edge in increasing ilay direction `[W/m2]` `(nlay, ncol)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_inc::FTA2D
    "Planck source at layer edge in decreasing ilay direction `[W/m2]` `(nlay, ncol)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source_dec::FTA2D
    "Surface source `[W/m2]` `(ncol)`"
    sfc_source::FTA1D
end
Adapt.@adapt_structure SourceLWNoScat

function SourceLWNoScat(param_set::RP.ARP, ::Type{FT}, ::Type{DA}, nlay::Int, ncol::Int) where {FT <: AbstractFloat, DA}
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}

    lay_source = FTA2D(undef, nlay, ncol) # lay_source
    lev_source_inc = FTA2D(undef, nlay, ncol) # lev_source_inc
    lev_source_dec = FTA2D(undef, nlay, ncol) # lev_source_dec
    sfc_source = FTA1D(undef, ncol)       # sfc_source

    return SourceLWNoScat{eltype(lay_source), typeof(sfc_source), typeof(lay_source), typeof(param_set)}(
        param_set,
        lay_source,
        lev_source_inc,
        lev_source_dec,
        sfc_source,
    )
end

"""
    SourceLW2Str{FT,FTA1D,FTA2D} <: AbstractSourceLW{FT,FTA1D,FTA2D}

Longwave sources: computed at layer center, layer edges, 
and at the surface for 2-stream calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceLW2Str{FT <: AbstractFloat, FTA1D <: AbstractArray{FT, 1}, FTA2D <: AbstractArray{FT, 2}, PS <: RP.ARP} <:
       AbstractSourceLW{FT, FTA1D, FTA2D}
    "Parameter set"
    param_set::PS
    "Surface source `[W/m2]` `(ncol)`"
    sfc_source::FTA1D
    "level source `[W/m2]` `(nlay+1, ncol)`, used in 2 stream calculations"
    lev_source::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    albedo::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    src::FTA2D
end
Adapt.@adapt_structure SourceLW2Str

function SourceLW2Str(param_set::RP.ARP, ::Type{FT}, ::Type{DA}, nlay::Int, ncol::Int) where {FT <: AbstractFloat, DA}
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}

    sfc_source = FTA1D(undef, ncol) # sfc_source
    lev_source = FTA2D(undef, nlay + 1, ncol) # lev_source
    albedo = FTA2D(undef, nlay + 1, ncol) # albedo
    src = FTA2D(undef, nlay + 1, ncol) # src

    return SourceLW2Str{eltype(lev_source), typeof(sfc_source), typeof(lev_source), typeof(param_set)}(
        param_set,
        sfc_source,
        lev_source,
        albedo,
        src,
    )
end

"""
    source_func_longwave(
        param_set,
        ::Type{FT},
        ncol::Int,
        nlay::Int,
        OPC::Symbol,
        ::Type{DA},
    ) where {FT<:AbstractFloat,DA}

Initializes the longwave source for one scalar and two stream simulations.
"""
function source_func_longwave(
    param_set::RP.ARP,
    ::Type{FT},
    ncol::Int,
    nlay::Int,
    OPC::Symbol,
    ::Type{DA},
) where {FT <: AbstractFloat, DA}
    sfc_source = DA{FT, 1}(undef, ncol)

    if OPC === :OneScalar
        lay_source = DA{FT, 2}(undef, nlay, ncol)
        lev_source_inc = DA{FT, 2}(undef, nlay, ncol)
        lev_source_dec = DA{FT, 2}(undef, nlay, ncol)
        return SourceLWNoScat{eltype(lay_source), typeof(sfc_source), typeof(lay_source), typeof(param_set)}(
            param_set,
            lay_source,
            lev_source_inc,
            lev_source_dec,
            sfc_source,
        )
    else
        PS = typeof(param_set)
        lev_source = DA{FT, 2}(undef, nlay + 1, ncol)
        albedo = DA{FT, 2}(undef, nlay + 1, ncol)
        src = DA{FT, 2}(undef, nlay + 1, ncol)
        return SourceLW2Str{FT, typeof(sfc_source), typeof(lev_source), PS}(
            param_set,
            sfc_source,
            lev_source,
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
struct SourceSW2Str{FT <: AbstractFloat, FTA1D <: AbstractArray{FT, 1}, FTA2D <: AbstractArray{FT, 2}}
    "surface source `(ncol)`"
    sfc_source::FTA1D
    "albedo `(nlay + 1, ncol)`"
    albedo::FTA2D
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    src::FTA2D
end
Adapt.@adapt_structure SourceSW2Str

function SourceSW2Str(::Type{FT}, ::Type{DA}, nlay::Int, ncol::Int) where {FT <: AbstractFloat, DA}
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}

    sfc_source = FTA1D(undef, ncol)
    albedo = FTA2D(undef, nlay + 1, ncol)
    src = FTA2D(undef, nlay + 1, ncol)

    return SourceSW2Str{eltype(sfc_source), typeof(sfc_source), typeof(albedo)}(sfc_source, albedo, src)
end

"""
    source_func_shortwave(
        ::Type{FT},
        ncol::Int,
        nlay::Int,
        opc::Symbol,
        ::Type{DA},
    ) where {FT<:AbstractFloat,DA}

Initializes the shortwave source for one scalar and two stream simulations.
"""
function source_func_shortwave(
    ::Type{FT},
    ncol::Int,
    nlay::Int,
    opc::Symbol,
    ::Type{DA},
) where {FT <: AbstractFloat, DA}
    if opc == :OneScalar
        return nothing
    else
        sfc_source = DA{FT}(undef, ncol)
        albedo = DA{FT}(undef, nlay + 1, ncol)
        src = DA{FT}(undef, nlay + 1, ncol)
        return SourceSW2Str{FT, typeof(sfc_source), typeof(src)}(sfc_source, albedo, src)
    end
end

end
