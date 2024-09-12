module Sources

using Adapt
using DocStringExtensions

import ..Parameters as RP

export AbstractSourceLW, SourceLWNoScat, SourceLW2Str, SourceSW2Str, source_func_longwave, source_func_shortwave

"""
    AbstractSourceLW

Abstract longwave source for no-scattering and two stream longwave solvers.
"""
abstract type AbstractSourceLW end

"""
    SourceLWNoScat{FT,FTA1D,FTA2D} <: AbstractSourceLW

Longwave sources: computed at layer center, layer edges, 
and at the surface for no scattering calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceLWNoScat{S, D, PS} <: AbstractSourceLW
    "Parameter set"
    param_set::PS
    "Surface source `[W/m2]` `(ncol)`"
    sfc_source::S
    "Planck source at layer average temperature `[W/m2]` `(ncol, nlay)`"
    lay_source::D
    "Planck level source at layer edges `[W/m2]` `(ncol,  nlay+1)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source::D
end
Adapt.@adapt_structure SourceLWNoScat

function SourceLWNoScat(param_set::RP.ARP, ::Type{FT}, ::Type{DA}, nlay::Int, ncol::Int) where {FT <: AbstractFloat, DA}
    sfc_source = DA{FT, 1}(undef, ncol)
    lay_source = DA{FT, 2}(undef, ncol, nlay)
    lev_source = DA{FT, 2}(undef, ncol, nlay + 1)

    return SourceLWNoScat{typeof(sfc_source), typeof(lay_source), typeof(param_set)}(
        param_set,
        sfc_source,
        lay_source,
        lev_source,
    )
end

"""
    SourceLW2Str{S, D, V, PS} <: AbstractSourceLW

Longwave sources: computed at layer center, layer edges, 
and at the surface for 2-stream calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceLW2Str{S, D, V, PS} <: AbstractSourceLW
    "Parameter set"
    param_set::PS
    "Surface source `[W/m2]` `(ncol)`"
    sfc_source::S
    "storage for level source, albedo and src `(3, ncol, nlay+1)`"
    leveldata::D
    "level source `[W/m2]` `(ncol, nlay+1)`, used in 2 stream calculations"
    lev_source::V
    "temporary storage array, used in 2 stream calculations `(ncol, nlay + 1)`"
    albedo::V
    "temporary storage array, used in 2 stream calculations `(ncol, nlay + 1)`"
    src::V
end
Adapt.@adapt_structure SourceLW2Str

function SourceLW2Str(param_set::RP.ARP, ::Type{FT}, ::Type{DA}, nlay::Int, ncol::Int) where {FT <: AbstractFloat, DA}
    sfc_source = DA{FT, 1}(undef, ncol) # sfc_source
    leveldata = DA{FT, 3}(undef, 3, ncol, nlay + 1)
    lev_source = view(leveldata, 1, :, :) # lev_source
    albedo = view(leveldata, 2, :, :) # albedo
    src = view(leveldata, 3, :, :) # src

    return SourceLW2Str{typeof(sfc_source), typeof(leveldata), typeof(lev_source), typeof(param_set)}(
        param_set,
        sfc_source,
        leveldata,
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
    ) where {FT, DA}

Initializes the longwave source for one scalar and two stream simulations.
"""
source_func_longwave(param_set::RP.ARP, ::Type{FT}, ncol::Int, nlay::Int, OPC::Symbol, ::Type{DA}) where {FT, DA} =
    (OPC === :OneScalar) ? SourceLWNoScat(param_set, FT, DA, nlay, ncol) : SourceLW2Str(param_set, FT, DA, nlay, ncol)

"""
    SourceSW2Str{S,D,V}

Shortwave sources: computed at layer center, layer edges, 
and at the surface for 2-stream calculations

# Fields

$(DocStringExtensions.FIELDS)
"""
struct SourceSW2Str{S, D, V}
    "surface source `(ncol)`"
    sfc_source::S
    "storage for albedo and src `(2, ncol, nlay + 1)`"
    leveldata::D
    "albedo `(ncol, nlay + 1)`"
    albedo::V
    "temporary storage array, used in 2 stream calculations `(ncol, nlay + 1)`"
    src::V
end
Adapt.@adapt_structure SourceSW2Str

function SourceSW2Str(::Type{FT}, ::Type{DA}, nlay::Int, ncol::Int) where {FT <: AbstractFloat, DA}
    sfc_source = DA{FT, 1}(undef, ncol)
    leveldata = DA{FT, 3}(undef, 2, ncol, nlay + 1)
    albedo = view(leveldata, 1, :, :)
    src = view(leveldata, 2, :, :)

    return SourceSW2Str{typeof(sfc_source), typeof(leveldata), typeof(albedo)}(sfc_source, leveldata, albedo, src)
end

"""
    source_func_shortwave(
        ::Type{FT},
        ncol::Int,
        nlay::Int,
        opc::Symbol,
        ::Type{DA},
    ) where {FT,DA}

Initializes the shortwave source for one scalar and two stream simulations.
"""
source_func_shortwave(::Type{FT}, ncol::Int, nlay::Int, opc::Symbol, ::Type{DA}) where {FT, DA} =
    (opc == :OneScalar) ? nothing : SourceSW2Str(FT, DA, nlay, ncol)

end
