module Sources

using Adapt
using DocStringExtensions

using ClimaComms
import ..Parameters as RP
import ..RRTMGPGridParams

export AbstractSourceLW, SourceLWNoScat, SourceLW2Str, SourceSW2Str

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
    "Planck source at layer average temperature `[W/m2]` `(nlay, ncol)`"
    lay_source::D
    "Planck level source at layer edges `[W/m2]` `(nlay+1, ncol)`, includes spectral weighting that accounts for state-dependent frequency to g-space mapping"
    lev_source::D
end
Adapt.@adapt_structure SourceLWNoScat

function SourceLWNoScat(grid_params::RRTMGPGridParams; params::RP.ARP)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    sfc_source = DA{FT, 1}(undef, ncol)
    lay_source = DA{FT, 2}(undef, nlay, ncol)
    lev_source = DA{FT, 2}(undef, nlay + 1, ncol)

    return SourceLWNoScat{
        typeof(sfc_source),
        typeof(lay_source),
        typeof(params),
    }(
        params,
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
    "storage for level source, albedo and src `(3, nlay+1, ncol)`"
    leveldata::D
    "level source `[W/m2]` `(nlay+1, ncol)`, used in 2 stream calculations"
    lev_source::V
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    albedo::V
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    src::V
end

function Adapt.adapt_structure(to, src::SourceLW2Str)
    sfc_source = Adapt.adapt(to, src.sfc_source)
    leveldata = Adapt.adapt(to, src.leveldata)
    lev_source = view(leveldata, 1, :, :)
    albedo = view(leveldata, 2, :, :)
    src_view = view(leveldata, 3, :, :)
    return SourceLW2Str{
        typeof(sfc_source),
        typeof(leveldata),
        typeof(lev_source),
        typeof(src.param_set),
    }(
        src.param_set,
        sfc_source,
        leveldata,
        lev_source,
        albedo,
        src_view,
    )
end

function SourceLW2Str(grid_params::RRTMGPGridParams; params::RP.ARP)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    sfc_source = DA{FT, 1}(undef, ncol)
    leveldata = DA{FT, 3}(undef, 3, nlay + 1, ncol)
    lev_source = view(leveldata, 1, :, :)
    albedo = view(leveldata, 2, :, :)
    src = view(leveldata, 3, :, :)

    return SourceLW2Str{
        typeof(sfc_source),
        typeof(leveldata),
        typeof(lev_source),
        typeof(params),
    }(
        params,
        sfc_source,
        leveldata,
        lev_source,
        albedo,
        src,
    )
end


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
    "storage for albedo and src `(2, nlay + 1, ncol)`"
    leveldata::D
    "albedo `(nlay + 1, ncol)`"
    albedo::V
    "temporary storage array, used in 2 stream calculations `(nlay + 1, ncol)`"
    src::V
end

function Adapt.adapt_structure(to, src::SourceSW2Str)
    sfc_source = Adapt.adapt(to, src.sfc_source)
    leveldata = Adapt.adapt(to, src.leveldata)
    albedo = view(leveldata, 1, :, :)
    src_view = view(leveldata, 2, :, :)
    return SourceSW2Str{typeof(sfc_source), typeof(leveldata), typeof(albedo)}(
        sfc_source,
        leveldata,
        albedo,
        src_view,
    )
end

function SourceSW2Str(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    sfc_source = DA{FT, 1}(undef, ncol)
    leveldata = DA{FT, 3}(undef, 2, nlay + 1, ncol)
    albedo = view(leveldata, 1, :, :)
    src = view(leveldata, 2, :, :)

    return SourceSW2Str{typeof(sfc_source), typeof(leveldata), typeof(albedo)}(
        sfc_source,
        leveldata,
        albedo,
        src,
    )
end

end
