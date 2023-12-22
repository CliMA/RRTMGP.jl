module Optics

using DocStringExtensions
using CUDA
using Adapt
using Random
import ClimaComms
#---------------------------------------
using ..Vmrs
import ..pow_fast
using ..LookUpTables
using ..AtmosphericStates
using ..Sources
using ..AngularDiscretizations
import ..Parameters as RP
#---------------------------------------

export AbstractOpticalProps, OneScalar, TwoStream, compute_col_gas!, compute_optical_props!, OneScalarc, TwoStreamc

"""
    AbstractOpticalProps{FT,FTA2D}

Optical properties for one scalar and two stream calculations.
"""
abstract type AbstractOpticalProps{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} end

"""
    OneScalar{FT,FTA1D,FTA2D,AD} <: AbstractOpticalProps{FT,FTA2D}

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OneScalar{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}, AD <: AngularDiscretization} <:
       AbstractOpticalProps{FT, FTA2D}
    "Optical Depth"
    τ::FTA2D
    "Angular discretization"
    angle_disc::AD
end
Adapt.@adapt_structure OneScalar

function OneScalar(::Type{FT}, ncol::Int, nlay::Int, ::Type{DA}) where {FT <: AbstractFloat, DA}
    τ = DA{FT, 2}(undef, nlay, ncol)
    ad = AngularDiscretization(FT, DA, 1)

    return OneScalar{eltype(τ), typeof(τ), typeof(ad)}(τ, ad)
end

struct OneScalarc{FTA1D, AD}
    τ::FTA1D
    angle_disc::AD
end
Adapt.@adapt_structure OneScalarc

OneScalarc(op::OneScalar, gcol::Int) = OneScalarc(view(op.τ, :, gcol), op.angle_disc)

"""
    TwoStream{FT,FTA2D} <: AbstractOpticalProps{FT,FTA2D}

Two stream approximation for optical properties, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStream{FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}} <: AbstractOpticalProps{FT, FTA2D}
    "Optical depth"
    τ::FTA2D
    "Single-scattering albedo"
    ssa::FTA2D
    "Asymmetry parameter"
    g::FTA2D
end
Adapt.@adapt_structure TwoStream

function TwoStream(::Type{FT}, ncol::Int, nlay::Int, ::Type{DA}) where {FT <: AbstractFloat, DA}
    τ = DA{FT, 2}(zeros(nlay, ncol))
    ssa = DA{FT, 2}(zeros(nlay, ncol))
    g = DA{FT, 2}(zeros(nlay, ncol))
    return TwoStream{eltype(τ), typeof(τ)}(τ, ssa, g)
end

struct TwoStreamc{FTA1D}
    τ::FTA1D
    ssa::FTA1D
    g::FTA1D
end
Adapt.@adapt_structure TwoStreamc
TwoStreamc(op::TwoStream, gcol::Int) = TwoStreamc(view(op.τ, :, gcol), view(op.ssa, :, gcol), view(op.g, :, gcol))

"""
    compute_col_gas!(
        context,
        p_lev,
        col_dry,
        param_set,
        vmr_h2o,
        lat,
        max_threads = Int(256),
    )

This function computes the column amounts of dry or moist air.

"""
function compute_col_gas!(
    context,
    p_lev::FTA2D,
    col_dry::FTA2D,
    param_set::RP.ARP,
    vmr_h2o::Union{AbstractArray{FT, 2}, Nothing} = nothing,
    lat::Union{AbstractArray{FT, 1}, Nothing} = nothing,
    max_threads::Int = Int(256),
) where {FT <: AbstractFloat, FTA2D <: AbstractArray{FT, 2}}
    nlay, ncol = size(col_dry)
    mol_m_dry = FT(RP.molmass_dryair(param_set))
    mol_m_h2o = FT(RP.molmass_water(param_set))
    avogadro = FT(RP.avogad(param_set))
    helmert1 = FT(RP.grav(param_set))
    args = (p_lev, mol_m_dry, mol_m_h2o, avogadro, helmert1, vmr_h2o, lat)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        tx = min(nlay * ncol, max_threads)
        bx = cld(nlay * ncol, tx)
        @cuda threads = (tx) blocks = (bx) compute_col_gas_CUDA!(col_dry, args...)
    else
        @inbounds begin
            ClimaComms.@threaded device for icnt in 1:(nlay * ncol)
                gcol = cld(icnt, nlay)
                glay = (icnt % nlay == 0) ? nlay : (icnt % nlay)
                compute_col_gas_kernel!(col_dry, args..., glay, gcol)
            end
        end
    end
    return nothing
end

function compute_col_gas_CUDA!(col_dry, args...)
    glx = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlay, ncol = size(col_dry)
    if glx ≤ nlay * ncol
        glay = (glx % nlay == 0) ? nlay : (glx % nlay)
        gcol = cld(glx, nlay)
        compute_col_gas_kernel!(col_dry, args..., glay, gcol)
    end
    return nothing
end
#-----------------------------------------------------------------------------
"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::AtmosphericState{FT},
        sf::AbstractSourceLW{FT},
        gcol::Int,
        igpt::Int,
        lkp::LookUpLW{FT},
        lkp_cld::Union{LookUpCld,PadeCld,Nothing} = nothing,
    ) where {FT<:AbstractFloat}

Computes optical properties for the longwave problem.
"""
function compute_optical_props!(
    op::AbstractOpticalProps{FT},
    as::AtmosphericState{FT},
    sf::AbstractSourceLW{FT},
    gcol::Int,
    igpt::Int,
    lkp::LookUpLW,
    ::Nothing,
) where {FT <: AbstractFloat}
    (; nlay, vmr) = as
    @inbounds ibnd = lkp.major_gpt2bnd[igpt]
    @inbounds t_sfc = as.t_sfc[gcol]
    is2stream = op isa TwoStream
    planck_args = (lkp.t_planck, lkp.totplnk, ibnd)
    @inbounds for glay in 1:nlay
        col_dry = as.col_dry[glay, gcol]
        p_lay = as.p_lay[glay, gcol]
        t_lay = as.t_lay[glay, gcol]
        # gas optics
        τ, ssa, g = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
        op.τ[glay, gcol] = τ
        if is2stream
            op.ssa[glay, gcol] = ssa
            op.g[glay, gcol] = g
        end
        # compute Planck sources
        p_frac = compute_lw_planck_fraction(lkp, vmr, p_lay, t_lay, igpt, ibnd, glay, gcol)

        sf.lay_source[glay, gcol] = interp1d(t_lay, planck_args...) * p_frac
        sf.lev_source_inc[glay, gcol] = interp1d(as.t_lev[glay + 1, gcol], planck_args...) * p_frac
        sf.lev_source_dec[glay, gcol] = interp1d(as.t_lev[glay, gcol], planck_args...) * p_frac

        if glay == 1
            sf.sfc_source[gcol] = interp1d(t_sfc, planck_args...) * p_frac
        end
    end
    return nothing
end

function compute_optical_props!(
    op::TwoStream,
    as::AtmosphericState,
    sf::AbstractSourceLW,
    gcol::Int,
    igpt::Int,
    lkp::LookUpLW,
    lkp_cld::Union{LookUpCld, PadeCld},
)
    (; nlay, vmr, ice_rgh) = as
    @inbounds ibnd = lkp.major_gpt2bnd[igpt]
    @inbounds t_sfc = as.t_sfc[gcol]
    planck_args = (lkp.t_planck, lkp.totplnk, ibnd)
    op_col = TwoStreamc(op, gcol)
    col_dry_col = view(as.col_dry, :, gcol)
    p_lay_col = view(as.p_lay, :, gcol)
    t_lay_col = view(as.t_lay, :, gcol)
    re_liq_col = view(as.cld_r_eff_liq, :, gcol)
    re_ice_col = view(as.cld_r_eff_ice, :, gcol)
    cld_path_liq_col = view(as.cld_path_liq, :, gcol)
    cld_path_ice_col = view(as.cld_path_ice, :, gcol)
    cld_mask_lw_col = view(as.cld_mask_lw, :, gcol)
    (; τ, ssa, g) = op_col
    @inbounds for glay in 1:nlay
        col_dry = col_dry_col[glay]
        p_lay = p_lay_col[glay]
        t_lay = t_lay_col[glay]
        cld_mask = cld_mask_lw_col[glay]
        # compute gas optics
        τ[glay], ssa[glay], g[glay] = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
        # compute Planck sources
        p_frac = compute_lw_planck_fraction(lkp, vmr, p_lay, t_lay, igpt, ibnd, glay, gcol)

        sf.lay_source[glay, gcol] = interp1d(t_lay, planck_args...) * p_frac
        sf.lev_source_inc[glay, gcol] = interp1d(as.t_lev[glay + 1, gcol], planck_args...) * p_frac
        sf.lev_source_dec[glay, gcol] = interp1d(as.t_lev[glay, gcol], planck_args...) * p_frac

        if glay == 1
            sf.sfc_source[gcol] = interp1d(t_sfc, planck_args...) * p_frac
        end
    end
    @inbounds for glay in 1:nlay
        cld_mask = cld_mask_lw_col[glay]
        if cld_mask
            re_liq = re_liq_col[glay]
            re_ice = re_ice_col[glay]
            cld_path_liq = cld_path_liq_col[glay]
            cld_path_ice = cld_path_ice_col[glay]
            # add cloud optics
            τ_cl, ssa_cl, g_cl = compute_cld_props(lkp_cld, re_liq, re_ice, ice_rgh, cld_path_liq, cld_path_ice, ibnd)
            τ_cl, ssa_cl, g_cl = delta_scale(τ_cl, ssa_cl, g_cl)
            τ[glay], ssa[glay], g[glay] = increment_2stream(τ[glay], ssa[glay], g[glay], τ_cl, ssa_cl, g_cl)
        end
    end
    return nothing
end

"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::AtmosphericState{FT},
        gcol::Int,
        igpt::Int,
        lkp::LookUpSW{FT},
        lkp_cld::Union{LookUpCld,PadeCld,Nothing} = nothing,
    ) where {FT<:AbstractFloat}

Computes optical properties for the shortwave problem.
"""
function compute_optical_props!(
    op::AbstractOpticalProps,
    as::AtmosphericState,
    gcol::Int,
    igpt::Int,
    lkp::LookUpSW,
    ::Nothing,
)
    (; nlay, vmr) = as
    is2stream = op isa TwoStream
    @inbounds ibnd = lkp.major_gpt2bnd[igpt]
    @inbounds for glay in 1:nlay
        col_dry = as.col_dry[glay, gcol]
        p_lay = as.p_lay[glay, gcol]
        t_lay = as.t_lay[glay, gcol]
        τ, ssa, g = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
        op.τ[glay, gcol] = τ
        if is2stream
            op.ssa[glay, gcol] = ssa
            op.g[glay, gcol] = g
        end
    end
    return nothing
end

function compute_optical_props!(
    op::TwoStream,
    as::AtmosphericState,
    gcol::Int,
    igpt::Int,
    lkp::LookUpSW,
    lkp_cld::Union{LookUpCld, PadeCld},
)
    (; nlay, vmr, ice_rgh) = as
    @inbounds ibnd = lkp.major_gpt2bnd[igpt]
    op_col = TwoStreamc(op, gcol)
    col_dry_col = view(as.col_dry, :, gcol)
    p_lay_col = view(as.p_lay, :, gcol)
    t_lay_col = view(as.t_lay, :, gcol)
    re_liq_col = view(as.cld_r_eff_liq, :, gcol)
    re_ice_col = view(as.cld_r_eff_ice, :, gcol)
    cld_path_liq_col = view(as.cld_path_liq, :, gcol)
    cld_path_ice_col = view(as.cld_path_ice, :, gcol)
    cld_mask_sw_col = view(as.cld_mask_sw, :, gcol)
    (; τ, ssa, g) = op_col
    @inbounds for glay in 1:nlay
        col_dry = col_dry_col[glay]
        p_lay = p_lay_col[glay]
        t_lay = t_lay_col[glay]
        cld_mask = cld_mask_sw_col[glay]
        # compute gas optics
        τ[glay], ssa[glay], g[glay] = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
    end
    @inbounds for glay in 1:nlay
        cld_mask = cld_mask_sw_col[glay]
        if cld_mask
            re_liq = re_liq_col[glay]
            re_ice = re_ice_col[glay]
            cld_path_liq = cld_path_liq_col[glay]
            cld_path_ice = cld_path_ice_col[glay]
            # add cloud optics
            τ_cl, ssa_cl, g_cl = compute_cld_props(lkp_cld, re_liq, re_ice, ice_rgh, cld_path_liq, cld_path_ice, ibnd)
            τ_cl, ssa_cl, g_cl = delta_scale(τ_cl, ssa_cl, g_cl)
            τ[glay], ssa[glay], g[glay] = increment_2stream(τ[glay], ssa[glay], g[glay], τ_cl, ssa_cl, g_cl)
        end
    end
    return nothing
end
#-----------------------------------------------------------------------------
"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::GrayAtmosphericState{FT},
        sf::AbstractSourceLW{FT},
        gcol::Int,
        igpt::Int = 1,
    ) where {FT<:AbstractFloat}

Computes optical properties for the longwave gray radiation problem.
"""
function compute_optical_props!(
    op::AbstractOpticalProps{FT},
    as::GrayAtmosphericState{FT},
    sf::AbstractSourceLW{FT},
    gcol::Int,
    igpt::Int = 1,
    lkp::Union{AbstractLookUp, Nothing} = nothing,
    lkp_cld::Union{AbstractLookUp, Nothing} = nothing,
) where {FT <: AbstractFloat}
    nlay = as.nlay
    @inbounds for ilay in 1:nlay
        compute_optical_props_kernel!(op, as, ilay, gcol, sf)
    end
    return nothing
end

"""
    compute_optical_props!(
        op::AbstractOpticalProps{FT},
        as::GrayAtmosphericState{FT},
        gcol::Int,
        igpt::Int = 1,
    ) where {FT<:AbstractFloat}

Computes optical properties for the shortwave gray radiation problem.
"""
function compute_optical_props!(
    op::AbstractOpticalProps{FT},
    as::GrayAtmosphericState{FT},
    gcol::Int,
    igpt::Int = 1,
    lkp::Union{AbstractLookUp, Nothing} = nothing,
    lkp_cld::Union{AbstractLookUp, Nothing} = nothing,
) where {FT <: AbstractFloat}
    nlay = as.nlay
    @inbounds for ilay in 1:nlay
        compute_optical_props_kernel!(op, as, ilay, gcol)
    end
    return nothing
end

include("OpticsUtils.jl")
include("GasOptics.jl")
include("CloudOptics.jl")
include("GrayOpticsKernels.jl")

end
