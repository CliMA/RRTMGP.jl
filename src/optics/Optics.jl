module Optics

using DocStringExtensions
using CUDA
using Adapt
using Random
import ClimaComms

using ..Vmrs
import ..pow_fast
using ..LookUpTables
using ..AtmosphericStates
using ..Sources
using ..AngularDiscretizations
import ..Parameters as RP

export AbstractOpticalProps, OneScalar, TwoStream, compute_col_gas!, compute_optical_props!

"""
    AbstractOpticalProps

Optical properties for one scalar and two stream calculations.
"""
abstract type AbstractOpticalProps end

"""
    OneScalar{FTA2D,AD} <: AbstractOpticalProps

Single scalar approximation for optical depth, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OneScalar{D, V, AD <: AngularDiscretization} <: AbstractOpticalProps
    "storage for optical thickness"
    layerdata::D
    "view into optical depth"
    τ::V
    "Angular discretization"
    angle_disc::AD
end
Adapt.@adapt_structure OneScalar

function OneScalar(::Type{FT}, ncol::Int, nlay::Int, ::Type{DA}) where {FT <: AbstractFloat, DA}
    layerdata = DA{FT, 3}(undef, 1, nlay, ncol)
    τ = view(layerdata, 1, :, :)
    ad = AngularDiscretization(FT, DA, 1)
    V = typeof(τ)
    return OneScalar{typeof(layerdata), V, typeof(ad)}(layerdata, τ, ad)
end

"""
    TwoStream{FTA2D} <: AbstractOpticalProps

Two stream approximation for optical properties, used in
calculations accounting for extinction and emission

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TwoStream{D, V} <: AbstractOpticalProps
    "storage for optical depth, single scattering albedo and asymmerty parameter"
    layerdata::D
    "view into optical depth"
    τ::V
    "view into single scattering albedo"
    ssa::V
    "view into asymmetry parameter"
    g::V
end
Adapt.@adapt_structure TwoStream

function TwoStream(::Type{FT}, ncol::Int, nlay::Int, ::Type{DA}) where {FT <: AbstractFloat, DA}
    layerdata = DA{FT, 3}(zeros(3, nlay, ncol))
    V = typeof(view(layerdata, 1, :, :))
    return TwoStream{typeof(layerdata), V}(
        layerdata,
        view(layerdata, 1, :, :),
        view(layerdata, 2, :, :),
        view(layerdata, 3, :, :),
    )
end

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
        @cuda always_inline = true threads = (tx) blocks = (bx) compute_col_gas_CUDA!(col_dry, args...)
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

"""
    compute_optical_props!(
        op::AbstractOpticalProps,
        as::AtmosphericState{FT},
        sf::AbstractSourceLW{FT},
        gcol::Int,
        igpt::Int,
        lkp::LookUpLW{FT},
        lkp_cld::Union{LookUpCld,PadeCld,Nothing} = nothing,
    ) where {FT<:AbstractFloat}

Computes optical properties for the longwave problem.
"""
@inline function compute_optical_props!(
    op::OneScalar,
    as::AtmosphericState,
    sf::SourceLWNoScat,
    gcol::Int,
    igpt::Int,
    lkp::LookUpLW,
    lkp_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    nlay = AtmosphericStates.get_nlay(as)
    (; vmr) = as
    (; t_planck) = lkp.planck
    (; lay_source, lev_source_inc, lev_source_dec, sfc_source) = sf
    @inbounds begin
        t_sfc = as.t_sfc[gcol]
        ibnd = lkp.band_data.major_gpt2bnd[igpt]
        totplnk = view(lkp.planck.tot_planck, :, ibnd)
        col_dry_col = view(as.col_dry, :, gcol)
        p_lay_col = view(as.p_lay, :, gcol)
        t_lay_col = view(as.t_lay, :, gcol)
        t_lev_col = view(as.t_lev, :, gcol)
        τ = view(op.τ, :, gcol)

        t_lev_dec = t_lev_col[1]
        for glay in 1:nlay
            col_dry = col_dry_col[glay]
            p_lay = p_lay_col[glay]
            t_lay = t_lay_col[glay]
            # compute gas optics
            τ[glay], _, _, planckfrac = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
            # compute longwave source terms
            t_lev_inc = t_lev_col[glay + 1]

            lay_source[glay, gcol] = interp1d(t_lay, t_planck, totplnk) * planckfrac
            lev_source_inc[glay, gcol] = interp1d(t_lev_inc, t_planck, totplnk) * planckfrac
            lev_source_dec[glay, gcol] = interp1d(t_lev_dec, t_planck, totplnk) * planckfrac
            if glay == 1
                sfc_source[gcol] = interp1d(t_sfc, t_planck, totplnk) * planckfrac
            end
            t_lev_dec = t_lev_inc
        end
    end
    return nothing
end

@inline function compute_optical_props!(
    op::TwoStream,
    as::AtmosphericState,
    sf::SourceLW2Str,
    gcol::Int,
    igpt::Int,
    lkp::LookUpLW,
    lkp_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    nlay = AtmosphericStates.get_nlay(as)
    (; vmr) = as
    (; t_planck) = lkp.planck
    (; lev_source, sfc_source) = sf
    @inbounds begin
        t_sfc = as.t_sfc[gcol]
        ibnd = lkp.band_data.major_gpt2bnd[igpt]
        totplnk = view(lkp.planck.tot_planck, :, ibnd)
        col_dry_col = view(as.col_dry, :, gcol)
        p_lay_col = view(as.p_lay, :, gcol)
        t_lay_col = view(as.t_lay, :, gcol)
        t_lev_col = view(as.t_lev, :, gcol)
        τ = view(op.τ, :, gcol)
        ssa = view(op.ssa, :, gcol)
        g = view(op.g, :, gcol)
    end

    lev_src_inc_prev = zero(t_sfc)
    lev_src_dec_prev = zero(t_sfc)

    @inbounds begin
        t_lev_dec = t_lev_col[1]
        for glay in 1:nlay
            col_dry = col_dry_col[glay]
            p_lay = p_lay_col[glay]
            t_lay = t_lay_col[glay]
            # compute gas optics
            τ[glay], ssa[glay], g[glay], planckfrac =
                compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
            # compute longwave source terms
            t_lev_inc = t_lev_col[glay + 1]

            lev_src_inc = interp1d(t_lev_inc, t_planck, totplnk) * planckfrac
            lev_src_dec = interp1d(t_lev_dec, t_planck, totplnk) * planckfrac
            if glay == 1
                sfc_source[gcol] = interp1d(t_sfc, t_planck, totplnk) * planckfrac
                lev_source[glay, gcol] = lev_src_dec
            else
                lev_source[glay, gcol] = sqrt(lev_src_inc_prev * lev_src_dec)
            end
            lev_src_dec_prev = lev_src_dec
            lev_src_inc_prev = lev_src_inc
            t_lev_dec = t_lev_inc
        end
        lev_source[nlay + 1, gcol] = lev_src_inc_prev
    end
    if !isnothing(lkp_cld) # clouds need TwoStream optics
        cld_r_eff_liq = view(as.cld_r_eff_liq, :, gcol)
        cld_r_eff_ice = view(as.cld_r_eff_ice, :, gcol)
        cld_path_liq = view(as.cld_path_liq, :, gcol)
        cld_path_ice = view(as.cld_path_ice, :, gcol)
        cld_mask = view(as.cld_mask_lw, :, gcol)

        add_cloud_optics_2stream!(
            τ,
            ssa,
            g,
            cld_mask,
            cld_r_eff_liq,
            cld_r_eff_ice,
            cld_path_liq,
            cld_path_ice,
            as.ice_rgh,
            lkp_cld,
            ibnd;
            delta_scaling = false,
        )
    end
    return nothing
end

"""
    compute_optical_props!(
        op::AbstractOpticalProps,
        as::AtmosphericState,
        gcol::Int,
        igpt::Int,
        lkp::LookUpSW,
        lkp_cld::Union{LookUpCld,PadeCld,Nothing} = nothing,
    )

Computes optical properties for the shortwave problem.
"""
@inline function compute_optical_props!(
    op::OneScalar,
    as::AtmosphericState,
    gcol::Int,
    igpt::Int,
    lkp::LookUpSW,
    ::Nothing,
)
    nlay = AtmosphericStates.get_nlay(as)
    (; vmr) = as
    @inbounds ibnd = lkp.band_data.major_gpt2bnd[igpt]
    @inbounds t_sfc = as.t_sfc[gcol]
    col_dry_col = view(as.col_dry, :, gcol)
    p_lay_col = view(as.p_lay, :, gcol)
    t_lay_col = view(as.t_lay, :, gcol)
    τ = view(op.τ, :, gcol)

    @inbounds for glay in 1:nlay
        col_dry = col_dry_col[glay]
        p_lay = p_lay_col[glay]
        t_lay = t_lay_col[glay]
        # compute gas optics
        τ[glay], _, _ = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
    end
    return nothing
end

@inline function compute_optical_props!(
    op::TwoStream,
    as::AtmosphericState,
    gcol::Int,
    igpt::Int,
    lkp::LookUpSW,
    lkp_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    nlay = AtmosphericStates.get_nlay(as)
    (; vmr) = as
    @inbounds ibnd = lkp.band_data.major_gpt2bnd[igpt]
    @inbounds t_sfc = as.t_sfc[gcol]
    col_dry_col = view(as.col_dry, :, gcol)
    p_lay_col = view(as.p_lay, :, gcol)
    t_lay_col = view(as.t_lay, :, gcol)
    τ = view(op.τ, :, gcol)
    ssa = view(op.ssa, :, gcol)
    g = view(op.g, :, gcol)

    @inbounds for glay in 1:nlay
        col_dry = col_dry_col[glay]
        p_lay = p_lay_col[glay]
        t_lay = t_lay_col[glay]
        # compute gas optics
        τ[glay], ssa[glay], g[glay] = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
    end
    if !isnothing(lkp_cld) # clouds need TwoStream optics
        cld_r_eff_liq = view(as.cld_r_eff_liq, :, gcol)
        cld_r_eff_ice = view(as.cld_r_eff_ice, :, gcol)
        cld_path_liq = view(as.cld_path_liq, :, gcol)
        cld_path_ice = view(as.cld_path_ice, :, gcol)
        cld_mask = view(as.cld_mask_sw, :, gcol)

        add_cloud_optics_2stream!(
            τ,
            ssa,
            g,
            cld_mask,
            cld_r_eff_liq,
            cld_r_eff_ice,
            cld_path_liq,
            cld_path_ice,
            as.ice_rgh,
            lkp_cld,
            ibnd;
            delta_scaling = true,
        )
    end
    return nothing
end

include("OpticsUtils.jl")
include("GasOptics.jl")
include("CloudOptics.jl")
include("GrayOpticsKernels.jl")

end
