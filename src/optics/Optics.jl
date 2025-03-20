module Optics

using DocStringExtensions
using Adapt
using Random
import ClimaComms

using ..Vmrs
import ..pow_fast
using ..LookUpTables
using ..AtmosphericStates
import ..RRTMGPGridParams
using ..Sources
using ..AngularDiscretizations
import ..Parameters as RP

export AbstractOpticalProps, OneScalar, TwoStream, compute_col_gas!, compute_relative_humidity!, compute_optical_props!

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
struct OneScalar{D, V} <: AbstractOpticalProps
    "storage for optical thickness"
    layerdata::D
    "view into optical depth"
    τ::V
end
Adapt.@adapt_structure OneScalar

function OneScalar(::Type{FT}, ncol::Int, nlay::Int, ::Type{DA}) where {FT <: AbstractFloat, DA}
    layerdata = DA{FT, 3}(undef, 1, nlay, ncol)
    τ = view(layerdata, 1, :, :)
    V = typeof(τ)
    @warn "Please use OneScalar with RRTMGPGridParams instead."
    return OneScalar{typeof(layerdata), V}(layerdata, τ)
end

function OneScalar(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
    layerdata = DA{FT, 3}(undef, 1, nlay, ncol)
    τ = view(layerdata, 1, :, :)
    V = typeof(τ)
    return OneScalar{typeof(layerdata), V}(layerdata, τ)
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
    @warn "Please use TwoStream with RRTMGPGridParams instead"
    return TwoStream{typeof(layerdata), V}(
        layerdata,
        view(layerdata, 1, :, :),
        view(layerdata, 2, :, :),
        view(layerdata, 3, :, :),
    )
end

function TwoStream(grid_params::RRTMGPGridParams)
    (; ncol, nlay) = grid_params
    DA = ClimaComms.array_type(grid_params)
    FT = eltype(grid_params)
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
    device::ClimaComms.AbstractCPUDevice,
    p_lev::AbstractArray{FT, 2},
    col_dry::AbstractArray{FT, 2},
    param_set::RP.ARP,
    vmr_h2o::Union{AbstractArray{FT, 2}, Nothing} = nothing,
    lat::Union{AbstractArray{FT, 1}, Nothing} = nothing,
    max_threads::Int = Int(256),
) where {FT}
    nlay, ncol = size(col_dry)
    mol_m_dry = RP.molmass_dryair(param_set)
    mol_m_h2o = RP.molmass_water(param_set)
    avogadro = RP.avogad(param_set)
    helmert1 = RP.grav(param_set)
    args = (p_lev, mol_m_dry, mol_m_h2o, avogadro, helmert1, vmr_h2o, lat)
    @inbounds begin
        ClimaComms.@threaded device for icnt in 1:(nlay * ncol)
            gcol = cld(icnt, nlay)
            glay = (icnt % nlay == 0) ? nlay : (icnt % nlay)
            compute_col_gas_kernel!(col_dry, args..., glay, gcol)
        end
    end
    return nothing
end

"""
    compute_relative_humidity!(
        device::ClimaComms.AbstractCPUDevice,
        rh::AbstractArray{FT, 2},
        p_lay::AbstractArray{FT, 2},
        t_lay::AbstractArray{FT, 2},
        param_set::RP.ARP,
        vmr_h2o::Union{AbstractArray{FT, 2}, Nothing} = nothing,
    ) where {FT}

This function computes the relative humidity.

"""
function compute_relative_humidity!(
    device::ClimaComms.AbstractCPUDevice,
    rh::AbstractArray{FT, 2},
    p_lay::AbstractArray{FT, 2},
    t_lay::AbstractArray{FT, 2},
    param_set::RP.ARP,
    vmr_h2o::AbstractArray{FT, 2},
) where {FT}
    nlay, ncol = size(p_lay)
    # ratio of water to dry air molecular weights
    mwd = RP.molmass_water(param_set) / RP.molmass_dryair(param_set)
    t_ref = FT(273.16) # reference temperature (K)
    q_lay_min = FT(1e-7) # minimum water mass mixing ratio

    args = (rh, p_lay, t_lay, vmr_h2o, mwd, t_ref, q_lay_min)
    @inbounds begin
        ClimaComms.@threaded device for icnt in 1:(nlay * ncol)
            gcol = cld(icnt, nlay)
            glay = (icnt % nlay == 0) ? nlay : (icnt % nlay)
            compute_relative_humidity_kernel!(args..., glay, gcol)
        end
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
        lkp_cld::Union{LookUpCld,Nothing} = nothing,
        lkp_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
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
    lkp_cld::Union{LookUpCld, Nothing} = nothing,
    lkp_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay = AtmosphericStates.get_nlay(as)
    (; vmr) = as
    (; t_planck) = lkp.planck
    (; lay_source, lev_source, sfc_source) = sf
    @inbounds begin
        t_sfc = as.t_sfc[gcol]
        ibnd = lkp.band_data.major_gpt2bnd[igpt]
        totplnk = view(lkp.planck.tot_planck, :, ibnd)
        as_layerdata = AtmosphericStates.getview_layerdata(as, gcol)
        t_lev_col = view(as.t_lev, :, gcol)
        τ = view(op.τ, :, gcol)

        lev_src_inc_prev = zero(t_sfc)
        lev_src_dec_prev = zero(t_sfc)
        t_lev_dec = t_lev_col[1]

        for glay in 1:nlay
            col_dry, p_lay, t_lay = as_layerdata[1, glay], as_layerdata[2, glay], as_layerdata[3, glay]
            # compute gas optics
            τ[glay], _, _, planckfrac = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
            # compute longwave source terms
            t_lev_inc = t_lev_col[glay + 1]

            lay_source[glay, gcol] = interp1d_equispaced(t_lay, t_planck, totplnk) * planckfrac
            lev_src_inc = interp1d_equispaced(t_lev_inc, t_planck, totplnk) * planckfrac
            lev_src_dec = interp1d_equispaced(t_lev_dec, t_planck, totplnk) * planckfrac
            if glay == 1
                sfc_source[gcol] = interp1d_equispaced(t_sfc, t_planck, totplnk) * planckfrac
                lev_source[glay, gcol] = lev_src_dec
            else
                lev_source[glay, gcol] = sqrt(lev_src_inc_prev * lev_src_dec)
            end
            lev_src_dec_prev = lev_src_dec
            lev_src_inc_prev = lev_src_inc
            t_lev_dec = t_lev_inc
        end
        lev_source[nlay + 1, gcol] = lev_src_inc_prev
        if !isnothing(lkp_cld)
            cloud_state = as.cloud_state
            cld_r_eff_liq = view(cloud_state.cld_r_eff_liq, :, gcol)
            cld_r_eff_ice = view(cloud_state.cld_r_eff_ice, :, gcol)
            cld_path_liq = view(cloud_state.cld_path_liq, :, gcol)
            cld_path_ice = view(cloud_state.cld_path_ice, :, gcol)
            cld_mask = view(cloud_state.mask_lw, :, gcol)

            add_cloud_optics_1scalar!(
                τ,
                cld_mask,
                cld_r_eff_liq,
                cld_r_eff_ice,
                cld_path_liq,
                cld_path_ice,
                cloud_state.ice_rgh,
                lkp_cld,
                ibnd;
            )
        end
        if !isnothing(lkp_aero)
            aod_sw_ext = nothing
            aod_sw_sca = nothing
            iband_550nm = nothing
            aero_mask = view(as.aerosol_state.aero_mask, :, gcol)
            aero_size = view(as.aerosol_state.aero_size, :, :, gcol)
            aero_mass = view(as.aerosol_state.aero_mass, :, :, gcol)
            rel_hum = AtmosphericStates.getview_rel_hum(as, gcol)

            add_aerosol_optics_1scalar!(
                τ,
                aod_sw_ext,
                aod_sw_sca,
                aero_mask,
                aero_size,
                aero_mass,
                rel_hum,
                lkp_aero,
                ibnd,
                iband_550nm,
            )
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
    lkp_cld::Union{LookUpCld, Nothing} = nothing,
    lkp_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay = AtmosphericStates.get_nlay(as)
    (; vmr) = as
    (; t_planck) = lkp.planck
    (; lev_source, sfc_source) = sf
    @inbounds begin
        t_sfc = as.t_sfc[gcol]
        ibnd = lkp.band_data.major_gpt2bnd[igpt]
        totplnk = view(lkp.planck.tot_planck, :, ibnd)
        as_layerdata = AtmosphericStates.getview_layerdata(as, gcol)
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
            col_dry, p_lay, t_lay = as_layerdata[1, glay], as_layerdata[2, glay], as_layerdata[3, glay]
            # compute gas optics
            τ[glay], ssa[glay], g[glay], planckfrac =
                compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
            # compute longwave source terms
            t_lev_inc = t_lev_col[glay + 1]

            lev_src_inc = interp1d_equispaced(t_lev_inc, t_planck, totplnk) * planckfrac
            lev_src_dec = interp1d_equispaced(t_lev_dec, t_planck, totplnk) * planckfrac
            if glay == 1
                sfc_source[gcol] = interp1d_equispaced(t_sfc, t_planck, totplnk) * planckfrac
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
        cloud_state = as.cloud_state
        cld_r_eff_liq = view(cloud_state.cld_r_eff_liq, :, gcol)
        cld_r_eff_ice = view(cloud_state.cld_r_eff_ice, :, gcol)
        cld_path_liq = view(cloud_state.cld_path_liq, :, gcol)
        cld_path_ice = view(cloud_state.cld_path_ice, :, gcol)
        cld_mask = view(cloud_state.mask_lw, :, gcol)

        add_cloud_optics_2stream!(
            τ,
            ssa,
            g,
            cld_mask,
            cld_r_eff_liq,
            cld_r_eff_ice,
            cld_path_liq,
            cld_path_ice,
            cloud_state.ice_rgh,
            lkp_cld,
            ibnd;
            delta_scaling = false,
        )
    end
    if !isnothing(lkp_aero)
        aod_sw_ext = nothing
        aod_sw_sca = nothing
        iband_550nm = nothing
        aero_mask = view(as.aerosol_state.aero_mask, :, gcol)
        aero_size = view(as.aerosol_state.aero_size, :, :, gcol)
        aero_mass = view(as.aerosol_state.aero_mass, :, :, gcol)
        rel_hum = AtmosphericStates.getview_rel_hum(as, gcol)

        add_aerosol_optics_2stream!(
            τ,
            ssa,
            g,
            aod_sw_ext,
            aod_sw_sca,
            aero_mask,
            aero_size,
            aero_mass,
            rel_hum,
            lkp_aero,
            ibnd,
            iband_550nm,
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
        lkp_cld::Union{LookUpCld,Nothing} = nothing,
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
    ::Nothing,
)
    nlay = AtmosphericStates.get_nlay(as)
    (; vmr) = as
    @inbounds begin
        ibnd = lkp.band_data.major_gpt2bnd[igpt]
        t_sfc = as.t_sfc[gcol]
        as_layerdata = AtmosphericStates.getview_layerdata(as, gcol)
        τ = view(op.τ, :, gcol)
    end
    @inbounds for glay in 1:nlay
        col_dry, p_lay, t_lay = as_layerdata[1, glay], as_layerdata[2, glay], as_layerdata[3, glay]
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
    lkp_cld::Union{LookUpCld, Nothing} = nothing,
    lkp_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay = AtmosphericStates.get_nlay(as)
    (; vmr) = as
    @inbounds begin
        ibnd = lkp.band_data.major_gpt2bnd[igpt]
        t_sfc = as.t_sfc[gcol]
        as_layerdata = AtmosphericStates.getview_layerdata(as, gcol)
        τ = view(op.τ, :, gcol)
        ssa = view(op.ssa, :, gcol)
        g = view(op.g, :, gcol)
    end
    @inbounds for glay in 1:nlay
        col_dry, p_lay, t_lay = as_layerdata[1, glay], as_layerdata[2, glay], as_layerdata[3, glay]
        # compute gas optics
        τ[glay], ssa[glay], g[glay] = compute_gas_optics(lkp, vmr, col_dry, igpt, ibnd, p_lay, t_lay, glay, gcol)
    end
    if !isnothing(lkp_cld) # clouds need TwoStream optics
        cloud_state = as.cloud_state
        cld_r_eff_liq = view(cloud_state.cld_r_eff_liq, :, gcol)
        cld_r_eff_ice = view(cloud_state.cld_r_eff_ice, :, gcol)
        cld_path_liq = view(cloud_state.cld_path_liq, :, gcol)
        cld_path_ice = view(cloud_state.cld_path_ice, :, gcol)
        cld_mask = view(cloud_state.mask_sw, :, gcol)

        add_cloud_optics_2stream!(
            τ,
            ssa,
            g,
            cld_mask,
            cld_r_eff_liq,
            cld_r_eff_ice,
            cld_path_liq,
            cld_path_ice,
            cloud_state.ice_rgh,
            lkp_cld,
            ibnd;
            delta_scaling = true,
        )
    end
    if !isnothing(lkp_aero)
        (; iband_550nm) = lkp_aero
        aod_sw_ext = view(as.aerosol_state.aod_sw_ext, gcol)
        aod_sw_sca = view(as.aerosol_state.aod_sw_sca, gcol)
        aero_mask = view(as.aerosol_state.aero_mask, :, gcol)
        aero_size = view(as.aerosol_state.aero_size, :, :, gcol)
        aero_mass = view(as.aerosol_state.aero_mass, :, :, gcol)
        rel_hum = AtmosphericStates.getview_rel_hum(as, gcol)

        add_aerosol_optics_2stream!(
            τ,
            ssa,
            g,
            aod_sw_ext,
            aod_sw_sca,
            aero_mask,
            aero_size,
            aero_mass,
            rel_hum,
            lkp_aero,
            ibnd,
            iband_550nm,
            delta_scaling = true,
        )
    end
    return nothing
end

include("optics_utils.jl")
include("gas_optics.jl")
include("cloud_optics.jl")
include("aerosol_optics.jl")
include("gray_optics_kernels.jl")

end
