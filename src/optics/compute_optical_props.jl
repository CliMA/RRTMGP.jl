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

Compute optical properties for the longwave problem.
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
        τ = view(op.τ, gcol, :)

        lev_src_inc_prev = zero(t_sfc)
        lev_src_dec_prev = zero(t_sfc)
        t_lev_dec = t_lev_col[1]

        for glay in 1:nlay
            col_dry, p_lay, t_lay = as_layerdata[1, glay],
            as_layerdata[2, glay],
            as_layerdata[3, glay]
            # compute gas optics
            τ[glay], _, _, planckfrac = compute_gas_optics(
                lkp,
                vmr,
                col_dry,
                igpt,
                ibnd,
                p_lay,
                t_lay,
                glay,
                gcol,
            )
            # compute longwave source terms
            t_lev_inc = t_lev_col[glay + 1]

            lay_source[gcol, glay] =
                interp1d_equispaced(t_lay, t_planck, totplnk) * planckfrac
            lev_src_inc =
                interp1d_equispaced(t_lev_inc, t_planck, totplnk) * planckfrac
            lev_src_dec =
                interp1d_equispaced(t_lev_dec, t_planck, totplnk) * planckfrac
            if glay == 1
                sfc_source[gcol] =
                    interp1d_equispaced(t_sfc, t_planck, totplnk) * planckfrac
                lev_source[gcol, glay] = lev_src_dec
            else
                lev_source[gcol, glay] = sqrt(lev_src_inc_prev * lev_src_dec)
            end
            lev_src_dec_prev = lev_src_dec
            lev_src_inc_prev = lev_src_inc
            t_lev_dec = t_lev_inc
        end
        lev_source[gcol, nlay + 1] = lev_src_inc_prev
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
        τ = view(op.τ, gcol, :)
        ssa = view(op.ssa, gcol, :)
        g = view(op.g, gcol, :)
    end

    lev_src_inc_prev = zero(t_sfc)
    lev_src_dec_prev = zero(t_sfc)

    @inbounds begin
        t_lev_dec = t_lev_col[1]
        for glay in 1:nlay
            col_dry, p_lay, t_lay = as_layerdata[1, glay],
            as_layerdata[2, glay],
            as_layerdata[3, glay]
            # compute gas optics
            τ[glay], ssa[glay], g[glay], planckfrac = compute_gas_optics(
                lkp,
                vmr,
                col_dry,
                igpt,
                ibnd,
                p_lay,
                t_lay,
                glay,
                gcol,
            )
            # compute longwave source terms
            t_lev_inc = t_lev_col[glay + 1]

            lev_src_inc =
                interp1d_equispaced(t_lev_inc, t_planck, totplnk) * planckfrac
            lev_src_dec =
                interp1d_equispaced(t_lev_dec, t_planck, totplnk) * planckfrac
            if glay == 1
                sfc_source[gcol] =
                    interp1d_equispaced(t_sfc, t_planck, totplnk) * planckfrac
                lev_source[gcol, glay] = lev_src_dec
            else
                lev_source[gcol, glay] = sqrt(lev_src_inc_prev * lev_src_dec)
            end
            lev_src_dec_prev = lev_src_dec
            lev_src_inc_prev = lev_src_inc
            t_lev_dec = t_lev_inc
        end
        lev_source[gcol, nlay + 1] = lev_src_inc_prev
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
        lkp_cld::Union{LookUpCld, Nothing} = nothing,
        lkp_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
    )

Compute optical properties for the shortwave problem.
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
    @inbounds begin
        ibnd = lkp.band_data.major_gpt2bnd[igpt]
        t_sfc = as.t_sfc[gcol]
        as_layerdata = AtmosphericStates.getview_layerdata(as, gcol)
        τ = view(op.τ, gcol, :)
    end
    @inbounds for glay in 1:nlay
        col_dry, p_lay, t_lay =
            as_layerdata[1, glay], as_layerdata[2, glay], as_layerdata[3, glay]
        # compute gas optics
        τ[glay], _, _ = compute_gas_optics(
            lkp,
            vmr,
            col_dry,
            igpt,
            ibnd,
            p_lay,
            t_lay,
            glay,
            gcol,
        )
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
        τ = view(op.τ, gcol, :)
        ssa = view(op.ssa, gcol, :)
        g = view(op.g, gcol, :)
    end
    @inbounds for glay in 1:nlay
        col_dry, p_lay, t_lay =
            as_layerdata[1, glay], as_layerdata[2, glay], as_layerdata[3, glay]
        # compute gas optics
        τ[glay], ssa[glay], g[glay] = compute_gas_optics(
            lkp,
            vmr,
            col_dry,
            igpt,
            ibnd,
            p_lay,
            t_lay,
            glay,
            gcol,
        )
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
