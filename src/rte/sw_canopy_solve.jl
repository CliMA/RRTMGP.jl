using ..Canopy: CanopyState, fill_canopy_sw_optics!, canopy_sw_2stream_coeffs

"""
    rte_sw_2stream!(op, src_sw, bcs_sw, flux, solar_frac, igpt, n_gpt, ibnd, nlev, gcol, canopy)

Canopy-aware two-stream solver for the shortwave problem.

Identical to the standard `rte_sw_2stream!` except that layers `1:ncanlay` use
`canopy_sw_2stream_coeffs` with physically-computed `beta_dir`, while layers
above the canopy use the standard `sw_2stream_coeffs`.
"""
function rte_sw_2stream!(
    (; τ, ssa, g)::TwoStream,
    (; albedo, src)::SourceSW2Str,
    bcs_sw::SwBCs,
    (; flux_up, flux_dn, flux_dn_dir)::FluxSW,
    solar_frac::FT,
    igpt::Int,
    n_gpt::Int,
    ibnd::Int,
    nlev::Int,
    gcol::Int,
    canopy::CanopyState,
) where {FT}
    nlay = nlev - 1
    ncanlay = canopy.ncanlay
    β_dir_val = FT(canopy.beta_dir[ibnd, gcol])

    @inbounds begin
        toa_flux = bcs_sw.toa_flux[gcol]
        sfc_alb_direct = bcs_sw.sfc_alb_direct[ibnd, gcol]
        μ₀ = bcs_sw.cos_zenith[gcol]
    end
    τ_sum = FT(0)
    for ilay in 1:nlay
        @inbounds τ_sum += τ[ilay, gcol]
    end
    flux_dn_dir_top = toa_flux * solar_frac * μ₀
    flux_dn_dir_bot = flux_dn_dir_top * exp(-τ_sum / max(μ₀, eps(FT)))
    @inbounds flux_dn_dir[1, gcol] = flux_dn_dir_bot
    sfc_source = flux_dn_dir_bot * sfc_alb_direct

    @inbounds flux_dn[nlev, gcol] = FT(0)
    @inbounds surface_albedo = albedo[1, gcol] = bcs_sw.sfc_alb_diffuse[ibnd, gcol]
    @inbounds src[1, gcol] = sfc_source
    τ_cum = τ_sum
    albedo_ilev, src_ilev = surface_albedo, sfc_source
    @inbounds for ilev in 1:nlay
        τ_ilev, ssa_ilev, g_ilev = τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        (Rdir, Tdir, _, Rdif, Tdif) = if ilev ≤ ncanlay
            canopy_sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀, β_dir_val)
        else
            sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀)
        end
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)
        albedo_ilevplus1 = Rdif + Tdif * Tdif * albedo_ilev * denom
        τ_cum -= τ_ilev
        τ_cum = max(τ_cum, FT(0))
        flux_dn_dir_ilevplus1 = flux_dn_dir_top * exp(-τ_cum / max(μ₀, eps(FT)))
        src_up_ilev = Rdir * flux_dn_dir_ilevplus1
        src_dn_ilev = Tdir * flux_dn_dir_ilevplus1
        src_ilevplus1 = src_up_ilev + Tdif * denom * (src_ilev + albedo_ilev * src_dn_ilev)
        albedo[ilev + 1, gcol], src[ilev + 1, gcol] = albedo_ilevplus1, src_ilevplus1
        albedo_ilev = albedo_ilevplus1
        src_ilev = src_ilevplus1
    end
    @inbounds flux_up[nlev, gcol] =
        flux_dn[nlev, gcol] * albedo[nlev, gcol] + src[nlev, gcol]

    @inbounds flux_dn_ilevplus1 = flux_dn[nlev, gcol]
    @inbounds flux_dn[nlev, gcol] += flux_dn_dir_top
    τ_cum = FT(0)
    ilev = nlay
    @inbounds while ilev ≥ 1
        τ_ilev, ssa_ilev, g_ilev = τ[ilev, gcol], ssa[ilev, gcol], g[ilev, gcol]
        albedo_ilev, src_ilev = albedo[ilev, gcol], src[ilev, gcol]
        (_, Tdir, _, Rdif, Tdif) = if ilev ≤ ncanlay
            canopy_sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀, β_dir_val)
        else
            sw_2stream_coeffs(τ_ilev, ssa_ilev, g_ilev, μ₀)
        end
        denom = FT(1) / (FT(1) - Rdif * albedo_ilev)
        src_dn_ilev = Tdir * flux_dn_dir_top * exp(-τ_cum / max(μ₀, eps(FT)))
        τ_cum += τ_ilev
        flux_dn_ilev = (Tdif * flux_dn_ilevplus1 + Rdif * src_ilev + src_dn_ilev) * denom
        flux_up[ilev, gcol] = flux_dn_ilev * albedo_ilev + src_ilev
        flux_dn[ilev, gcol] = flux_dn_ilev + flux_dn_dir_top * exp(-τ_cum / max(μ₀, eps(FT)))
        flux_dn_ilevplus1 = flux_dn_ilev
        ilev -= 1
    end
    return nothing
end

"""
    rte_sw_2stream_canopy_solve!(device, flux, flux_sw, op_atm, op_exp, bcs_sw, src_sw, as, lookup_sw, lookup_sw_cld, lookup_sw_aero, canopy)

CPU dispatch for the coupled atmosphere-canopy shortwave 2-stream solver.

For each g-point and column:
1. Computes atmospheric optical properties via `compute_optical_props!`
2. Fills expanded column (atmosphere + canopy) via `fill_canopy_sw_optics!`
3. Solves using canopy-aware `rte_sw_2stream!` with `beta_dir` dispatch
"""
function rte_sw_2stream_canopy_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op_atm::TwoStream,
    op_exp::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing},
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing},
    canopy::CanopyState,
    rad_output::Union{CanopyRadiationOutput, Nothing} = nothing,
)
    nlay_atm, ncol = AtmosphericStates.get_dims(as)
    ncanlay = canopy.ncanlay
    nlay_total = nlay_atm + ncanlay
    nlev_total = nlay_total + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    @inbounds begin
        (; cloud_state, aerosol_state) = as
        bld_cld_mask = cloud_state isa CloudState
        flux_up_sw = flux_sw.flux_up
        flux_dn_sw = flux_sw.flux_dn
        flux_dn_dir_sw = flux_sw.flux_dn_dir
        flux_net_sw = flux_sw.flux_net
        (; flux_up, flux_dn, flux_dn_dir) = flux
        cos_zenith = bcs_sw.cos_zenith
        FT = eltype(flux_up)
        if aerosol_state isa AerosolState
            ClimaComms.@threaded device for gcol in 1:ncol
                Optics.compute_aero_mask!(
                    view(aerosol_state.aero_mask, :, gcol),
                    view(aerosol_state.aero_mass, :, :, gcol),
                )
            end
        end
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                bld_cld_mask && Optics.build_cloud_mask!(
                    view(cloud_state.mask_sw, :, gcol),
                    view(cloud_state.cld_frac, :, gcol),
                    cloud_state.mask_type,
                )
                # Compute atmospheric optical properties
                compute_optical_props!(op_atm, as, gcol, igpt, lookup_sw, lookup_sw_cld, lookup_sw_aero)
                if cos_zenith[gcol] > 0
                    solar_frac = lookup_sw.solar_src_scaled[igpt]
                    ibnd = lookup_sw.band_data.major_gpt2bnd[igpt]
                    # Fill expanded column (atm + canopy)
                    fill_canopy_sw_optics!(op_exp, op_atm, canopy, ibnd, gcol, nlay_atm)
                    # Solve with canopy-aware coefficient dispatch
                    rte_sw_2stream!(op_exp, src_sw, bcs_sw, flux, solar_frac, igpt, n_gpt, ibnd, nlev_total, gcol, canopy)
                    # Accumulate fluxes
                    if igpt == 1
                        map!(x -> x, view(flux_up_sw, :, gcol), view(flux_up, :, gcol))
                        map!(x -> x, view(flux_dn_sw, :, gcol), view(flux_dn, :, gcol))
                        map!(x -> x, view(flux_dn_dir_sw, :, gcol), view(flux_dn_dir, :, gcol))
                    else
                        for ilev in 1:nlev_total
                            @inbounds flux_up_sw[ilev, gcol] += flux_up[ilev, gcol]
                            @inbounds flux_dn_sw[ilev, gcol] += flux_dn[ilev, gcol]
                        end
                        @inbounds flux_dn_dir_sw[1, gcol] += flux_dn_dir[1, gcol]
                    end
                    # Per-band APAR accumulation for fAPAR/green-APAR diagnostics
                    if !isnothing(rad_output)
                        for ilay in 1:ncanlay
                            @inbounds begin
                                flux_in = flux_up[ilay, gcol] + flux_dn[ilay + 1, gcol]
                                flux_out = flux_up[ilay + 1, gcol] + flux_dn[ilay, gcol]
                                rad_output.apar_by_band[ibnd, ilay, gcol] += flux_in - flux_out
                            end
                        end
                    end
                else
                    set_flux_to_zero!(flux_sw, gcol)
                end
            end
        end

        ClimaComms.@threaded device for gcol in 1:ncol
            if cos_zenith[gcol] > 0
                for ilev in 1:nlev_total
                    flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
                end
            end
        end
    end
    return nothing
end
