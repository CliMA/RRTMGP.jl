function rte_sw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt, igpt = 1, 1
    cos_zenith = bcs_sw.cos_zenith
    FT = eltype(cos_zenith)
    solar_frac = FT(1)
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            if cos_zenith[gcol] > 0
                compute_optical_props!(op, as, gcol)
                rte_sw_noscat!(flux_sw, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
            else
                set_flux_to_zero!(flux_sw, gcol)
            end
        end
    end
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    flux_up_sw = flux_sw.flux_up
    flux_dn_sw = flux_sw.flux_dn
    flux_net_sw = flux_sw.flux_net
    (; flux_up, flux_dn) = flux
    cos_zenith = bcs_sw.cos_zenith
    @inbounds begin
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                if cos_zenith[gcol] > 0
                    compute_optical_props!(op, as, gcol, igpt, lookup_sw, nothing)
                    solar_frac = lookup_sw.solar_src_scaled[igpt]
                    rte_sw_noscat!(flux, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
                    if igpt == 1
                        map!(x -> x, view(flux_up_sw, :, gcol), view(flux_up, :, gcol))
                        map!(x -> x, view(flux_dn_sw, :, gcol), view(flux_dn, :, gcol))
                    else
                        for ilev in 1:nlev
                            flux_up_sw[ilev, gcol] += flux_up[ilev, gcol]
                            flux_dn_sw[ilev, gcol] += flux_dn[ilev, gcol]
                        end
                    end
                else
                    set_flux_to_zero!(flux_sw, gcol)
                end
            end
        end
        ClimaComms.@threaded device for gcol in 1:ncol
            if cos_zenith[gcol] > 0
                for ilev in 1:nlev
                    flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
                end
            end
        end
    end
    return nothing
end

"""
    rte_sw_noscat!(
        flux::FluxSW,
        op::OneScalar,
        bcs_sw::SwBCs,
        solar_frac::AbstractFloat,
        gcol,
        nlev,
    )

No-scattering solver for the shortwave problem.
(Extinction-only i.e. solar direct beam)
"""
function rte_sw_noscat!(
    flux::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    igpt::Int,
    n_gpt::Int,
    solar_frac::AbstractFloat,
    gcol::Int,
    nlev::Int,
)
    (; toa_flux, cos_zenith) = bcs_sw
    τ = op.τ
    (; flux_dn_dir, flux_net) = flux
    # downward propagation
    @inbounds flux_dn_dir[nlev, gcol] = toa_flux[gcol] * solar_frac * cos_zenith[gcol]
    ilev = nlev - 1
    @inbounds while ilev ≥ 1
        flux_dn_dir[ilev, gcol] = flux_dn_dir[ilev + 1, gcol] * exp(-τ[ilev, gcol] / cos_zenith[gcol])
        flux_net[ilev, gcol] = -flux_dn_dir[ilev, gcol]
        ilev -= 1
    end
end
