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
                rte_sw_noscat!(
                    flux_sw,
                    op,
                    bcs_sw,
                    igpt,
                    n_gpt,
                    solar_frac,
                    gcol,
                    nlev,
                )
                compute_net_flux!(flux_sw, gcol, nlev)
            else
                set_flux_to_zero!(flux_sw, gcol, nlev)
            end
        end
    end
    return nothing
end

# Device-agnostic per-(g-point, column) body, shared by the CPU driver below
# and the CUDA kernel in ext/cuda. Callers guard on cos_zenith > 0 (night
# columns are zeroed once, after the g-point loop).
@inline function sw_noscat_gpt_col!(
    igpt,
    gcol,
    flux,
    flux_sw,
    op,
    bcs_sw,
    as,
    state_cache,
    lookup_sw,
    n_gpt,
    nlev,
)
    compute_optical_props!(op, as, state_cache, gcol, igpt, lookup_sw, nothing)
    @inbounds solar_frac = lookup_sw.solar_src_scaled[igpt]
    rte_sw_noscat!(flux, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
    _accumulate_fluxes!(flux_sw, flux, gcol, nlev, igpt)
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.AbstractCPUDevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    cos_zenith = bcs_sw.cos_zenith
    @inbounds begin
        for igpt in 1:n_gpt
            ClimaComms.@threaded device for gcol in 1:ncol
                if cos_zenith[gcol] > 0
                    sw_noscat_gpt_col!(
                        igpt,
                        gcol,
                        flux,
                        flux_sw,
                        op,
                        bcs_sw,
                        as,
                        state_cache,
                        lookup_sw,
                        n_gpt,
                        nlev,
                    )
                end
            end
        end
        ClimaComms.@threaded device for gcol in 1:ncol
            if cos_zenith[gcol] > 0
                compute_net_flux!(flux_sw, gcol, nlev)
            else
                set_flux_to_zero!(flux_sw, gcol, nlev)
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
        igpt::Int,
        n_gpt::Int,
        solar_frac::AbstractFloat,
        gcol,
        nlev,
    )

No-scattering solver for the shortwave problem.
(Extinction-only i.e., solar direct beam)
"""
@inline function rte_sw_noscat!(
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
    (; flux_up, flux_dn, flux_dn_dir) = flux
    FT = eltype(toa_flux)
    # downward propagation
    @inbounds flux_dn_dir[gcol, nlev] =
        toa_flux[gcol] * solar_frac * cos_zenith[gcol]
    @inbounds flux_dn[gcol, nlev] = flux_dn_dir[gcol, nlev]
    @inbounds flux_up[gcol, nlev] = FT(0)
    ilev = nlev - 1
    @inbounds while ilev ≥ 1
        flux_dn_dir[gcol, ilev] =
            flux_dn_dir[gcol, ilev + 1] *
            exp(-τ[gcol, ilev] / max(cos_zenith[gcol], Numerics.μ₀_min(FT)))
        flux_dn[gcol, ilev] = flux_dn_dir[gcol, ilev]
        flux_up[gcol, ilev] = FT(0)
        ilev -= 1
    end
end
