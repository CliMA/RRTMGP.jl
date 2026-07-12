#####
##### Device-agnostic pieces shared by the CPU and CUDA drivers of the four
##### spectral solvers, so the per-(g-point, column) orchestration exists
##### exactly once.
#####

# Aerosol mask for one column (no-op without aerosols).
@inline _compute_aero_mask!(::Nothing, gcol) = nothing
@inline _compute_aero_mask!(aerosol_state::AerosolState, gcol) =
    Optics.compute_aero_mask!(
        view(aerosol_state.aero_mask, :, gcol),
        view(aerosol_state.aero_mass, :, :, gcol),
    )

# McICA cloud mask for one column and g-point; returns whether any layer is
# cloudy (for the cloud-cover diagnostic). `S` selects `:mask_lw`/`:mask_sw`.
@inline _build_cloud_mask!(::Nothing, ::Val, gcol) = false
@inline function _build_cloud_mask!(
    cloud_state::CloudState,
    ::Val{S},
    gcol,
) where {S}
    mask = view(getproperty(cloud_state, S), :, gcol)
    Optics.build_cloud_mask!(
        mask,
        view(cloud_state.cld_frac, :, gcol),
        cloud_state.mask_type,
    )
    return any(mask)
end

# Accumulate one g-point's fluxes into the broadband accumulator for one
# column: copy on the first g-point, add afterwards. `flux_net` is not
# accumulated — the drivers recompute it from the broadband up/dn once, after
# the g-point loop. For shortwave, the direct beam is accumulated at every
# level, so the broadband `flux_dn_dir` is a true per-level profile.
@inline function _accumulate_fluxes!(
    flux_acc::FluxLW,
    flux::FluxLW,
    gcol,
    nlev,
    igpt,
)
    acc_up, acc_dn = flux_acc.flux_up, flux_acc.flux_dn
    up, dn = flux.flux_up, flux.flux_dn
    if igpt == 1
        @inbounds for ilev in 1:nlev
            acc_up[gcol, ilev] = up[gcol, ilev]
            acc_dn[gcol, ilev] = dn[gcol, ilev]
        end
    else
        @inbounds for ilev in 1:nlev
            acc_up[gcol, ilev] += up[gcol, ilev]
            acc_dn[gcol, ilev] += dn[gcol, ilev]
        end
    end
    return nothing
end

@inline function _accumulate_fluxes!(
    flux_acc::FluxSW,
    flux::FluxSW,
    gcol,
    nlev,
    igpt,
)
    acc_up, acc_dn, acc_dir =
        flux_acc.flux_up, flux_acc.flux_dn, flux_acc.flux_dn_dir
    up, dn, dir = flux.flux_up, flux.flux_dn, flux.flux_dn_dir
    if igpt == 1
        @inbounds for ilev in 1:nlev
            acc_up[gcol, ilev] = up[gcol, ilev]
            acc_dn[gcol, ilev] = dn[gcol, ilev]
            acc_dir[gcol, ilev] = dir[gcol, ilev]
        end
    else
        @inbounds for ilev in 1:nlev
            acc_up[gcol, ilev] += up[gcol, ilev]
            acc_dn[gcol, ilev] += dn[gcol, ilev]
            acc_dir[gcol, ilev] += dir[gcol, ilev]
        end
    end
    return nothing
end
