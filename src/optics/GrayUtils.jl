module GrayUtils

import ClimaComms
using ..AtmosphericStates
using ..Fluxes

import ..Parameters as RP

export update_profile_lw!, compute_gray_heating_rate!

"""
    update_profile_lw!(
        ::ClimaComms.AbstractDevice,
        sbc,
        t_lay,
        t_lev,
        flux_dn,
        flux_net,
        hr_lay,
        flux_grad,
        T_ex_lev,
        Δt,
        nlay,
        nlev,
        ncol,
    )

Updates t_lay and t_lev based on heating rate.
"""
function update_profile_lw!(
    device::ClimaComms.AbstractCPUDevice,
    sbc,
    t_lay,
    t_lev,
    flux_dn,
    flux_net,
    hr_lay,
    flux_grad,
    T_ex_lev,
    Δt,
    nlay,
    nlev,
    ncol,
)
    # updating t_lay and t_lev based on heating rate
    # TODO: it looks like we should probably reorder these args and eliminate one method?
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            update_profile_lw_kernel!(
                sbc,
                t_lay,
                t_lev,
                hr_lay,
                flux_dn,
                flux_net,
                flux_grad,
                T_ex_lev,
                Δt,
                nlay,
                nlev,
                gcol,
            )
        end
    end
end

function update_profile_lw_kernel!(
    sbc::FT,
    t_lay,
    t_lev,
    hr_lay,
    flux_dn,
    flux_net,
    flux_grad,
    T_ex_lev,
    Δt,
    nlay,
    nlev,
    gcol,
) where {FT}
    ft_5by6 = FT(5 / 6)
    ft_1by3 = FT(1 / 3)
    ft_1by6 = FT(1 / 6)
    ft_1by2 = FT(0.5)
    @inbounds begin
        # updating t_lay based on heating rate
        for glay in 1:nlay
            t_lay[glay, gcol] += Δt * hr_lay[glay, gcol]
        end
        # compute t_lev from t_lay
        for glev in 2:(nlay - 1)
            t_lev[glev, gcol] =
                ft_1by3 * t_lay[glev - 1, gcol] + ft_5by6 * t_lay[glev, gcol] - ft_1by6 * t_lay[glev + 1, gcol]
        end

        t_lev[nlay, gcol] =
            ft_1by3 * t_lay[nlay, gcol] + ft_5by6 * t_lay[nlay - 1, gcol] - ft_1by6 * t_lay[nlay - 2, gcol]

        t_lev[1, gcol] = FT(2) * t_lay[1, gcol] - t_lev[2, gcol]

        t_lev[nlay + 1, gcol] = FT(2) * t_lay[nlay, gcol] - t_lev[nlay, gcol]

        for glev in 1:nlev
            T_ex_lev[glev, gcol] = sqrt(sqrt((flux_dn[gcol, glev] + (flux_net[gcol, glev] * ft_1by2)) / sbc))
        end

        for glev in 2:nlev
            flux_grad[gcol, glev - 1] = abs(flux_net[gcol, glev] - flux_net[gcol, glev - 1])
        end
    end
end

"""
    compute_gray_heating_rate!(
        context,
        hr_lay,
        p_lev,
        ncol,
        nlay,
        flux_net,
        cp_d_,
        grav_,
    )

Computes the heating rate for the gray radiation simulation.
"""
function compute_gray_heating_rate!(
    device::ClimaComms.AbstractCPUDevice,
    hr_lay,
    p_lev,
    ncol,
    nlay,
    flux_net,
    cp_d_,
    grav_,
)
    args = (hr_lay, flux_net, p_lev, grav_, cp_d_, nlay)
    @inbounds begin
        ClimaComms.@threaded device for gcol in 1:ncol
            compute_gray_heating_rate_kernel!(args..., gcol)
        end
    end
    return nothing
end

function compute_gray_heating_rate_kernel!(hr_lay, flux_net, p_lev, grav_, cp_d_, nlay, gcol)
    for ilay in 1:nlay
        @inbounds hr_lay[ilay, gcol] =
            grav_ * (flux_net[gcol, ilay + 1] - flux_net[gcol, ilay]) / (p_lev[ilay + 1, gcol] - p_lev[ilay, gcol]) /
            cp_d_
    end
    return nothing
end
end
