module GrayUtils

using CUDA
import ClimaComms
using ..AtmosphericStates
using ..Fluxes

import ..Parameters as RP

export update_profile_lw!, compute_gray_heating_rate!

"""
    update_profile_lw!(
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
    context,
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
    args = (t_lay, t_lev, hr_lay, flux_dn, flux_net, flux_grad, T_ex_lev, Δt, nlay, nlev)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        max_threads = 256
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        @cuda always_inline = true threads = (tx) blocks = (bx) update_profile_lw_CUDA!(sbc, ncol, args...)
    else
        @inbounds begin
            ClimaComms.@threaded device for gcol in 1:ncol
                update_profile_lw_kernel!(sbc, args..., gcol)
            end
        end
    end
end

function update_profile_lw_CUDA!(sbc, ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        update_profile_lw_kernel!(sbc, args..., gcol)
    end
    return nothing
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
            T_ex_lev[glev, gcol] = sqrt(sqrt((flux_dn[glev, gcol] + (flux_net[glev, gcol] * ft_1by2)) / sbc))
        end

        for glev in 2:nlev
            flux_grad[glev - 1, gcol] = abs(flux_net[glev, gcol] - flux_net[glev - 1, gcol])
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
function compute_gray_heating_rate!(context, hr_lay, p_lev, ncol, nlay, flux_net, cp_d_, grav_)
    args = (hr_lay, flux_net, p_lev, grav_, cp_d_, nlay)
    device = ClimaComms.device(context)
    if device isa ClimaComms.CUDADevice
        max_threads = 256
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        @cuda always_inline = true threads = (tx) blocks = (bx) compute_gray_heating_rate_CUDA!(ncol, args...)
    else
        @inbounds begin
            ClimaComms.@threaded device for gcol in 1:ncol
                compute_gray_heating_rate_kernel!(args..., gcol)
            end
        end
    end
    return nothing
end

function compute_gray_heating_rate_CUDA!(ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        compute_gray_heating_rate_kernel!(args..., gcol)
    end
    return nothing
end

function compute_gray_heating_rate_kernel!(hr_lay, flux_net, p_lev, grav_, cp_d_, nlay, gcol)
    for ilay in 1:nlay
        @inbounds hr_lay[ilay, gcol] =
            grav_ * (flux_net[ilay + 1, gcol] - flux_net[ilay, gcol]) / (p_lev[ilay + 1, gcol] - p_lev[ilay, gcol]) /
            cp_d_
    end
    return nothing
end
end
