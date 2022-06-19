module GrayUtils

using CUDA
using ..Device: array_type, array_device, CUDADevice, CPU
using ..AtmosphericStates
using ..Fluxes

import ..Parameters as RP

export update_profile_lw!, compute_gray_heating_rate!

"""
    update_profile_lw!(
        param_set,
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
        ::Type{FT},
    )

Updates t_lay and t_lev based on heating rate.
"""
function update_profile_lw!(
    param_set,
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
    ::Type{FT},
) where {FT <: AbstractFloat}
    # updating t_lay and t_lev based on heating rate
    args = (t_lay, t_lev, hr_lay, flux_dn, flux_net, flux_grad, T_ex_lev, Δt, nlay, nlev)
    device = array_device(t_lev)
    if device === CUDADevice()
        max_threads = 256
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        @cuda threads = (tx) blocks = (bx) update_profile_lw_CUDA!(param_set, ncol, args...)
    else
        # using julia native multithreading
        Threads.@threads for gcol in 1:ncol
            update_profile_lw_kernel!(param_set, args..., gcol)
        end
    end
    #----------------------------------------
end

function update_profile_lw_CUDA!(param_set::RP.ARP, ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        update_profile_lw_kernel!(param_set, args..., gcol)
    end
    return nothing
end

function update_profile_lw_kernel!(
    param_set::RP.ARP,
    t_lay,
    t_lev,
    hr_lay,
    flux_dn,
    flux_net,
    flux_grad,
    T_ex_lev,
    Δt::FT,
    nlay,
    nlev,
    gcol,
) where {FT <: AbstractFloat}
    # updating t_lay based on heating rate
    for glay in 1:nlay
        @inbounds t_lay[glay, gcol] += Δt * hr_lay[glay, gcol]
    end
    #--------compute t_lev from t_lay--------------------------------------
    for glev in 2:(nlay - 1)
        @inbounds t_lev[glev, gcol] =
            FT(1 / 3) * t_lay[glev - 1, gcol] + FT(5 / 6) * t_lay[glev, gcol] - FT(1 / 6) * t_lay[glev + 1, gcol]
    end

    @inbounds t_lev[nlay, gcol] =
        FT(1 / 3) * t_lay[nlay, gcol] + FT(5 / 6) * t_lay[nlay - 1, gcol] - FT(1 / 6) * t_lay[nlay - 2, gcol]

    @inbounds t_lev[1, gcol] = FT(2) * t_lay[1, gcol] - t_lev[2, gcol]

    @inbounds t_lev[nlay + 1, gcol] = FT(2) * t_lay[nlay, gcol] - t_lev[nlay, gcol]
    #-----------------------------------------------------------------------
    sbc = FT(RP.Stefan(param_set))
    for glev in 1:nlev
        @inbounds T_ex_lev[glev, gcol] = ((flux_dn[glev, gcol] + (flux_net[glev, gcol] / FT(2))) / sbc)^FT(0.25)
    end

    for glev in 2:nlev
        @inbounds flux_grad[glev - 1, gcol] = abs(flux_net[glev, gcol] - flux_net[glev - 1, gcol])
    end
    #-----------------------------------------------------------------------
end

"""
    compute_gray_heating_rate!(
        hr_lay,
        p_lev,
        ncol,
        nlay,
        flux_net,
        param_set,
        ::Type{FT},
    )

Computes the heating rate for the gray radiation simulation.
"""
function compute_gray_heating_rate!(hr_lay, p_lev, ncol, nlay, flux_net, param_set, ::Type{FT}) where {FT}
    cp_d_ = FT(RP.cp_d(param_set))
    grav_ = FT(RP.grav(param_set))
    args = (hr_lay, flux_net, p_lev, grav_, cp_d_, nlay)
    device = array_device(p_lev)
    if device === CUDADevice()
        max_threads = 256
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        @cuda threads = (tx) blocks = (bx) compute_gray_heating_rate_CUDA!(ncol, args...)
    else
        # Launcing native Julia multithreading Kernel
        Threads.@threads for gcol in 1:ncol
            compute_gray_heating_rate_kernel!(args..., gcol)
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
