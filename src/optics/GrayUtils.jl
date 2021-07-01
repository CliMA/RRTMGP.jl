module GrayUtils

using CUDA
using UnPack
using ..Device: array_type, array_device, CUDADevice, CPU
using ..AtmosphericStates
using ..Fluxes

using CLIMAParameters
using CLIMAParameters.Planet: grav, cp_d

export update_profile_lw!, compute_gray_heating_rate!, GPU_minmax

function update_profile_lw!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I},
    flux::AbstractFlux{FT,FTA2D},
    hr_lay::FTA2D,
    flux_grad::FTA2D,
    T_ex_lev::FTA2D,
    Δt::FT,
    nlay::I,
    nlev::I,
    ncol::I,
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
    B<:Bool,
}
    # updating t_lay and t_lev based on heating rate
    @unpack p_lay, p_lev, t_lay, t_lev = as
    @unpack flux_dn, flux_net = flux
    args = (
        p_lay,
        p_lev,
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
    )
    device = array_device(as.p_lev)
    if device === CUDADevice()
        max_threads = 256
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        @cuda threads = (tx) blocks = (bx) update_profile_lw_CUDA!(
            ncol,
            args...,
        )
    else # using julia native multithreading
        Threads.@threads for gcol = 1:ncol
            update_profile_lw_kernel!(args..., gcol)
        end
    end
    #----------------------------------------
end

function update_profile_lw_CUDA!(ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        update_profile_lw_kernel!(args..., gcol)
    end
    return nothing
end

function update_profile_lw_kernel!(
    p_lay::FTA2D,
    p_lev::FTA2D,
    t_lay::FTA2D,
    t_lev::FTA2D,
    hr_lay::FTA2D,
    flux_dn::FTA2D,
    flux_net::FTA2D,
    flux_grad::FTA2D,
    T_ex_lev::FTA2D,
    Δt::FT,
    nlay::Int,
    nlev::Int,
    gcol::Int,
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    # updating t_lay based on heating rate
    for glay = 1:nlay
        t_lay[glay, gcol] += Δt * hr_lay[glay, gcol]
    end
    #--------compute t_lev from t_lay--------------------------------------
    for glev = 2:nlay-1
        t_lev[glev, gcol] =
            FT(1 / 3) * t_lay[glev-1, gcol] + FT(5 / 6) * t_lay[glev, gcol] -
            FT(1 / 6) * t_lay[glev+1, gcol]
    end

    t_lev[nlay, gcol] =
        FT(1 / 3) * t_lay[nlay, gcol] + FT(5 / 6) * t_lay[nlay-1, gcol] -
        FT(1 / 6) * t_lay[nlay-2, gcol]

    t_lev[1, gcol] = FT(2) * t_lay[1, gcol] - t_lev[2, gcol]

    t_lev[nlay+1, gcol] = FT(2) * t_lay[nlay, gcol] - t_lev[nlay, gcol]
    #-----------------------------------------------------------------------
    sbc = FT(Stefan())
    for glev = 1:nlev
        T_ex_lev[glev, gcol] =
            (
                (flux_dn[glev, gcol] + (flux_net[glev, gcol] / FT(2))) / sbc
            )^FT(0.25)
    end

    for glev = 2:nlev
        flux_grad[glev-1, gcol] =
            abs(flux_net[glev, gcol] - flux_net[glev-1, gcol])
    end
    #-----------------------------------------------------------------------
end

function compute_gray_heating_rate!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I},
    flux::AbstractFlux{FT,FTA2D},
    hr_lay::FTA2D,
    param_set::AbstractEarthParameterSet,
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
    B<:Bool,
}
    @unpack ncol, nlay, p_lev = as
    flux_net = flux.flux_net

    cp_d_ = FT(cp_d(param_set))
    grav_ = FT(grav(param_set))
    args = (flux_net, p_lev, hr_lay, grav_, cp_d_, nlay)
    device = array_device(as.p_lev)
    if device === CUDADevice()
        max_threads = 256
        tx = min(ncol, max_threads)
        bx = cld(ncol, tx)
        @cuda threads = (tx) blocks = (bx) compute_gray_heating_rate_CUDA!(
            ncol,
            args...,
        )
    else # Launcing native Julia multithreading Kernel
        Threads.@threads for gcol = 1:ncol
            compute_gray_heating_rate_kernel!(args..., gcol)
        end
    end
    #----------------------------------
    return nothing
end

function compute_gray_heating_rate_CUDA!(ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        compute_gray_heating_rate_kernel!(args..., gcol)
    end
    return nothing
end

function compute_gray_heating_rate_kernel!(
    flux_net::FTA2D,
    p_lev::FTA2D,
    hr_lay::FTA2D,
    grav_::FT,
    cp_d_::FT,
    nlay::Int,
    gcol::Int,
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}
    for ilay = 1:nlay
        hr_lay[ilay, gcol] =
            grav_ * (flux_net[ilay+1, gcol] - flux_net[ilay, gcol]) /
            (p_lev[ilay+1, gcol] - p_lev[ilay, gcol]) / cp_d_
    end
    return nothing
end
end
