module GrayUtils

using KernelAbstractions
using CUDA
using ..Device: array_type, array_device
using ..GrayAtmosphericStates
using ..GrayFluxes

using CLIMAParameters
using CLIMAParameters.Planet: grav, cp_d

export update_profile_lw!, compute_gray_heating_rate!, GPU_minmax

function update_profile_lw!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I},
    flux::GrayFlux{FT,FTA2D},
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
    max_threads = 256
    thr_y = min(32, ncol)
    thr_x = min(Int(floor(FT(max_threads / thr_y))), nlev)

    device = array_device(as.t_lay)
    workgroup = (thr_x, thr_y)
    kernel = update_profile_lw_kernel!(device, workgroup)
    event = kernel(
        as.p_lay,
        as.p_lev,
        as.t_lay,
        as.t_lev,
        hr_lay,
        flux.flux_dn,
        flux.flux_net,
        flux_grad,
        T_ex_lev,
        Δt,
        nlay,
        nlev,
        ncol,
        ndrange = (nlev, ncol),
    )
    wait(event)
    #----------------------------------------
end


@kernel function update_profile_lw_kernel!(
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
    ncol::Int,
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2}}

    glev, gcol = @index(Global, NTuple)  # global col & lay ids
    glay = glev - 1

    # updating t_lay based on heating rate
    if glay > 0
        t_lay[glay, gcol] += Δt * hr_lay[glay, gcol]
    end
    @synchronize
    #--------compute t_lev from t_lay--------------------------------------
    if glev > 1 && glev < nlay
        t_lev[glev, gcol] =
            FT(1 / 3) * t_lay[glev-1, gcol] + FT(5 / 6) * t_lay[glev, gcol] -
            FT(1 / 6) * t_lay[glev+1, gcol]
    end

    if glev == nlay
        t_lev[nlay, gcol] =
            FT(1 / 3) * t_lay[nlay, gcol] + FT(5 / 6) * t_lay[nlay-1, gcol] -
            FT(1 / 6) * t_lay[nlay-2, gcol]
    end

    @synchronize

    if glev == 1
        t_lev[1, gcol] = FT(2) * t_lay[1, gcol] - t_lev[2, gcol]
    end

    if glev == nlev
        t_lev[nlay+1, gcol] = FT(2) * t_lay[nlay, gcol] - t_lev[nlay, gcol]
    end
    @synchronize
    #-----------------------------------------------------------------------
    sbc = FT(Stefan())

    T_ex_lev[glev, gcol] =
        ((flux_dn[glev, gcol] + (flux_net[glev, gcol] / FT(2))) / sbc)^FT(0.25)

    if glev > 1
        flux_grad[glev-1, gcol] =
            abs(flux_net[glev, gcol] - flux_net[glev-1, gcol])
    end
    @synchronize
    #-----------------------------------------------------------------------
end

function compute_gray_heating_rate!(
    as::GrayAtmosphericState{FT,FTA1D,FTA2D,I},
    flux::GrayFlux{FT,FTA2D},
    hr_lay::FTA2D,
    param_set::AbstractEarthParameterSet,
) where {
    FT<:AbstractFloat,
    FTA1D<:AbstractArray{FT,1},
    FTA2D<:AbstractArray{FT,2},
    I<:Int,
    B<:Bool,
}
    ncol = as.ncol
    nlay = as.nlay
    p_lev = as.p_lev
    flux_net = flux.flux_net

    cp_d_ = FT(cp_d(param_set))
    grav_ = FT(grav(param_set))
    #----Launcing KA Kernel---------------------------
    thr_x = min(256, ncol)

    device = array_device(p_lev)
    workgroup = (thr_x)
    ndrange = (ncol)

    kernel = compute_gray_heating_rate_kernel!(device, workgroup)

    event = kernel(
        flux_net,
        p_lev,
        hr_lay,
        grav_,
        cp_d_,
        Val(nlay),
        Val(ncol),
        ndrange = ndrange,
    )
    wait(event)
    #----------------------------------
    return nothing
end

@kernel function compute_gray_heating_rate_kernel!(
    flux_net::FTA2D,
    p_lev::FTA2D,
    hr_lay::FTA2D,
    grav_::FT,
    cp_d_::FT,
    ::Val{nlay},
    ::Val{ncol},
) where {FT<:AbstractFloat,FTA2D<:AbstractArray{FT,2},nlay,ncol}
    gcol = @index(Global, Linear)  # global col & lay ids

    for ilay = 1:nlay
        hr_lay[ilay, gcol] =
            grav_ * (flux_net[ilay+1, gcol] - flux_net[ilay, gcol]) /
            (p_lev[ilay+1, gcol] - p_lev[ilay, gcol]) / cp_d_
    end

end

end
