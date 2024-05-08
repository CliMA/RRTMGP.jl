
function update_profile_lw!(
    device::ClimaComms.CUDADevice,
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
    max_threads = 256
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    @cuda always_inline = true threads = (tx) blocks = (bx) _update_profile_lw_kernel!(
        sbc,
        ncol,
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
end

function _update_profile_lw_kernel!(sbc, ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        update_profile_lw_kernel!(sbc, args..., gcol)
    end
    return nothing
end

function compute_gray_heating_rate!(device::ClimaComms.CUDADevice, hr_lay, p_lev, ncol, nlay, flux_net, cp_d_, grav_)
    args = (hr_lay, flux_net, p_lev, grav_, cp_d_, nlay)
    max_threads = 256
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    @cuda always_inline = true threads = (tx) blocks = (bx) compute_gray_heating_rate_CUDA!(ncol, args...)
    return nothing
end

function compute_gray_heating_rate_CUDA!(ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        compute_gray_heating_rate_kernel!(args..., gcol)
    end
    return nothing
end
