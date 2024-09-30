
function compute_col_gas!(
    device::ClimaComms.CUDADevice,
    p_lev::AbstractArray{FT, 2},
    col_dry::AbstractArray{FT, 2},
    param_set::RP.ARP,
    vmr_h2o::Union{AbstractArray{FT, 2}, Nothing} = nothing,
    lat::Union{AbstractArray{FT, 1}, Nothing} = nothing,
) where {FT}
    nlay, ncol = size(col_dry)
    mol_m_dry = RP.molmass_dryair(param_set)
    mol_m_h2o = RP.molmass_water(param_set)
    avogadro = RP.avogad(param_set)
    helmert1 = RP.grav(param_set)
    args = (p_lev, mol_m_dry, mol_m_h2o, avogadro, helmert1, vmr_h2o, lat)
    tx, bx = _configure_threadblock(nlay * ncol)
    @cuda always_inline = true threads = (tx) blocks = (bx) compute_col_gas_CUDA!(col_dry, args...)
    return nothing
end

function compute_col_gas_CUDA!(col_dry, args...)
    glx = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlay, ncol = size(col_dry)
    if glx ≤ nlay * ncol
        glay = (glx % nlay == 0) ? nlay : (glx % nlay)
        gcol = cld(glx, nlay)
        compute_col_gas_kernel!(col_dry, args..., glay, gcol)
    end
    return nothing
end

function compute_relative_humidity!(
    device::ClimaComms.CUDADevice,
    rh::AbstractArray{FT, 2},
    p_lay::AbstractArray{FT, 2},
    t_lay::AbstractArray{FT, 2},
    param_set::RP.ARP,
    vmr_h2o::AbstractArray{FT, 2},
) where {FT}
    nlay, ncol = size(p_lay)
    # ratio of water to dry air molecular weights
    mwd = RP.molmass_water(param_set) / RP.molmass_dryair(param_set)
    t_ref = FT(273.16) # reference temperature (K)
    q_lay_min = FT(1e-7) # minimum water mass mixing ratio

    args = (p_lay, t_lay, vmr_h2o, mwd, t_ref, q_lay_min)
    tx, bx = _configure_threadblock(nlay * ncol)
    @cuda always_inline = true threads = (tx) blocks = (bx) compute_relative_humidity_CUDA!(rh, args...)
    return nothing
end

function compute_relative_humidity_CUDA!(rh, args...)
    glx = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlay, ncol = size(rh)
    if glx ≤ nlay * ncol
        glay = (glx % nlay == 0) ? nlay : (glx % nlay)
        gcol = cld(glx, nlay)
        compute_relative_humidity_kernel!(rh, args..., glay, gcol)
    end
end
