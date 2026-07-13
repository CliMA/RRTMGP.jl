"""
    compute_col_gas!(
        device,
        p_lev,
        col_dry,
        param_set,
        vmr_h2o,
        lat,
        max_threads = Int(256),
    )

Compute the column amounts of dry or moist air.
"""
function compute_col_gas!(
    device::ClimaComms.AbstractCPUDevice,
    p_lev::AbstractArray{FT, 2},
    col_dry::AbstractArray{FT, 2},
    param_set::RP.ARP,
    vmr_h2o::Union{AbstractArray{FT, 2}, Nothing} = nothing,
    lat::Union{AbstractArray{FT, 1}, Nothing} = nothing,
    max_threads::Int = Int(256),
) where {FT}
    nlay, ncol = size(col_dry)
    mol_m_dry = RP.molmass_dryair(param_set)
    mol_m_h2o = RP.molmass_water(param_set)
    avogadro = RP.avogad(param_set)
    helmert1 = RP.grav(param_set)
    args = (p_lev, mol_m_dry, mol_m_h2o, avogadro, helmert1, vmr_h2o, lat)
    @inbounds begin
        ClimaComms.@threaded device for icnt in 1:(nlay * ncol)
            gcol = cld(icnt, nlay)
            glay = (icnt % nlay == 0) ? nlay : (icnt % nlay)
            compute_col_gas_kernel!(col_dry, args..., glay, gcol)
        end
    end
    return nothing
end

"""
    compute_relative_humidity!(
        device::ClimaComms.AbstractCPUDevice,
        rh::AbstractArray{FT, 2},
        p_lay::AbstractArray{FT, 2},
        t_lay::AbstractArray{FT, 2},
        param_set::RP.ARP,
        vmr_h2o::Union{AbstractArray{FT, 2}, Nothing} = nothing,
    ) where {FT}

Compute the relative humidity.

"""
function compute_relative_humidity!(
    device::ClimaComms.AbstractCPUDevice,
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

    args = (rh, p_lay, t_lay, vmr_h2o, mwd, t_ref, q_lay_min)
    @inbounds begin
        ClimaComms.@threaded device for icnt in 1:(nlay * ncol)
            gcol = cld(icnt, nlay)
            glay = (icnt % nlay == 0) ? nlay : (icnt % nlay)
            compute_relative_humidity_kernel!(args..., glay, gcol)
        end
    end
    return nothing
end

