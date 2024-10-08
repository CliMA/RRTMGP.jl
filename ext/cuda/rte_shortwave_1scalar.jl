function rte_sw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx, bx = _configure_threadblock(ncol)
    args = (flux_sw, op, bcs_sw, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    return nothing
end

function rte_sw_noscat_solve_CUDA!(flux_sw::FluxSW, op::OneScalar, bcs_sw::SwBCs, nlay, ncol, as::GrayAtmosphericState)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt = 1, 1
    FT = eltype(bcs_sw.cos_zenith)
    solar_frac = FT(1)
    if gcol ≤ ncol
        flux_up_sw = flux_sw.flux_up
        flux_dn_sw = flux_sw.flux_dn
        flux_net_sw = flux_sw.flux_net
        @inbounds begin
            μ₀ = bcs_sw.cos_zenith[gcol]
            if μ₀ > 0
                compute_optical_props!(op, as, gcol)
                rte_sw_noscat!(flux_sw, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
            else
                set_flux_to_zero!(flux_sw, gcol)
            end
        end
    end
    return nothing
end

function rte_sw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx, bx = _configure_threadblock(ncol)
    args = (flux, flux_sw, op, bcs_sw, nlay, ncol, as, lookup_sw)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_noscat_solve_CUDA!(args...)
    return nothing
end

function rte_sw_noscat_solve_CUDA!(
    flux::FluxSW,
    flux_sw::FluxSW,
    op::OneScalar,
    bcs_sw::SwBCs,
    nlay,
    ncol,
    as::AtmosphericState,
    lookup_sw::LookUpSW,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.solar_src_scaled)
    if gcol ≤ ncol
        flux_up_sw = flux_sw.flux_up
        flux_dn_sw = flux_sw.flux_dn
        flux_net_sw = flux_sw.flux_net
        flux_dn_dir_sw = flux_sw.flux_dn_dir
        flux_up = flux.flux_up
        flux_dn = flux.flux_dn
        flux_dn_dir = flux.flux_dn_dir
        μ₀ = bcs_sw.cos_zenith[gcol]
        @inbounds begin
            for igpt in 1:n_gpt
                compute_optical_props!(op, as, gcol, igpt, lookup_sw, nothing)
                solar_frac = lookup_sw.solar_src_scaled[igpt]
                rte_sw_noscat!(flux, op, bcs_sw, igpt, n_gpt, solar_frac, gcol, nlev)
                if igpt == 1
                    map!(x -> x, view(flux_up_sw, :, gcol), view(flux_up, :, gcol))
                    map!(x -> x, view(flux_dn_sw, :, gcol), view(flux_dn, :, gcol))
                    map!(x -> x, view(flux_dn_dir_sw, :, gcol), view(flux_dn_dir, :, gcol))
                else
                    for ilev in 1:nlev
                        flux_up_sw[ilev, gcol] += flux_up[ilev, gcol]
                        flux_dn_sw[ilev, gcol] += flux_dn[ilev, gcol]
                    end
                    flux_dn_dir_sw[1, gcol] += flux_dn_dir[1, gcol]
                end
            end
            if μ₀ <= 0
                for ilev in 1:nlev
                    flux_up_sw[ilev, gcol] = FT(0)
                end
                for ilev in 1:nlev
                    flux_dn_sw[ilev, gcol] = FT(0)
                end
                for ilev in 1:nlev
                    flux_net_sw[ilev, gcol] = FT(0)
                end
            else
                for ilev in 1:nlev
                    flux_net_sw[ilev, gcol] = flux_up_sw[ilev, gcol] - flux_dn_sw[ilev, gcol]
                end
            end
        end
    end
    return nothing
end
