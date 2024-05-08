function rte_lw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux_lw::FluxLW,
    src_lw::SourceLWNoScat,
    bcs_lw::LwBCs,
    op::OneScalar,
    max_threads,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux_lw, src_lw, bcs_lw, op, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_noscat_solve_CUDA!(args...)
    return nothing
end

function rte_lw_noscat_solve_CUDA!(
    flux_lw::FluxLW,
    src_lw::SourceLWNoScat,
    bcs_lw::LwBCs,
    op::OneScalar,
    nlay,
    ncol,
    as::GrayAtmosphericState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    igpt, ibnd = 1, 1
    τ = op.τ
    Ds = op.angle_disc.gauss_Ds
    (; flux_up, flux_dn, flux_net) = flux_lw
    if gcol ≤ ncol
        compute_optical_props!(op, as, src_lw, gcol)
        rte_lw_noscat!(src_lw, bcs_lw, op, gcol, flux_lw, igpt, ibnd, nlay, nlev)
        @inbounds for ilev in 1:nlev
            flux_net[ilev, gcol] = flux_up[ilev, gcol] - flux_dn[ilev, gcol]
        end
    end
    return nothing
end

function rte_lw_noscat_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxLW,
    flux_lw::FluxLW,
    src_lw::SourceLWNoScat,
    bcs_lw::LwBCs,
    op::OneScalar,
    max_threads,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    args = (flux, flux_lw, src_lw, bcs_lw, op, nlay, ncol, as, lookup_lw, lookup_lw_cld)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_lw_noscat_solve_CUDA!(args...)
    return nothing
end

function rte_lw_noscat_solve_CUDA!(
    flux::FluxLW,
    flux_lw::FluxLW,
    src_lw::SourceLWNoScat,
    bcs_lw::LwBCs,
    op::OneScalar,
    nlay,
    ncol,
    as::AtmosphericState,
    lookup_lw::LookUpLW,
    lookup_lw_cld::Union{LookUpCld, PadeCld, Nothing} = nothing,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    (; major_gpt2bnd) = lookup_lw.band_data
    n_gpt = length(major_gpt2bnd)
    τ = op.τ
    Ds = op.angle_disc.gauss_Ds
    if gcol ≤ ncol
        flux_up_lw = flux_lw.flux_up
        flux_dn_lw = flux_lw.flux_dn
        flux_net_lw = flux_lw.flux_net
        @inbounds for igpt in 1:n_gpt
            ibnd = major_gpt2bnd[igpt]
            igpt == 1 && set_flux_to_zero!(flux_lw, gcol)
            compute_optical_props!(op, as, src_lw, gcol, igpt, lookup_lw, lookup_lw_cld)
            rte_lw_noscat!(src_lw, bcs_lw, op, gcol, flux, igpt, ibnd, nlay, nlev)
            add_to_flux!(flux_lw, flux, gcol)
        end
        @inbounds begin
            for ilev in 1:nlev
                flux_net_lw[ilev, gcol] = flux_up_lw[ilev, gcol] - flux_dn_lw[ilev, gcol]
            end
        end
    end
    return nothing
end
