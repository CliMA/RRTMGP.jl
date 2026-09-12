function rte_sw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::GrayAtmosphericState,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    tx, bx = _configure_threadblock(ncol)
    args = (flux_sw, op, bcs_sw, src_sw, nlay, ncol, as)
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux_sw::FluxSW,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay,
    ncol,
    as::GrayAtmosphericState,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt, igpt, ibnd = 1, 1, 1
    FT = eltype(bcs_sw.cos_zenith)
    solar_frac = FT(1)
    if gcol ≤ ncol
        μ₀ = bcs_sw.cos_zenith[gcol]
        @inbounds begin
            compute_optical_props!(op, as, gcol)
            # call shortwave rte solver
            rte_sw_2stream!(
                op,
                src_sw,
                bcs_sw,
                flux_sw,
                solar_frac,
                igpt,
                n_gpt,
                ibnd,
                nlev,
                gcol,
            )
            compute_net_flux!(flux_sw, gcol, nlev)
        end
        if μ₀ ≤ 0 # zero out columns with zenith angle ≥ π/2
            set_flux_to_zero!(flux_sw, gcol, nlev)
        end
    end
    return nothing
end

function rte_sw_2stream_solve!(
    device::ClimaComms.CUDADevice,
    flux::FluxSW,
    flux_sw::FluxSW,
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    nlay, ncol = AtmosphericStates.get_dims(as)
    set_band_flux_to_zero!(band_flux)
    tx, bx = _configure_threadblock(ncol)
    args = (
        flux,
        flux_sw,
        band_flux,
        op,
        bcs_sw,
        src_sw,
        nlay,
        ncol,
        as,
        state_cache,
        lookup_sw,
        lookup_sw_cld,
        lookup_sw_aero,
    )
    @cuda always_inline = true threads = (tx) blocks = (bx) rte_sw_2stream_solve_CUDA!(
        args...,
    )
    return nothing
end

function rte_sw_2stream_solve_CUDA!(
    flux::FluxSW,
    flux_sw::FluxSW,
    band_flux,
    op::TwoStream,
    bcs_sw::SwBCs,
    src_sw::SourceSW2Str,
    nlay,
    ncol,
    as::AtmosphericState,
    state_cache::Union{TransposedStateCache, Nothing},
    lookup_sw::LookUpSW,
    lookup_sw_cld::Union{LookUpCld, Nothing} = nothing,
    lookup_sw_aero::Union{LookUpAerosolMerra, Nothing} = nothing,
)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    nlev = nlay + 1
    n_gpt = length(lookup_sw.band_data.major_gpt2bnd)
    if gcol ≤ ncol
        FT = eltype(flux_sw.flux_up)
        (; cloud_state, aerosol_state) = as
        μ₀ = @inbounds bcs_sw.cos_zenith[gcol]
        n_cloudy_gpts = 0  # thread-local counter for cloud cover
        @inbounds begin
            # EXPERIMENT 2026-09-11: skip the solve entirely for night columns.
            #
            # The test is `> eps(FT)`, NOT `> 0`. ClimaAtmos writes
            # `cos_zenith = max(insolation.μ, eps(FT))` because RRTMGP requires a
            # non-zero μ, so `μ₀ ≤ 0` never occurs and RRTMGP's own
            # `set_flux_to_zero!` branch is dead code under this configuration.
            # A first attempt testing `μ₀ > 0` therefore skipped nothing and
            # measured -0.59%, inside the +-0.4% noise.
            #
            # The g-point loop below runs ~224 full optics-plus-vertical-solve
            # evaluations per column, and it ran for EVERY column -- including
            # those with the sun below the horizon, whose result was then
            # discarded by `set_flux_to_zero!`. At any instant roughly half the
            # globe is dark, so about half this kernel's work was thrown away.
            #
            # This should be warp-coherent, which is what decides whether a GPU
            # early-out pays: columns are geographically ordered, so day and
            # night cluster rather than alternating. A criterion true at
            # scattered points saves nothing, because a warp costs the maximum
            # over its 32 lanes -- that is why the microphysics clear-air
            # early-out returned only +1.79% off a 77.7% point fraction with a
            # 21.5% warp fraction. ncu measures 24.72 of 32 active threads per
            # warp here, consistent with a large coherent inactive fraction.
            #
            # BEHAVIOURAL NOTE: `cld_cover_sw` becomes 0 for night columns
            # instead of the fraction diagnosed from a solve whose fluxes were
            # discarded anyway. Shortwave cloud cover at night is arguably
            # undefined; it is still a diagnostic change and must be checked.
            if μ₀ > eps(FT)
            _compute_aero_mask!(aerosol_state, gcol)
            for igpt in 1:n_gpt
                cloudy = sw_2stream_gpt_col!(
                    igpt,
                    gcol,
                    flux,
                    flux_sw,
                    band_flux,
                    op,
                    bcs_sw,
                    src_sw,
                    as,
                    state_cache,
                    lookup_sw,
                    lookup_sw_cld,
                    lookup_sw_aero,
                    μ₀,
                    lookup_sw.band_data.major_gpt2bnd[igpt],
                    n_gpt,
                    nlev,
                )
                n_cloudy_gpts += cloudy ? 1 : 0
            end
            compute_net_flux!(flux_sw, gcol, nlev)
            else # sun below the horizon: nothing to solve
                set_flux_to_zero!(flux_sw, gcol, nlev)
            end
            # write out SW cloud cover
            if cloud_state isa CloudState &&
               !isnothing(cloud_state.cld_cover_sw)
                cloud_state.cld_cover_sw[gcol] = FT(n_cloudy_gpts) / n_gpt
            end
        end
    end
    return nothing
end
