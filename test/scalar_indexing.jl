using Test
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import RRTMGP
import CUDA

# The `(nlev, ncol)` flux buffers are exposed through `_domain_view` (src/api/getters.jl)
# as plain domain-masked `SubArray` views of the device array. Bringing one to the CPU
# with `Array(getter(solver))` must not fall back to element-wise GPU scalar indexing
# (CUDA.jl is expected to convert these strided views natively). This test pins that
# contract for the flux getters; if it fails on GPU, an RRTMGP-owned `to_cpu` is the
# fix -- not a `Base.Array` overload on foreign types.

# A gray two-stream lw+sw solver with an optional isothermal boundary layer.
function gray_solver(
    ::Type{FT},
    context,
    isothermal_boundary_layer::Bool,
) where {FT}
    DA = ClimaComms.array_type(ClimaComms.device(context))
    domain_nlay, ncol = 40, 4
    gp = RRTMGP.RRTMGPGridParams(
        FT;
        context,
        domain_nlay,
        ncol,
        isothermal_boundary_layer,
    )
    nlay = gp.nlay # total (boundary-extended) layer count
    params = RRTMGP.default_parameters(FT)
    lat = DA{FT}(range(FT(-80), FT(80); length = ncol))
    as = RRTMGP.AtmosphericStates.setup_gray_as_pr_grid(
        context,
        nlay,
        lat,
        FT(1.0e5),
        FT(9.0e3),
        RRTMGP.AtmosphericStates.GrayOpticalThicknessOGorman2008(FT),
        params,
        DA,
    )
    sfc_emis = fill!(DA{FT}(undef, 1, ncol), FT(1))
    bcs_lw = RRTMGP.BCs.LwBCs(sfc_emis, nothing)
    cos_zen = fill!(DA{FT}(undef, ncol), FT(0.5))
    toa = fill!(DA{FT}(undef, ncol), FT(1361))
    alb = fill!(DA{FT}(undef, 1, ncol), FT(0.2))
    bcs_sw = RRTMGP.BCs.SwBCs(cos_zen, toa, alb, nothing, copy(alb))
    solver = RRTMGP.RRTMGPSolver(
        gp,
        RRTMGP.GrayRadiation(),
        params,
        bcs_lw,
        bcs_sw,
        as,
    )
    RRTMGP.update_fluxes!(solver)
    return solver, domain_nlay, ncol
end

@testset "flux getters materialize to CPU without scalar indexing" begin
    context = ClimaComms.context()
    if ClimaComms.device(context) isa ClimaComms.CUDADevice
        # Disallow GPU scalar indexing (CUDA.jl's default) so that any `getter` that
        # materialized through the generic element-wise `Array` fallback -- rather than
        # the ext's bulk-copy `Base.Array` overload -- errors here instead of silently
        # running slowly. Left set for the remainder of the suite (the strict default).
        CUDA.allowscalar(false)
        getters = (
            RRTMGP.net_flux,
            RRTMGP.lw_flux_up,
            RRTMGP.lw_flux_dn,
            RRTMGP.lw_flux_net,
            RRTMGP.sw_flux_up,
            RRTMGP.sw_flux_dn,
            RRTMGP.sw_flux_net,
        )
        for FT in (Float32, Float64), boundary in (false, true)
            solver, domain_nlay, ncol = gray_solver(FT, context, boundary)
            for getter in getters
                arr = Array(getter(solver))  # must not scalar-index the GPU
                @test arr isa Array{FT, 2}
                @test size(arr) == (domain_nlay + 1, ncol)
                @test all(isfinite, arr)
            end
        end
    else
        @info "scalar-indexing guard skipped: no CUDA device"
    end
end
