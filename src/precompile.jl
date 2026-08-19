import PrecompileTools

# Precompile the gray-radiation solve at both working precisions. It touches
# the whole stack the standalone (Layer-3) API drives (profile construction,
# the Layer-2 `RRTMGPSolver`, `update_fluxes!`, the longwave/shortwave kernels,
# and the flux getters) while needing no lookup tables, so it runs without the
# NCDatasets extension or an artifact download. The spectral path shares those
# same layers, so its first call warms up quickly after this.
#
# The CPU device is pinned explicitly: the workload runs in the precompile
# subprocess, which inherits CLIMACOMMS_DEVICE, and `ClimaComms.context()`
# errors under CLIMACOMMS_DEVICE=CUDA because CUDA.jl is never loaded here.
# Only CPU code can precompile in any case.
PrecompileTools.@compile_workload begin
    context = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    for FT in (Float32, Float64)
        out = solve_gray(FT; context, nlay = 12, ncol = 2)
        net_flux(out.solver)
        heating_rate(out.solver)
    end
end
