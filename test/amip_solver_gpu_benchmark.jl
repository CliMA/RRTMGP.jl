#####
##### The call a coupled model actually makes, at the shape it makes it.
#####
##### The other high-resolution benchmarks here time `solve_lw!` and `solve_sw!`
##### directly, on a column count that already fills the device. ClimaAtmos does
##### not call them that way: it builds one `RRTMGPSolver` and calls
##### `update_fluxes!`, which solves both bands, samples the cloud optics and
##### carries the clear-sky diagnostics -- and the solver's own keyword defaults
##### are reached only through that constructor, so a benchmark that bypasses it
##### cannot see them. This times `update_fluxes!` at the grid and the flags
##### ClimaAtmos's AMIP benchmark configuration runs.
#####

# `_resample_benchmark_as`, the lookup readers and the state setup, with the
# `using` block they need. The sweep loops in these files are guarded behind
# `PROGRAM_FILE == @__FILE__`, so including one only defines functions.
if !@isdefined(_resample_benchmark_as)
    include(
        joinpath(@__DIR__, "all_sky_with_aerosols_highres_gpu_benchmark.jl"),
    )
end

# The AMIP benchmark grid: `h_elem: 16` on a cubed sphere at fourth order, which
# is 16^2 elements over 6 panels at 4^2 quadrature points, and `z_elem: 79`.
amip_ncol(; helems = 16, nq = 4) = Int(helems * helems * 6 * nq * nq)
const amip_domain_nlay = 79

"""
    benchmark_amip_solver(context, FT; ncol, domain_nlay, cldfrac)

Time one `update_fluxes!` on an `RRTMGPSolver` built the way ClimaAtmos builds
it for AMIP, returning the `BenchmarkTools` trial.

Everything here is a configuration ClimaAtmos sets: `rad: "allskywithclear"` is
`AllSkyRadiationWithClearSkyDiagnostics`, `aerosol_radiation: true` is its first
flag, two-stream optics in both bands, and the isothermal layer that makes
`z_elem: 79` into the 80 layers the solve works over.

The state comes from the all-sky-with-aerosols test file, resampled to the
grid's layer count and repeated across columns, so it is this grid's shape
rather than this grid's weather. Cost is what is being measured, and cloud
fraction is held at one so every column sweeps its cloudy increment.
"""
function benchmark_amip_solver(
    context,
    ::Type{FT};
    ncol = amip_ncol(),
    domain_nlay = amip_domain_nlay,
    cldfrac = FT(1),
) where {FT <: AbstractFloat}
    overrides =
        (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)
    device = ClimaComms.device(context)

    radiation_method = RRTMGP.AllSkyRadiationWithClearSkyDiagnostics(
        true, # aerosol_radiation
        true, # reset_rng_seed
    )
    grid_params = RRTMGPGridParams(
        FT;
        context,
        domain_nlay,
        ncol,
        isothermal_boundary_layer = true,
    )
    # The public loader, which is how a coupled model gets its tables. It also
    # carries the name-to-index maps the state setup below asks for, so the
    # eight lookup files are read once and by one call.
    lookups = RRTMGP.lookup_tables(grid_params, radiation_method)

    ds_in = Dataset(get_input_filename(:gas_clouds_aerosols, :lw), "r")
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, cos_zenith, toa_flux, _ =
        setup_allsky_with_aerosols_as(
            context,
            ds_in,
            lookups.idx_gases_lw,
            lookups.idx_aerosol_lw,
            lookups.idx_aerosize_lw,
            lookups.lookup_lw,
            lookups.lookup_sw,
            lookups.lookup_lw_cld,
            lookups.lookup_sw_cld,
            cldfrac,
            ncol,
            FT,
            param_set,
        )
    close(ds_in)
    as = _resample_benchmark_as(context, as, grid_params.nlay, FT)

    inc_flux = nothing
    inc_flux_diffuse = nothing
    bcs_lw = BCs.LwBCs(sfc_emis, inc_flux)
    bcs_sw = BCs.SwBCs(
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
    solver = RRTMGP.RRTMGPSolver(
        grid_params,
        radiation_method,
        param_set,
        bcs_lw,
        bcs_sw,
        as;
        op_lw = Optics.TwoStream(grid_params),
        op_sw = Optics.TwoStream(grid_params),
        lookups,
    )

    # One seed for every sample, so the McICA draw is fixed and what varies
    # between samples is the machine rather than the cloud field.
    seed = UInt32(42)
    RRTMGP.update_fluxes!(solver, seed) # warm up / compile
    return if device isa ClimaComms.CUDADevice
        @benchmark CUDA.@sync RRTMGP.update_fluxes!($solver, $seed)
    else
        @benchmark RRTMGP.update_fluxes!($solver, $seed)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    context = ClimaComms.context()
    ClimaComms.init(context)
    for FT in (Float32, Float64)
        trial = benchmark_amip_solver(context, FT)
        println("amip update_fluxes! ($FT): ", minimum(trial))
    end
end
