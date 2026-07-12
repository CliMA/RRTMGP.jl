# How to run on GPUs

RRTMGP.jl runs the same code on CPUs and CUDA GPUs; the device is selected
through [ClimaComms](https://github.com/CliMA/ClimaComms.jl), and every array
the solver owns lives on that device.

## Select the device

Set the device before constructing anything, either with an environment
variable:

```
CLIMACOMMS_DEVICE=CUDA julia --project my_script.jl
```

or explicitly in code:

```julia
import ClimaComms
ClimaComms.@import_required_backends  # loads CUDA.jl when a GPU is requested

context = ClimaComms.context()        # honors CLIMACOMMS_DEVICE
device = ClimaComms.device(context)
DA = ClimaComms.array_type(device)    # Array or CuArray
```

Pass `context` to [`RRTMGPGridParams`](@ref RRTMGP.RRTMGPGridParams) (or to
`solve_gray`/`solve`) and build your input arrays with `DA{FT}(...)`; from
there, everything (states, workspaces, solves) stays on the GPU.

## Reading results without scalar indexing

Getters return device-array views. Indexing them elementwise from the host
(`net_flux(solver)[1, 1]`) would trigger CUDA scalar indexing; instead
materialize what you need:

```julia
F = Array(RRTMGP.net_flux(solver))   # one device→host copy
```

Broadcasted writes into getters (`layer_temperature(solver) .= T_dev`) are
device-side and safe. Two getters do not return views:
`volume_mixing_ratio(solver, name)` for well-mixed gases returns a host
scalar, and `heating_rate(solver)` allocates a fresh device array; read them
outside hot loops.

## Checkpointing and device transfer

`RRTMGPSolver` supports [Adapt](https://github.com/JuliaGPU/Adapt.jl), so
`Adapt.adapt(Array, solver)` produces a host copy (e.g., for serialization) and
`Adapt.adapt(CuArray, solver)` moves it back. RRTMGP's custom `adapt_structure`
methods preserve the internal view topology (scratch views into packed
buffers), so an adapted solver keeps the zero-allocation property.

## Reproducibility caveat: McICA cloud sampling

With partial cloud fractions, the all-sky methods sample cloud overlap
stochastically (McICA) using the global RNG. On a single CPU thread,
`update_fluxes!(solver, seed)` with `reset_rng_seed = true` makes results
reproducible. In multi-threaded settings and on the GPU (at least on CUDA 6.x),
the per-thread RNG state produces statistically but not deterministically
reproducible results for partially cloudy columns (overcast and clear columns
remain deterministic).

## Parallelization strategy

Each column and g-point is an independent monochromatic transfer problem; only
the vertical sweeps within a column are sequential. RRTMGP parallelizes over
this independence, with the same per-(g-point, column) kernel bodies shared by
both backends:

- On the GPU, each CUDA thread owns one column and loops over its g-points,
  accumulating the column's broadband fluxes thread-locally — no atomics and no
  synchronization between threads. Parallelism therefore scales with `ncol`,
  and a GPU pays off with many columns, the regime the benchmarks below target.
- On the CPU, the loop order is inverted: for each g-point, the columns are
  distributed across threads (`ClimaComms.@threaded`), which keeps that
  g-point's lookup-table slices hot in cache while sweeping the columns. With
  a single-threaded context, the same code runs serially.

The memory layout matches the thread mapping: the internal scratch and output
buffers are stored column-first, `(ncol, nlev)`, so consecutive GPU threads
write consecutive addresses (coalesced access), and the frequently re-read
state arrays (layer pressures and temperatures, level temperatures) are
copied into a column-first `TransposedStateCache` once per solve so the
g-point loop also reads coalesced memory. The caller-owned
`AtmosphericState` and the getter views keep their vertical-first
`(nlev, ncol)` presentation.

## Performance expectations

CI benchmarks the kernels on GPUs (see `perf/benchmark_ratchet.jl`
and the Buildkite "Benchmarks" group). Indicative behavior:

- Work scales with `ncol × nlay × ngpt`; the g-point loop dominates.
- `update_fluxes!` allocates nothing on either device; data movement is the
  host's choice (the getters return views).
- `Float32` roughly halves memory traffic; see
  [Float32 and Float64](../precision.md) for the accuracy achieved at single
  precision and the numerics behind it.

## GPU versus CPU performance

Wall times per full radiation step (`solve_lw!` + `solve_sw!`) comparing the
multi-threaded CPU strategy against the GPU strategy across single (`Float32`) and double (`Float64`) precision.

### Benchmark problem dimensions

Every configuration evaluated by the benchmark (see [`perf/benchmark_ratchet.jl`](https://github.com/CliMA/RRTMGP.jl/blob/main/perf/benchmark_ratchet.jl)) operates on an identical, uniform spatial grid representative of a high-resolution atmospheric model:

- **Columns (`ncol`)**: **`86,400` horizontal columns** (`30² elements × 6 panels × 4² quadrature points`).
- **Vertical structure**: **`63` layers (`64` levels)** across every test profile (`5.53 × 10⁶` uniform 3D spatial coordinates, `ncol × nlev`), aligned with GPU warp sizes and memory boundaries.
- **Spectral quadrature points (`ngpt`)**: **`480` spectral g-points** (`256` longwave + `224` shortwave) across all problems, resulting in approximately `2.65 × 10⁹` total cell evaluations per `solve_lw!` + `solve_sw!` step.

### Single precision (`Float32`) — Primary target

At `Float32` precision, memory traffic across lookup-table interpolations and solver evaluations is halved compared to `Float64`, enabling massive GPU speedups when combined with contiguous `(ncol, nlev)` internal memory coalescing:

| Configuration | CPU (Intel Xeon, 12 threads) | GPU (NVIDIA A100 40GB) | Speedup |
|---|---|---|---|
| Clear-sky (two-stream) | 77.95 s | 0.95 s | 82.5× |
| All-sky, McICA clouds | 98.42 s | 1.39 s | 70.8× |
| All-sky with aerosols | 99.39 s | 1.52 s | 65.5× |

### Double precision (`Float64`)

For reference, running the identical uniform benchmark sweep in double precision (`Float64`) yields:

| Configuration | CPU (Intel Xeon, 12 threads) | GPU (NVIDIA A100 40GB) | Speedup |
|---|---|---|---|
| Clear-sky (two-stream) | 79.64 s | 1.32 s | 60.1× |
| All-sky, McICA clouds | 101.60 s | 1.93 s | 52.6× |
| All-sky with aerosols | 106.33 s | 2.07 s | 51.4× |
