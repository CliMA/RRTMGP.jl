using Test
using Pkg.Artifacts
using NCDatasets
using Statistics
using BenchmarkTools
using Printf

import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

using RRTMGP
using RRTMGP: RRTMGPGridParams
using RRTMGP.Vmrs
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.AngularDiscretizations
using RRTMGP.RTE
using RRTMGP.RTESolver
import RRTMGP.Parameters.RRTMGPParameters
import ClimaParams as CP
using RRTMGP.ArtifactPaths
# overriding some parameters to match with RRTMGP FORTRAN code

include("reference_files.jl")
include("read_clear_sky.jl")

function benchmark_clear_sky(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{VMR},
    ::Type{FT};
    ncol = 100,
) where {FT <: AbstractFloat, SLVLW, SLVSW, VMR}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)

    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    n_gauss_angles = 1
    expt_no = 1

    lw_file = get_lookup_filename(:gas, :lw) # lw lookup tables for gas optics
    sw_file = get_lookup_filename(:gas, :sw) # sw lookup tables for gas optics

    #reading longwave gas optics lookup data
    lookup_lw, idx_gases = Dataset(lw_file, "r") do ds
        LookUpLW(ds, FT, DA)
    end

    #reading shortwave gas optics lookup data
    lookup_sw, idx_gases = Dataset(sw_file, "r") do ds
        LookUpSW(ds, FT, DA)
    end

    # reading input file
    input_file = get_input_filename(:gas, :lw) # all-sky atmos state
    ds_lw_in = Dataset(input_file, "r")
    (as, sfc_emis, sfc_alb_direct, cos_zenith, toa_flux, bot_at_1) =
        setup_clear_sky_as(context, ds_lw_in, idx_gases, expt_no, lookup_lw, ncol, FT, VMR, param_set)
    close(ds_lw_in)

    nlay, _ = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)


    # setting up longwave problem
    inc_flux = nothing
    slv_lw = SLVLW(grid_params; params = param_set, sfc_emis, inc_flux)

    # setting up shortwave problem
    sfc_alb_diffuse = FTA2D(deepcopy(sfc_alb_direct))
    inc_flux_diffuse = nothing
    swbcs = (; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(grid_params; swbcs...)
    #------calling solvers
    solve_lw!(slv_lw, as, lookup_lw)
    trial_lw = @benchmark CUDA.@sync solve_lw!($slv_lw, $as, $lookup_lw)

    solve_sw!(slv_sw, as, lookup_sw)
    trial_sw = @benchmark CUDA.@sync solve_sw!($slv_sw, $as, $lookup_sw)
    return trial_lw, trial_sw
end

function generate_gpu_clear_sky_benchmarks(
    FT,
    npts,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{VMR},
) where {SLVLW, SLVSW, VMR}
    context = ClimaComms.context()
    # compute equivalent ncols for DYAMOND resolution
    helems, nlevels, nlev_test, nq = 30, 64, 61, 4
    ncols_dyamond = Int(ceil(helems * helems * 6 * nq * nq * (nlevels / nlev_test)))
    println("\n")
    printstyled("Running clear-sky benchmark on $(context.device) device with $FT precision\n", color = 130)
    printstyled("Longwave solver = $SLVLW; Shortwave solver = $SLVSW\n", color = 130)
    printstyled("==============|====================================|==================================\n", color = 130)
    printstyled(
        "  ncols       |   median time for longwave solver  | median time for shortwave solver \n",
        color = :green,
    )
    printstyled("==============|====================================|==================================\n", color = 130)
    for pts in 1:npts
        ncols = unsafe_trunc(Int, cld(ncols_dyamond, 2^(pts - 1)))
        ndof = ncols * nlev_test
        sz_per_fld_gb = ndof * sizeof(FT) / 1024 / 1024 / 1024
        trial_lw, trial_sw = benchmark_clear_sky(context, SLVLW, SLVSW, VMR, FT; ncol = ncols)
        Printf.@printf(
            "%10i    |           %25s|       %25s \n",
            ncols,
            Statistics.median(trial_lw),
            Statistics.median(trial_sw)
        )
    end
    printstyled("==============|====================================|==================================\n", color = 130)
    return nothing
end

for FT in (Float32, Float64)
    generate_gpu_clear_sky_benchmarks(FT, 4, NoScatLWRTE, TwoStreamSWRTE, VmrGM)
    generate_gpu_clear_sky_benchmarks(FT, 4, TwoStreamLWRTE, TwoStreamSWRTE, VmrGM)
end
