using Test
using Pkg.Artifacts
using NCDatasets
using Statistics
using BenchmarkTools
using Printf

import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

using RRTMGP
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
include("read_all_sky.jl")

function benchmark_all_sky(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{FT};
    ncol = 128,# repeats col#1 ncol times per RRTMGP example 
    use_lut::Bool = true,
    cldfrac = FT(1),
) where {FT <: AbstractFloat, SLVLW, SLVSW}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)

    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    n_gauss_angles = 1

    # read lookup_data
    lookup_lw, lookup_sw, lookup_lw_cld, lookup_sw_cld, idx_gases = read_all_sky_lookup_files(DA, FT, use_lut)

    # reading input file 
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, cos_zenith, toa_flux, bot_at_1 = setup_allsky_as(
        context,
        idx_gases,
        lookup_lw,
        lookup_sw,
        lookup_lw_cld,
        lookup_sw_cld,
        cldfrac,
        use_lut,
        ncol,
        FT,
        param_set,
    )
    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    #---reading comparison files -----------------------------------
    if bot_at_1
        flip_ind = 1:nlev
    else
        flip_ind = nlev:-1:1
    end

    # Setting up longwave problem---------------------------------------
    inc_flux = nothing
    slv_lw = SLVLW(FT, DA, context, param_set, nlay, ncol, sfc_emis, inc_flux)
    # Setting up shortwave problem---------------------------------------
    inc_flux_diffuse = nothing
    swbcs = (cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(FT, DA, context, nlay, ncol, swbcs...)
    #------calling solvers
    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld)
    trial_lw = @benchmark CUDA.@sync solve_lw!($slv_lw, $as, $lookup_lw, $lookup_lw_cld)

    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld)
    trial_sw = @benchmark CUDA.@sync solve_sw!($slv_sw, $as, $lookup_sw, $lookup_sw_cld)
    return trial_lw, trial_sw
end

function generate_gpu_allsky_benchmarks(FT, npts, ::Type{SLVLW}, ::Type{SLVSW}) where {SLVLW, SLVSW}
    context = ClimaComms.context()
    # compute equivalent ncols for DYAMOND resolution
    helems, nlevels, nlev_test, nq = 30, 64, 73, 4
    ncols_dyamond = Int(ceil(helems * helems * 6 * nq * nq * (nlevels / nlev_test)))
    println("\n")
    printstyled("Running DYAMOND all-sky benchmark on $(context.device) device with $FT precision\n", color = 130)
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
        trial_lw, trial_sw = benchmark_all_sky(context, SLVLW, SLVSW, FT; ncol = ncols, use_lut = true, cldfrac = FT(1))
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
    generate_gpu_allsky_benchmarks(FT, 4, NoScatLWRTE, TwoStreamSWRTE)
    generate_gpu_allsky_benchmarks(FT, 4, TwoStreamLWRTE, TwoStreamSWRTE)
end
