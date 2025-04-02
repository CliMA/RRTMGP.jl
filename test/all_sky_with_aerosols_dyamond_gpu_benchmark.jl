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
include("read_all_sky_with_aerosols.jl")

function benchmark_all_sky_with_aerosols(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{FT};
    ncol = 128,# repeats col#1 ncol times per RRTMGP example
    cldfrac = FT(1),
) where {FT <: AbstractFloat, SLVLW, SLVSW}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)

    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}
    n_gauss_angles = 1

    lw_file = get_lookup_filename(:gas, :lw)          # lw lookup tables for gas optics
    lw_cld_file = get_lookup_filename(:cloud, :lw)    # lw cloud lookup tables
    lw_aero_file = get_lookup_filename(:aerosol, :lw) # lw aerosol lookup tables
    sw_file = get_lookup_filename(:gas, :sw)          # sw lookup tables for gas optics
    sw_cld_file = get_lookup_filename(:cloud, :sw)    # sw cloud lookup tables
    sw_aero_file = get_lookup_filename(:aerosol, :sw) # sw aerosol lookup tables

    input_file = get_input_filename(:gas_clouds_aerosols, :lw) # all-sky atmos state

    #reading longwave gas optics lookup data
    lookup_lw, idx_gases = Dataset(lw_file, "r") do ds
        LookUpLW(ds, FT, DA)
    end
    # reading longwave cloud lookup data
    lookup_lw_cld = Dataset(lw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end
    # reading longwave aerosol lookup data
    lookup_lw_aero, idx_aerosol, idx_aerosize = Dataset(lw_aero_file, "r") do ds
        LookUpAerosolMerra(ds, FT, DA)
    end

    #reading shortwave gas optics lookup data
    lookup_sw, idx_gases = Dataset(sw_file, "r") do ds
        LookUpSW(ds, FT, DA)
    end
    # reading longwave cloud lookup data
    lookup_sw_cld = Dataset(sw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end
    # reading shortwave aerosol lookup data
    lookup_sw_aero, _, _ = Dataset(sw_aero_file, "r") do ds
        LookUpAerosolMerra(ds, FT, DA)
    end

    # reading input file
    ds_in = Dataset(input_file, "r")
    as, sfc_emis, sfc_alb_direct, sfc_alb_diffuse, cos_zenith, toa_flux, bot_at_1 = setup_allsky_with_aerosols_as(
        context,
        ds_in,
        idx_gases,
        idx_aerosol,
        idx_aerosize,
        lookup_lw,
        lookup_sw,
        lookup_lw_cld,
        lookup_sw_cld,
        cldfrac,
        ncol,
        FT,
        param_set,
    )
    close(ds_in)

    nlay, _ = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)
    # Setting up longwave problem---------------------------------------
    inc_flux = nothing
    slv_lw = SLVLW(grid_params; params = param_set, sfc_emis, inc_flux)
    # Setting up shortwave problem---------------------------------------
    inc_flux_diffuse = nothing
    swbcs = (; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(grid_params; swbcs...)
    #------calling solvers
    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, lookup_lw_aero)
    trial_lw = @benchmark CUDA.@sync solve_lw!($slv_lw, $as, $lookup_lw, $lookup_lw_cld, $lookup_lw_aero)

    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld, lookup_sw_aero)
    trial_sw = @benchmark CUDA.@sync solve_sw!($slv_sw, $as, $lookup_sw, $lookup_sw_cld, $lookup_sw_aero)
    return trial_lw, trial_sw
end

function generate_gpu_all_sky_with_aerosols_benchmarks(FT, npts, ::Type{SLVLW}, ::Type{SLVSW}) where {SLVLW, SLVSW}
    context = ClimaComms.context()
    # compute equivalent ncols for DYAMOND resolution
    helems, nlevels, nlev_test, nq = 30, 64, 73, 4
    ncols_dyamond = Int(ceil(helems * helems * 6 * nq * nq * (nlevels / nlev_test)))
    println("\n")
    printstyled(
        "Running DYAMOND all-sky with aerosols benchmark on $(context.device) device with $FT precision\n",
        color = 130,
    )
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
        trial_lw, trial_sw = benchmark_all_sky_with_aerosols(context, SLVLW, SLVSW, FT; ncol = ncols, cldfrac = FT(1))
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
    generate_gpu_all_sky_with_aerosols_benchmarks(FT, 4, NoScatLWRTE, TwoStreamSWRTE)
    generate_gpu_all_sky_with_aerosols_benchmarks(FT, 4, TwoStreamLWRTE, TwoStreamSWRTE)
end
