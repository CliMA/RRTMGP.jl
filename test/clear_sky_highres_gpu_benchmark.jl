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
using RRTMGP.VolumeMixingRatios
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

if !@isdefined(_resample_benchmark_as)
    function _resample_benchmark_as(
        context,
        as::AtmosphericStates.AtmosphericState,
        target_nlay,
        ::Type{FT},
    ) where {FT}
        nlay_orig, _ = AtmosphericStates.get_dims(as)
        nlay_orig == target_nlay && return as
        device = ClimaComms.device(context)
        DA = ClimaComms.array_type(device)
        nlev_orig = nlay_orig + 1
        target_nlev = target_nlay + 1

        lay_idx = [
            clamp(round(Int, k * nlay_orig / target_nlay), 1, nlay_orig) for
            k in 1:target_nlay
        ]
        lev_idx = [
            clamp(round(Int, l * nlev_orig / target_nlev), 1, nlev_orig) for
            l in 1:target_nlev
        ]

        new_layerdata = DA{FT}(Array(as.layerdata)[:, lay_idx, :])
        new_p_lev = DA{FT}(Array(as.p_lev)[lev_idx, :])
        new_t_lev = DA{FT}(Array(as.t_lev)[lev_idx, :])

        vmr = as.vmr
        new_vmr = if vmr isa VolumeMixingRatios.VmrGM
            VolumeMixingRatios.VmrGM(
                DA{FT}(Array(vmr.vmr_h2o)[lay_idx, :]),
                DA{FT}(Array(vmr.vmr_o3)[lay_idx, :]),
                vmr.vmr,
            )
        elseif vmr isa VolumeMixingRatios.Vmr
            VolumeMixingRatios.Vmr(DA{FT}(Array(vmr.vmr)[:, lay_idx, :]))
        else
            vmr
        end

        cld = as.cloud_state
        new_cld = if cld !== nothing
            AtmosphericStates.CloudState(
                DA{FT}(Array(cld.cld_r_eff_liq)[lay_idx, :]),
                DA{FT}(Array(cld.cld_r_eff_ice)[lay_idx, :]),
                DA{FT}(Array(cld.cld_path_liq)[lay_idx, :]),
                DA{FT}(Array(cld.cld_path_ice)[lay_idx, :]),
                DA{FT}(Array(cld.cld_frac)[lay_idx, :]),
                cld.cld_cover_sw,
                cld.cld_cover_lw,
                DA{Bool}(Array(cld.mask_lw)[lay_idx, :]),
                DA{Bool}(Array(cld.mask_sw)[lay_idx, :]),
                cld.mask_type,
                cld.ice_rgh,
            )
        else
            nothing
        end

        aer = as.aerosol_state
        new_aer = if aer !== nothing
            AtmosphericStates.AerosolState(
                aer.aod_sw_ext,
                aer.aod_sw_sca,
                DA{Bool}(Array(aer.aero_mask)[lay_idx, :]),
                DA{FT}(Array(aer.aero_size)[:, lay_idx, :]),
                DA{FT}(Array(aer.aero_mass)[:, lay_idx, :]),
            )
        else
            nothing
        end

        return AtmosphericStates.AtmosphericState(
            as.lon,
            as.lat,
            new_layerdata,
            new_p_lev,
            new_t_lev,
            as.t_sfc,
            new_vmr,
            new_cld,
            new_aer,
        )
    end
end

function benchmark_clear_sky(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{VMR},
    ::Type{FT};
    ncol = 100,
    nlay = nothing,
) where {FT <: AbstractFloat, SLVLW, SLVSW, VMR}
    overrides =
        (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
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
        setup_clear_sky_as(
            context,
            ds_lw_in,
            idx_gases,
            expt_no,
            lookup_lw,
            ncol,
            FT,
            VMR,
            param_set,
        )
    close(ds_lw_in)

    if nlay !== nothing
        as = _resample_benchmark_as(context, as, nlay, FT)
    end
    nlay, _ = AtmosphericStates.get_dims(as)
    nlev = nlay + 1
    grid_params = RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)


    # setting up longwave problem
    inc_flux = nothing
    slv_lw = SLVLW(grid_params; params = param_set, sfc_emis, inc_flux)

    # setting up shortwave problem
    sfc_alb_diffuse = FTA2D(deepcopy(sfc_alb_direct))
    inc_flux_diffuse = nothing
    swbcs = (;
        cos_zenith,
        toa_flux,
        sfc_alb_direct,
        inc_flux_diffuse,
        sfc_alb_diffuse,
    )
    slv_sw = SLVSW(grid_params; swbcs...)
    #------calling solvers
    metric_scaling = DA(one.(as.p_lev))
    solve_lw!(slv_lw, as, lookup_lw, nothing, nothing, metric_scaling)
    trial_lw = if device isa ClimaComms.CUDADevice
        @benchmark CUDA.@sync solve_lw!(
            $slv_lw,
            $as,
            $lookup_lw,
            nothing,
            nothing,
            $metric_scaling,
        )
    else
        @benchmark solve_lw!(
            $slv_lw,
            $as,
            $lookup_lw,
            nothing,
            nothing,
            $metric_scaling,
        )
    end

    solve_sw!(slv_sw, as, lookup_sw, nothing, nothing, metric_scaling)
    trial_sw = if device isa ClimaComms.CUDADevice
        @benchmark CUDA.@sync solve_sw!(
            $slv_sw,
            $as,
            $lookup_sw,
            nothing,
            nothing,
            $metric_scaling,
        )
    else
        @benchmark solve_sw!(
            $slv_sw,
            $as,
            $lookup_sw,
            nothing,
            nothing,
            $metric_scaling,
        )
    end
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
    # Column count matching the memory footprint of a high-resolution run: a
    # cubed sphere of 30² elements × 6 panels × 4² quadrature points (86,400
    # columns at roughly 80-km spacing) on 64 levels, rescaled to this test's
    # level count so the total degrees of freedom stay the same.
    helems, nlevels, nlev_test, nq = 30, 64, 61, 4
    ncols_highres =
        Int(ceil(helems * helems * 6 * nq * nq * (nlevels / nlev_test)))
    println("\n")
    printstyled(
        "Running clear-sky benchmark on $(context.device) device with $FT precision\n",
        color = 130,
    )
    printstyled(
        "Longwave solver = $SLVLW; Shortwave solver = $SLVSW\n",
        color = 130,
    )
    printstyled(
        "==============|====================================|==================================\n",
        color = 130,
    )
    printstyled(
        "  ncols       |   median time for longwave solver  | median time for shortwave solver \n",
        color = :green,
    )
    printstyled(
        "==============|====================================|==================================\n",
        color = 130,
    )
    for pts in 1:npts
        ncols = unsafe_trunc(Int, cld(ncols_highres, 2^(pts - 1)))
        ndof = ncols * nlev_test
        sz_per_fld_gb = ndof * sizeof(FT) / 1024 / 1024 / 1024
        trial_lw, trial_sw =
            benchmark_clear_sky(context, SLVLW, SLVSW, VMR, FT; ncol = ncols)
        Printf.@printf(
            "%10i    |           %25s|       %25s \n",
            ncols,
            Statistics.median(trial_lw),
            Statistics.median(trial_sw)
        )
    end
    printstyled(
        "==============|====================================|==================================\n",
        color = 130,
    )
    return nothing
end

# Run the resolution sweep only when executed as a script;
# perf/benchmark_ratchet.jl includes this file for `benchmark_clear_sky`.
if abspath(PROGRAM_FILE) == @__FILE__
    for FT in (Float32, Float64)
        generate_gpu_clear_sky_benchmarks(
            FT,
            4,
            NoScatLWRTE,
            TwoStreamSWRTE,
            VmrGM,
        )
        generate_gpu_clear_sky_benchmarks(
            FT,
            4,
            TwoStreamLWRTE,
            TwoStreamSWRTE,
            VmrGM,
        )
    end
end
