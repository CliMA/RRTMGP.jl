using Test
using Pkg.Artifacts
using NCDatasets
import JET
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

import Infiltrator

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

# overriding some parameters to match with RRTMGP FORTRAN code

include("reference_files.jl")
include("read_rfmip_clear_sky.jl")
#---------------------------------------------------------------
function clear_sky(
    context,
    ::Type{SLVLW},
    ::Type{SLVSW},
    ::Type{VMR},
    ::Type{FT},
    toler_lw,
    toler_sw;
    exfiltrate = false,
) where {FT, SLVLW, SLVSW, VMR}
    overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016)
    param_set = RRTMGPParameters(FT, overrides)
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)

    lw_file = get_lookup_filename(:gas, :lw) # lw lookup tables for gas optics
    sw_file = get_lookup_filename(:gas, :sw) # sw lookup tables for gas optics
    input_file = get_input_filename(:gas, :lw) # clear-sky atmos state
    # reference data files for comparison
    flux_up_file_lw = get_reference_filename(:gas, :lw, :flux_up)
    flux_dn_file_lw = get_reference_filename(:gas, :lw, :flux_dn)
    flux_up_file_sw = get_reference_filename(:gas, :sw, :flux_up)
    flux_dn_file_sw = get_reference_filename(:gas, :sw, :flux_dn)

    exp_no = 1
    n_gauss_angles = 1

    # reading longwave lookup data
    ds_lw = Dataset(lw_file, "r")
    lookup_lw, idx_gases = LookUpLW(ds_lw, FT, DA)
    close(ds_lw)

    # reading shortwave lookup data
    ds_sw = Dataset(sw_file, "r")
    lookup_sw, idx_gases = LookUpSW(ds_sw, FT, DA)
    close(ds_sw)

    # reading rfmip data to atmospheric state
    ds_lw_in = Dataset(input_file, "r")
    (as, sfc_emis, sfc_alb_direct, cos_zenith, toa_flux) =
        setup_rfmip_as(context, ds_lw_in, idx_gases, exp_no, lookup_lw, FT, VMR, param_set)
    close(ds_lw_in)

    nlay, ncol = AtmosphericStates.get_dims(as)
    nlev = nlay + 1

    # setting up longwave problem
    inc_flux = nothing
    slv_lw = SLVLW(FT, DA, context, param_set, nlay, ncol, sfc_emis, inc_flux)

    # setting up shortwave problem
    sfc_alb_diffuse = DA{FT, 2}(deepcopy(sfc_alb_direct))
    inc_flux_diffuse = nothing
    swbcs = (cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(FT, DA, context, nlay, ncol, swbcs...)
    #--------------------------------------------------
    # initializing RTE solver
    exfiltrate && Infiltrator.@exfiltrate
    solve_lw!(slv_lw, as, lookup_lw)
    if device isa ClimaComms.CPUSingleThreaded
        JET.@test_opt solve_lw!(slv_lw, as, lookup_lw)
        @test (@allocated solve_lw!(slv_lw, as, lookup_lw)) == 0
        @test (@allocated solve_lw!(slv_lw, as, lookup_lw)) ≤ 448
    end

    solve_sw!(slv_sw, as, lookup_sw)
    if device isa ClimaComms.CPUSingleThreaded
        JET.@test_opt solve_sw!(slv_sw, as, lookup_sw)
        @test (@allocated solve_sw!(slv_sw, as, lookup_sw)) == 0
        @test (@allocated solve_sw!(slv_sw, as, lookup_sw)) ≤ 448
    end

    # comparing longwave fluxes with data from RRTMGP FORTRAN code
    flip_ind = nlev:-1:1

    ds_flux_up_lw = Dataset(flux_up_file_lw, "r")
    comp_flux_up_lw = Array(ds_flux_up_lw["rlu"])[flip_ind, :, exp_no]
    close(ds_flux_up_lw)

    ds_flux_dn_lw = Dataset(flux_dn_file_lw, "r")
    comp_flux_dn_lw = Array(ds_flux_dn_lw["rld"])[flip_ind, :, exp_no]
    close(ds_flux_dn_lw)

    comp_flux_net_lw = comp_flux_up_lw .- comp_flux_dn_lw

    flux_up_lw = Array(slv_lw.flux.flux_up)
    flux_dn_lw = Array(slv_lw.flux.flux_dn)
    flux_net_lw = Array(slv_lw.flux.flux_net)

    max_err_flux_up_lw = maximum(abs.(flux_up_lw .- comp_flux_up_lw))
    max_err_flux_dn_lw = maximum(abs.(flux_dn_lw .- comp_flux_dn_lw))
    max_err_flux_net_lw = maximum(abs.(flux_net_lw .- comp_flux_net_lw))

    rel_err_flux_net_lw = abs.(flux_net_lw .- comp_flux_net_lw)

    for gcol in 1:ncol, glev in 1:nlev
        den = abs(comp_flux_net_lw[glev, gcol])
        if den > 10 * eps(FT)
            rel_err_flux_net_lw[glev, gcol] /= den
        end
    end
    max_rel_err_flux_net_lw = maximum(rel_err_flux_net_lw)

    color2 = :cyan
    printstyled("Clear-sky longwave test with ncol = $ncol, nlev = $nlev, Solver = $SLVLW, FT = $FT\n", color = color2)
    printstyled("device = $device\n\n", color = color2)
    println("L∞ error in flux_up           = $max_err_flux_up_lw")
    println("L∞ error in flux_dn           = $max_err_flux_dn_lw")
    println("L∞ error in flux_net          = $max_err_flux_net_lw")
    println("L∞ relative error in flux_net = $(max_rel_err_flux_net_lw * 100) %\n")

    # comparing shortwave fluxes with data from RRTMGP FORTRAN code
    ds_flux_up_sw = Dataset(flux_up_file_sw, "r")
    comp_flux_up_sw = Array(ds_flux_up_sw["rsu"])[flip_ind, :, exp_no]
    close(ds_flux_up_sw)

    ds_flux_dn_sw = Dataset(flux_dn_file_sw, "r")
    comp_flux_dn_sw = Array(ds_flux_dn_sw["rsd"])[flip_ind, :, exp_no]
    close(ds_flux_dn_sw)

    comp_flux_net_sw = comp_flux_up_sw .- comp_flux_dn_sw

    flux_up_sw = Array(slv_sw.flux.flux_up)
    flux_dn_sw = Array(slv_sw.flux.flux_dn)
    flux_net_sw = Array(slv_sw.flux.flux_net)

    # Test if shortwave fluxes are zero if zenith angle is ≥ π/2
    cos_zenith = Array(cos_zenith)
    test_night_cols = true
    nnightcol = 0
    for gcol in 1:ncol
        if cos_zenith[gcol] ≤ 0
            test_flux_up_sw = maximum(abs.(flux_up_sw[:, gcol])) ≈ FT(0)
            test_flux_dn_sw = maximum(abs.(flux_dn_sw[:, gcol])) ≈ FT(0)
            if !(test_flux_up_sw && test_flux_dn_sw)
                test_night_cols = false
            end
            nnightcol += 1
        end
    end
    @test test_night_cols

    max_err_flux_up_sw = maximum(abs.(flux_up_sw .- comp_flux_up_sw))
    max_err_flux_dn_sw = maximum(abs.(flux_dn_sw .- comp_flux_dn_sw))
    max_err_flux_net_sw = maximum(abs.(flux_net_sw .- comp_flux_net_sw))

    rel_err_flux_net_sw = abs.(flux_net_sw .- comp_flux_net_sw)

    for gcol in 1:ncol, glev in 1:nlev
        den = abs(comp_flux_net_sw[glev, gcol])
        if den > 10 * eps(FT)
            rel_err_flux_net_sw[glev, gcol] /= den
        end
    end

    max_rel_err_flux_net_sw = maximum(rel_err_flux_net_sw)

    printstyled("Clear-sky shortwave test with ncol = $ncol, nlev = $nlev, Solver = $SLVSW, FT = $FT\n", color = color2)
    printstyled("device = $device\n\n", color = color2)
    println("L∞ error in flux_up           = $max_err_flux_up_sw")
    println("L∞ error in flux_dn           = $max_err_flux_dn_sw")
    println("L∞ error in flux_net          = $max_err_flux_net_sw")
    println("L∞ relative error in flux_net = $(max_rel_err_flux_net_sw * 100) %\n")

    # Note: The reference results for the longwave solver are generated using a non-scattering solver with rescaling,
    # which differ from the results generated by the TwoStream currently used. Due to the above difference, 
    # the longwave tests using currently have higher error.
    @test max_err_flux_up_lw ≤ toler_lw[FT]
    @test max_err_flux_dn_lw ≤ toler_lw[FT]
    @test max_err_flux_net_lw ≤ toler_lw[FT]

    @test max_err_flux_up_sw ≤ toler_sw[FT]
    @test max_err_flux_dn_sw ≤ toler_sw[FT]
    @test max_err_flux_net_sw ≤ toler_sw[FT]
end
