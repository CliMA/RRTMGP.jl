function ncol_ds_clear_sky()
    flux_file = get_reference_filename(:gas, :lw, :flux_up)
    ds_comp = Dataset(flux_file, "r")
    return size(Array(ds_comp["rlu"]), 2)
end

function setup_rfmip_as(
    context,
    ds_lw_in,
    idx_gases,
    expt_no,
    lookup_lw,
    ncol,
    ::Type{FT},
    ::Type{VMR},
    param_set,
) where {FT <: AbstractFloat, VMR}
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    FTA1D = DA{FT, 1}
    FTA2D = DA{FT, 2}

    deg2rad = FT(π) / FT(180)
    nlay = Int(ds_lw_in.dim["layer"])
    ncol_ds = ncol_ds_clear_sky()
    nrepeat = cld(ncol, ncol_ds)
    nlev = nlay + 1
    ngas = LookUpTables.get_n_gases(lookup_lw)
    nbnd_lw = LookUpTables.get_n_bnd(lookup_lw)
    lon = FTA1D(repeat(Array(ds_lw_in["lon"]), nrepeat)[1:ncol])
    lat = FTA1D(repeat(Array(ds_lw_in["lat"]), nrepeat)[1:ncol])

    lon = nothing # This example skips latitude dependent gravity computation
    lat = nothing # to be consistent with the FORTRAN RRTMGP test case.

    # all bands use same emissivity
    sfc_emis = repeat(reshape(Array{FT}(Array(ds_lw_in["surface_emissivity"])), 1, :), nbnd_lw, 1)
    sfc_emis = FTA2D(repeat(sfc_emis, 1, nrepeat)[:, 1:ncol])
    # all bands use same albedo
    sfc_alb = repeat(reshape(Array{FT}(Array(ds_lw_in["surface_albedo"])), 1, :), nbnd_lw, 1)
    sfc_alb = FTA2D(repeat(sfc_alb, 1, nrepeat)[:, 1:ncol])
    #--------------------------------------------------------------
    zenith = Array{FT, 1}(deg2rad .* Array(ds_lw_in["solar_zenith_angle"]))
    irrad = Array{FT, 1}(Array(ds_lw_in["total_solar_irradiance"]))

    cos_zenith = FTA1D(repeat(cos.(zenith), nrepeat)[1:ncol])
    irrad = FTA1D(repeat(irrad, nrepeat)[1:ncol])
    #--------------------------------------------------------------

    p_lev = Array(ds_lw_in["pres_level"])
    bot_at_1 = p_lev[1, 1] > p_lev[end, 1]

    lev_ind = bot_at_1 ? (1:nlev) : (nlev:-1:1)
    lay_ind = bot_at_1 ? (1:nlay) : (nlay:-1:1)

    p_lev[lev_ind[end], :] .= lookup_lw.p_ref_min

    p_lev = p_lev[lev_ind, :]
    p_lay = Array(ds_lw_in["pres_layer"])[lay_ind, :]
    t_lev = Array(ds_lw_in["temp_level"])[lev_ind, :, expt_no]
    t_lay = Array(ds_lw_in["temp_layer"])[lay_ind, :, expt_no]

    p_lev = FTA2D(repeat(p_lev, 1, nrepeat)[:, 1:ncol])
    p_lay = FTA2D(repeat(p_lay, 1, nrepeat)[:, 1:ncol])
    t_lev = FTA2D(repeat(t_lev, 1, nrepeat)[:, 1:ncol])
    t_lay = FTA2D(repeat(t_lay, 1, nrepeat)[:, 1:ncol])

    t_sfc = FTA1D(repeat(ds_lw_in["surface_temperature"][:, expt_no], nrepeat)[1:ncol])
    col_dry = FTA2D(undef, nlay, ncol)
    rel_hum = FTA2D(undef, nlay, ncol)

    # Reading volume mixing ratios 

    vmr_h2o = (Array(ds_lw_in["water_vapor"])[lay_ind, :, expt_no]) # vmr of H2O and O3
    vmr_o3 = (Array(ds_lw_in["ozone"])[lay_ind, :, expt_no])       # vary with height
    vmr_h2o = FTA2D(repeat(vmr_h2o, 1, nrepeat)[:, 1:ncol])
    vmr_o3 = FTA2D(repeat(vmr_o3, 1, nrepeat)[:, 1:ncol])

    vmrat = zeros(FT, ngas)

    vmrat[idx_gases["co2"]] =
        FT(ds_lw_in["carbon_dioxide_GM"][expt_no]) * parse(FT, ds_lw_in["carbon_dioxide_GM"].attrib["units"])

    vmrat[idx_gases["n2o"]] =
        FT(ds_lw_in["nitrous_oxide_GM"][expt_no]) * parse(FT, ds_lw_in["nitrous_oxide_GM"].attrib["units"])

    vmrat[idx_gases["co"]] =
        FT(ds_lw_in["carbon_monoxide_GM"][expt_no]) * parse(FT, ds_lw_in["carbon_monoxide_GM"].attrib["units"])

    vmrat[idx_gases["ch4"]] = FT(ds_lw_in["methane_GM"][expt_no]) * parse(FT, ds_lw_in["methane_GM"].attrib["units"])

    vmrat[idx_gases["o2"]] = FT(ds_lw_in["oxygen_GM"][expt_no]) * parse(FT, ds_lw_in["oxygen_GM"].attrib["units"])

    vmrat[idx_gases["n2"]] = FT(ds_lw_in["nitrogen_GM"][expt_no]) * parse(FT, ds_lw_in["nitrogen_GM"].attrib["units"])

    vmrat[idx_gases["ccl4"]] =
        FT(ds_lw_in["carbon_tetrachloride_GM"][expt_no]) *
        parse(FT, ds_lw_in["carbon_tetrachloride_GM"].attrib["units"])

    vmrat[idx_gases["cfc11"]] = FT(ds_lw_in["cfc11_GM"][expt_no]) * parse(FT, ds_lw_in["cfc11_GM"].attrib["units"])

    vmrat[idx_gases["cfc12"]] = FT(ds_lw_in["cfc12_GM"][expt_no]) * parse(FT, ds_lw_in["cfc12_GM"].attrib["units"])

    vmrat[idx_gases["cfc22"]] = FT(ds_lw_in["hcfc22_GM"][expt_no]) * parse(FT, ds_lw_in["hcfc22_GM"].attrib["units"])

    vmrat[idx_gases["hfc143a"]] =
        FT(ds_lw_in["hfc143a_GM"][expt_no]) * parse(FT, ds_lw_in["hfc143a_GM"].attrib["units"])

    vmrat[idx_gases["hfc125"]] = FT(ds_lw_in["hfc125_GM"][expt_no]) * parse(FT, ds_lw_in["hfc125_GM"].attrib["units"])

    vmrat[idx_gases["hfc23"]] = FT(ds_lw_in["hfc23_GM"][expt_no]) * parse(FT, ds_lw_in["hfc23_GM"].attrib["units"])

    vmrat[idx_gases["hfc32"]] = FT(ds_lw_in["hfc32_GM"][expt_no]) * parse(FT, ds_lw_in["hfc32_GM"].attrib["units"])

    vmrat[idx_gases["hfc134a"]] =
        FT(ds_lw_in["hfc134a_GM"][expt_no]) * parse(FT, ds_lw_in["hfc134a_GM"].attrib["units"])

    vmrat[idx_gases["cf4"]] = FT(ds_lw_in["cf4_GM"][expt_no]) * parse(FT, ds_lw_in["hfc23_GM"].attrib["units"])

    #    vmr[idx_gases["no2"]] = FT(ds_lw_in["no2_GM"][expt_no]) *                # missing from input file
    #                                         parse(FT, ds_lw_in["hfc32_GM"].attrib["units"])

    # This example skips latitude dependent gravity compution to be consistent with the
    # FORTRAN RRTMGP test case.
    device = ClimaComms.device(context)
    compute_col_gas!(device, p_lev, col_dry, param_set, vmr_h2o, lat) # the example skips lat based gravity calculation
    compute_relative_humidity!(device, rel_hum, p_lay, t_lay, param_set, vmr_h2o) # compute relative humidity

    layerdata = similar(p_lay, 4, nlay, ncol)
    layerdata[1, :, :] .= col_dry
    layerdata[2, :, :] .= p_lay
    layerdata[3, :, :] .= t_lay
    layerdata[4, :, :] .= rel_hum

    vmr = VMR(vmr_h2o, vmr_o3, FTA1D(vmrat))
    #------------------
    return (
        AtmosphericState(lon, lat, layerdata, p_lev, t_lev, t_sfc, vmr, nothing, nothing),
        sfc_emis,
        sfc_alb,
        cos_zenith,
        irrad,
        bot_at_1,
    )
end

_orient_data(data, bot_at_1) = bot_at_1 ? data : reverse(data, dims = 1)

function load_comparison_data(expt_no, bot_at_1, ncol)
    flux_up_file_lw = get_reference_filename(:gas, :lw, :flux_up)
    flux_dn_file_lw = get_reference_filename(:gas, :lw, :flux_dn)
    flux_up_file_sw = get_reference_filename(:gas, :sw, :flux_up)
    flux_dn_file_sw = get_reference_filename(:gas, :sw, :flux_dn)

    ds_comp_lw_up = Dataset(flux_up_file_lw, "r")
    ds_comp_lw_dn = Dataset(flux_dn_file_lw, "r")
    ds_comp_sw_up = Dataset(flux_up_file_sw, "r")
    ds_comp_sw_dn = Dataset(flux_dn_file_sw, "r")

    ncol_ds = size(Array(ds_comp_lw_up["rlu"]), 2)
    nrepeat = cld(ncol, ncol_ds)

    comp_flux_up_lw = repeat(_orient_data(Array(ds_comp_lw_up["rlu"][:, :, expt_no]), bot_at_1), 1, nrepeat)
    comp_flux_dn_lw = repeat(_orient_data(Array(ds_comp_lw_dn["rld"][:, :, expt_no]), bot_at_1), 1, nrepeat)
    comp_flux_up_sw = repeat(_orient_data(Array(ds_comp_sw_up["rsu"][:, :, expt_no]), bot_at_1), 1, nrepeat)
    comp_flux_dn_sw = repeat(_orient_data(Array(ds_comp_sw_dn["rsd"][:, :, expt_no]), bot_at_1), 1, nrepeat)

    close(ds_comp_lw_up)
    close(ds_comp_lw_dn)
    close(ds_comp_sw_up)
    close(ds_comp_sw_dn)

    nlev, ncol_ds = size(comp_flux_up_lw)

    return Array(transpose(comp_flux_up_lw[:, 1:ncol])),
    Array(transpose(comp_flux_dn_lw[:, 1:ncol])),
    Array(transpose(comp_flux_up_sw[:, 1:ncol])),
    Array(transpose(comp_flux_dn_sw[:, 1:ncol]))
end
