

function setup_rfmip_as(
    ds_lw_in,
    idx_gases,
    exp_no,
    lookup,
    FT,
    DA,
    max_threads,
)

    nlay = Int(ds_lw_in.dim["layer"])
    ncol = Int(ds_lw_in.dim["site"])
    nlev = nlay + 1
    ngas = lookup.n_gases

    lon = DA{FT,1}(ds_lw_in["lon"][:])
    lat = DA{FT,1}(ds_lw_in["lat"][:])

    sfc_emis = DA{FT,1}(ds_lw_in["surface_emissivity"][:])
    sfc_alb = DA{FT,1}(ds_lw_in["surface_albedo"][:])

    zenith = DA{FT,1}(ds_lw_in["solar_zenith_angle"][:])
    irrad = DA{FT,1}(ds_lw_in["total_solar_irradiance"][:])


    p_lev = ds_lw_in["pres_level"][:, 1]

    lev_ind = p_lev[1] > p_lev[end] ? (1:nlev) : (nlev:-1:1)
    lay_ind = p_lev[1] > p_lev[end] ? (1:nlay) : (nlay:-1:1)

    p_lev = DA{FT,2}(ds_lw_in["pres_level"][:][lev_ind, :])
    t_lev = DA{FT,2}(ds_lw_in["temp_level"][:][lev_ind, :, exp_no])
    p_lay = DA{FT,2}(ds_lw_in["pres_layer"][:][lay_ind, :])
    t_lay = DA{FT,2}(ds_lw_in["temp_layer"][:][lay_ind, :, exp_no])

    t_sfc = DA{FT,1}(ds_lw_in["surface_temperature"][:, exp_no])

    col_dry = DA{FT,2}(undef, nlay, ncol)

    # Reading volume mixing ratios 

    vmr_h2o = DA{FT,2}(ds_lw_in["water_vapor"][:][lay_ind, :, exp_no]) # vmr of H2O and O3
    vmr_o3 = DA{FT,2}(ds_lw_in["ozone"][:][lay_ind, :, exp_no])       # vary with height

    vmrat = zeros(FT, ngas)

    vmrat[idx_gases["co2"]] =
        FT(ds_lw_in["carbon_dioxide_GM"][exp_no]) *
        parse(FT, ds_lw_in["carbon_dioxide_GM"].attrib["units"])

    vmrat[idx_gases["n2o"]] =
        FT(ds_lw_in["nitrous_oxide_GM"][exp_no]) *
        parse(FT, ds_lw_in["nitrous_oxide_GM"].attrib["units"])

    vmrat[idx_gases["co"]] =
        FT(ds_lw_in["carbon_monoxide_GM"][exp_no]) *
        parse(FT, ds_lw_in["carbon_monoxide_GM"].attrib["units"])

    vmrat[idx_gases["ch4"]] =
        FT(ds_lw_in["methane_GM"][exp_no]) *
        parse(FT, ds_lw_in["methane_GM"].attrib["units"])

    vmrat[idx_gases["o2"]] =
        FT(ds_lw_in["oxygen_GM"][exp_no]) *
        parse(FT, ds_lw_in["oxygen_GM"].attrib["units"])

    vmrat[idx_gases["n2"]] =
        FT(ds_lw_in["nitrogen_GM"][exp_no]) *
        parse(FT, ds_lw_in["nitrogen_GM"].attrib["units"])

    vmrat[idx_gases["ccl4"]] =
        FT(ds_lw_in["carbon_tetrachloride_GM"][exp_no]) *
        parse(FT, ds_lw_in["carbon_tetrachloride_GM"].attrib["units"])

    vmrat[idx_gases["cfc11"]] =
        FT(ds_lw_in["cfc11_GM"][exp_no]) *
        parse(FT, ds_lw_in["cfc11_GM"].attrib["units"])

    vmrat[idx_gases["cfc12"]] =
        FT(ds_lw_in["cfc12_GM"][exp_no]) *
        parse(FT, ds_lw_in["cfc12_GM"].attrib["units"])

    vmrat[idx_gases["cfc22"]] =
        FT(ds_lw_in["hcfc22_GM"][exp_no]) *
        parse(FT, ds_lw_in["hcfc22_GM"].attrib["units"])

    vmrat[idx_gases["hfc143a"]] =
        FT(ds_lw_in["hfc143a_GM"][exp_no]) *
        parse(FT, ds_lw_in["hfc143a_GM"].attrib["units"])

    vmrat[idx_gases["hfc125"]] =
        FT(ds_lw_in["hfc125_GM"][exp_no]) *
        parse(FT, ds_lw_in["hfc125_GM"].attrib["units"])

    vmrat[idx_gases["hfc23"]] =
        FT(ds_lw_in["hfc23_GM"][exp_no]) *
        parse(FT, ds_lw_in["hfc23_GM"].attrib["units"])

    vmrat[idx_gases["hfc32"]] =
        FT(ds_lw_in["hfc32_GM"][exp_no]) *
        parse(FT, ds_lw_in["hfc32_GM"].attrib["units"])

    vmrat[idx_gases["hfc134a"]] =
        FT(ds_lw_in["hfc134a_GM"][exp_no]) *
        parse(FT, ds_lw_in["hfc134a_GM"].attrib["units"])

    vmrat[idx_gases["cf4"]] =
        FT(ds_lw_in["cf4_GM"][exp_no]) *
        parse(FT, ds_lw_in["hfc23_GM"].attrib["units"])

    #    vmr[idx_gases["no2"]] = FT(ds_lw_in["no2_GM"][exp_no]) *                # missing from input file
    #                                         parse(FT, ds_lw_in["hfc32_GM"].attrib["units"])

    compute_col_dry!(p_lev, t_lay, col_dry, param_set, vmr_h2o, lat)

    vmr = Vmr(vmr_h2o, vmr_o3, DA(vmrat))
    #------------------
    return AtmosphericState{FT,DA{FT,1},DA{FT,2},typeof(vmr),Int}(
        lon,
        lat,
        sfc_emis,
        sfc_alb,
        zenith,
        irrad,
        p_lay,
        p_lev,
        t_lay,
        t_lev,
        t_sfc,
        col_dry,
        vmr,
        nlay,
        ncol,
        ngas,
    )
end
