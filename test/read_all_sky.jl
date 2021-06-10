
function setup_allsky_as(
    ds_in,
    idx_gases,
    lkp_lw,
    lkp_sw,
    lkp_lw_cld,
    lkp_sw_cld,
    use_lut,
    ncol,
    FT,
    DA,
    max_threads,
)

    deg2rad = FT(π) / FT(180)
    nlay = Int(ds_in.dim["lay"])
    #ncol = Int(ds_in.dim["col"]) # col#1 repeated 128 times, per RRTMGP example
    nlev = nlay + 1
    ngas = lkp_lw.n_gases
    nbnd_lw = lkp_lw.n_bnd
    nbnd_sw = lkp_sw.n_bnd

    #---no lat / long information for this
    lon = nothing
    lat = nothing
    # The example only reads the first column and 
    # replicates it ncol times.

    sfc_emis = DA{FT,2}(undef, nbnd_lw, ncol)
    sfc_alb_direct = DA{FT,2}(undef, nbnd_sw, ncol)
    sfc_alb_diffuse = DA{FT,2}(undef, nbnd_sw, ncol)
    zenith = DA{FT,1}(undef, ncol)
    irrad = DA{FT,1}(undef, ncol)
    # these values are taken from the example
    sfc_emis .= FT(0.98)
    sfc_alb_direct .= FT(0.06)
    sfc_alb_diffuse .= FT(0.06)
    zenith .= acos(FT(0.86))
    irrad .= FT(lkp_sw.solar_src_tot)


    p_lev = Array{FT}(reshape(ds_in["p_lev"][:][1, :], nlev, 1))

    bot_at_1 = p_lev[1, 1] > p_lev[end, 1]
    lev_ind = bot_at_1 ? (1:nlev) : (nlev:-1:1)
    lay_ind = bot_at_1 ? (1:nlay) : (nlay:-1:1)

    p_lev = Array{FT}(reshape(ds_in["p_lev"][:][1, lev_ind], nlev, 1))
    p_lay = Array{FT}(reshape(ds_in["p_lay"][:][1, lay_ind], nlay, 1))
    t_lev = Array{FT}(reshape(ds_in["t_lev"][:][1, lev_ind], nlev, 1))
    t_lay = Array{FT}(reshape(ds_in["t_lay"][:][1, lay_ind], nlay, 1))
    t_sfc = Array{FT}(reshape([t_lev[1, 1]], 1))

    p_lev = repeat(p_lev, 1, ncol)
    p_lay = repeat(p_lay, 1, ncol)
    t_lev = repeat(t_lev, 1, ncol)
    t_lay = repeat(t_lay, 1, ncol)
    t_sfc = repeat(t_sfc, ncol)
    #col_dry = DA{FT,2}(transpose(ds_in["col_dry"][:][:, lay_ind]))
    #col_dry from the dataset not used in the FORTRAN RRTMGP example

    # Reading volume mixing ratios 
    vmrat = zeros(FT, nlay, ncol, ngas)

    vmrat[:, 1, idx_gases["h2o"]] .= Array{FT}(ds_in["vmr_h2o"][:][1, lay_ind])
    vmrat[:, 1, idx_gases["o3"]] .= Array{FT}(ds_in["vmr_o3"][:][1, lay_ind])
    vmrat[:, 1, idx_gases["co2"]] .= Array{FT}(ds_in["vmr_co2"][:][1, lay_ind])
    vmrat[:, 1, idx_gases["n2o"]] .= Array{FT}(ds_in["vmr_n2o"][:][1, lay_ind])
    vmrat[:, 1, idx_gases["co"]] .= Array{FT}(ds_in["vmr_co"][:][1, lay_ind])
    vmrat[:, 1, idx_gases["ch4"]] .= Array{FT}(ds_in["vmr_ch4"][:][1, lay_ind])
    vmrat[:, 1, idx_gases["o2"]] .= Array{FT}(ds_in["vmr_o2"][:][1, lay_ind])
    vmrat[:, 1, idx_gases["n2"]] .= Array{FT}(ds_in["vmr_n2"][:][1, lay_ind])

    for icol = 2:ncol
        vmrat[:, icol, :] .= vmrat[:, 1, :]
    end
    vmr = Vmr(DA(vmrat))
    col_dry = DA{FT,2}(undef, nlay, ncol)
    vmr_h2o = DA(vmrat[:, :, idx_gases["h2o"]])
    compute_col_dry!(p_lev, t_lay, col_dry, param_set, vmr_h2o, lat)

    cld_mask = zeros(Bool, nlay, ncol)
    cld_r_eff_liq = zeros(FT, nlay, ncol)
    cld_r_eff_ice = zeros(FT, nlay, ncol)
    cld_path_liq = zeros(FT, nlay, ncol)
    cld_path_ice = zeros(FT, nlay, ncol)

    if use_lut
        r_eff_liq = (lkp_lw_cld.radliq_lwr + lkp_lw_cld.radliq_upr) / FT(2)
        r_eff_ice = (lkp_lw_cld.radice_lwr + lkp_lw_cld.radice_upr) / FT(2)
    else
        r_eff_liq = sum(extrema(lkp_lw_cld.pade_sizreg_extliq)) / FT(2)
        r_eff_ice = sum(extrema(lkp_lw_cld.pade_sizreg_extice)) / FT(2)
    end
    # Restrict clouds to troposphere (> 100 hPa = 100*100 Pa)
    # and not very close to the ground (< 900 hPa), and
    # put them in 2/3 of the columns since that's roughly the
    # total cloudiness of earth
    for icol = 1:ncol, ilay = 1:nlay
        if p_lay[ilay, icol] > FT(10000) &&
           p_lay[ilay, icol] < FT(90000) &&
           icol % 3 ≠ 0
            cld_mask[ilay, icol] = true
            if t_lay[ilay, icol] > FT(263)
                cld_path_liq[ilay, icol] = FT(10)
                cld_r_eff_liq[ilay, icol] = r_eff_liq
            end
            if t_lay[ilay, icol] < FT(273)
                cld_path_ice[ilay, icol] = FT(10)
                cld_r_eff_ice[ilay, icol] = r_eff_ice
            end
        end
    end


    p_lay = DA(p_lay)
    p_lev = DA(p_lev)
    t_lay = DA(t_lay)
    t_lev = DA(t_lev)
    t_sfc = DA(t_sfc)

    cld_mask = DA(cld_mask)
    cld_r_eff_liq = DA(cld_r_eff_liq)
    cld_r_eff_ice = DA(cld_r_eff_ice)
    cld_path_liq = DA(cld_path_liq)
    cld_path_ice = DA(cld_path_ice)
    ice_rgh = 2 # medium ice roughness

    CLDP = typeof(cld_r_eff_liq)
    CLDM = typeof(cld_mask)

    return (
        AtmosphericState{
            FT,
            DA{FT,1},
            typeof(lat),
            DA{FT,2},
            CLDP,
            CLDM,
            typeof(vmr),
            Int,
        }(
            lon,
            lat,
            p_lay,
            p_lev,
            t_lay,
            t_lev,
            t_sfc,
            col_dry,
            vmr,
            cld_r_eff_liq,
            cld_r_eff_ice,
            cld_path_liq,
            cld_path_ice,
            cld_mask,
            ice_rgh,
            nlay,
            ncol,
            ngas,
        ),
        sfc_emis,
        sfc_alb_direct,
        sfc_alb_diffuse,
        zenith,
        irrad,
        bot_at_1,
    )
end
