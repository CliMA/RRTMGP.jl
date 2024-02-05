function ncol_ds_all_sky(use_lut)
    lut_pade = use_lut ? :lut : :pade
    flux_file = get_ref_filename(:comparison, :cloudysky, lut_pade = lut_pade)
    ds_comp = Dataset(flux_file, "r")
    return size(Array(ds_comp["lw_flux_up"]), 1)
end

function setup_allsky_as(
    context,
    ds_in,
    idx_gases,
    lkp_lw,
    lkp_sw,
    lkp_lw_cld,
    lkp_sw_cld,
    cldfrac,
    use_lut,
    ncol,
    FT,
    max_threads,
    param_set,
)

    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)
    deg2rad = FT(π) / FT(180)
    nlay = Int(ds_in.dim["lay"])
    #ncol = Int(ds_in.dim["col"]) # col#1 repeated 128 times, per RRTMGP example
    nlev = nlay + 1
    ngas = lkp_lw.n_gases
    nbnd_lw = lkp_lw.n_bnd
    nbnd_sw = lkp_sw.n_bnd
    ngpt_lw = lkp_lw.n_gpt
    ngpt_sw = lkp_sw.n_gpt
    #---no lat / long information for this
    lon = nothing
    lat = nothing
    # The example only reads the first column and 
    # replicates it ncol times.

    sfc_emis = DA{FT, 2}(undef, nbnd_lw, ncol)
    sfc_alb_direct = DA{FT, 2}(undef, nbnd_sw, ncol)
    sfc_alb_diffuse = DA{FT, 2}(undef, nbnd_sw, ncol)
    cos_zenith = DA{FT, 1}(undef, ncol)
    irrad = DA{FT, 1}(undef, ncol)
    # these values are taken from the example
    sfc_emis .= FT(0.98)
    sfc_alb_direct .= FT(0.06)
    sfc_alb_diffuse .= FT(0.06)
    cos_zenith .= FT(0.86)
    irrad .= FT(lkp_sw.solar_src_tot)


    p_lev = Array{FT}(reshape(Array(ds_in["p_lev"])[1, :], nlev, 1))

    bot_at_1 = p_lev[1, 1] > p_lev[end, 1]
    lev_ind = bot_at_1 ? (1:nlev) : (nlev:-1:1)
    lay_ind = bot_at_1 ? (1:nlay) : (nlay:-1:1)

    p_lev = Array{FT}(reshape(Array(ds_in["p_lev"])[1, lev_ind], nlev, 1))
    p_lay = Array{FT}(reshape(Array(ds_in["p_lay"])[1, lay_ind], nlay, 1))
    t_lev = Array{FT}(reshape(Array(ds_in["t_lev"])[1, lev_ind], nlev, 1))
    t_lay = Array{FT}(reshape(Array(ds_in["t_lay"])[1, lay_ind], nlay, 1))
    t_sfc = Array{FT}(reshape([t_lev[1, 1]], 1))

    p_lev = repeat(p_lev, 1, ncol)
    p_lay = repeat(p_lay, 1, ncol)
    t_lev = repeat(t_lev, 1, ncol)
    t_lay = repeat(t_lay, 1, ncol)
    t_sfc = repeat(t_sfc, ncol)
    #col_dry = DA{FT,2}(transpose(Array(ds_in["col_dry"])[:, lay_ind]))
    #col_dry from the dataset not used in the FORTRAN RRTMGP example

    # Reading volume mixing ratios 
    vmrat = zeros(FT, ngas, nlay, ncol)

    vmrat[idx_gases["h2o"], :, 1] .= Array{FT}(Array(ds_in["vmr_h2o"])[1, lay_ind])
    vmrat[idx_gases["o3"], :, 1] .= Array{FT}(Array(ds_in["vmr_o3"])[1, lay_ind])
    vmrat[idx_gases["co2"], :, 1] .= Array{FT}(Array(ds_in["vmr_co2"])[1, lay_ind])
    vmrat[idx_gases["n2o"], :, 1] .= Array{FT}(Array(ds_in["vmr_n2o"])[1, lay_ind])
    vmrat[idx_gases["co"], :, 1] .= Array{FT}(Array(ds_in["vmr_co"])[1, lay_ind])
    vmrat[idx_gases["ch4"], :, 1] .= Array{FT}(Array(ds_in["vmr_ch4"])[1, lay_ind])
    vmrat[idx_gases["o2"], :, 1] .= Array{FT}(Array(ds_in["vmr_o2"])[1, lay_ind])
    vmrat[idx_gases["n2"], :, 1] .= Array{FT}(Array(ds_in["vmr_n2"])[1, lay_ind])

    for icol in 2:ncol
        vmrat[:, :, icol] .= vmrat[:, :, 1]
    end
    vmr = Vmr(DA(vmrat))
    col_dry = DA{FT, 2}(undef, nlay, ncol)
    vmr_h2o = view(vmr.vmr, idx_gases["h2o"], :, :)

    cld_frac = zeros(FT, nlay, ncol)
    cld_mask_lw = zeros(Bool, nlay, ncol)
    cld_mask_sw = zeros(Bool, nlay, ncol)
    cld_r_eff_liq = zeros(FT, nlay, ncol)
    cld_r_eff_ice = zeros(FT, nlay, ncol)
    cld_path_liq = zeros(FT, nlay, ncol)
    cld_path_ice = zeros(FT, nlay, ncol)

    if use_lut
        r_eff_liq = (lkp_lw_cld.radliq_lwr + lkp_lw_cld.radliq_upr) / FT(2)
        r_eff_ice = (lkp_lw_cld.radice_lwr + lkp_lw_cld.radice_upr) / FT(2)
    else
        pade_sizreg_extliq = Array(lkp_lw_cld.pade_sizreg_extliq) # TODO temp fix to avoid scalar indexing on GPU
        pade_sizreg_extice = Array(lkp_lw_cld.pade_sizreg_extice)
        r_eff_liq = sum(extrema(pade_sizreg_extliq)) / FT(2)
        r_eff_ice = sum(extrema(pade_sizreg_extice)) / FT(2)
    end
    # Restrict clouds to troposphere (> 100 hPa = 100*100 Pa)
    # and not very close to the ground (< 900 hPa), and
    # put them in 2/3 of the columns since that's roughly the
    # total cloudiness of earth
    ncol_ds = ncol_ds_all_sky(use_lut)
    for icol in 1:ncol, ilay in 1:nlay
        icol_ds = icol % ncol_ds == 0 ? ncol_ds : icol % ncol_ds
        if p_lay[ilay, icol] > FT(10000) && p_lay[ilay, icol] < FT(90000) && icol_ds % 3 ≠ 0
            cld_frac[ilay, icol] = cldfrac
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

    compute_col_gas!(context, p_lev, col_dry, param_set, vmr_h2o, lat) # the example skips lat based gravity calculation

    t_sfc = DA(t_sfc)

    cld_frac = DA(cld_frac)
    cld_mask_lw = DA(cld_mask_lw)
    cld_mask_sw = DA(cld_mask_sw)
    cld_r_eff_liq = DA(cld_r_eff_liq)
    cld_r_eff_ice = DA(cld_r_eff_ice)
    cld_path_liq = DA(cld_path_liq)
    cld_path_ice = DA(cld_path_ice)
    ice_rgh = 2 # medium ice roughness

    return (
        AtmosphericState(
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
            cld_frac,
            cld_mask_lw,
            cld_mask_sw,
            MaxRandomOverlap(),
            ice_rgh,
            nlay,
            ncol,
            ngas,
        ),
        sfc_emis,
        sfc_alb_direct,
        sfc_alb_diffuse,
        cos_zenith,
        irrad,
        bot_at_1,
    )
end

_orient_data(data, bot_at_1) = bot_at_1 ? data : reverse(data, dims = 1)

function load_comparison_data(use_lut, bot_at_1, ncol)
    lut_pade = use_lut ? :lut : :pade
    flux_file = get_ref_filename(:comparison, :cloudysky, lut_pade = lut_pade) # flux files for comparison (LUT)
    ds_comp = Dataset(flux_file, "r")
    ncol_ds = size(Array(ds_comp["lw_flux_up"]), 1)
    nrepeat = cld(ncol, ncol_ds)
    comp_flux_up_lw = repeat(_orient_data(transpose(Array(ds_comp["lw_flux_up"])), bot_at_1), 1, nrepeat)
    comp_flux_dn_lw = repeat(_orient_data(transpose(Array(ds_comp["lw_flux_dn"])), bot_at_1), 1, nrepeat)
    comp_flux_up_sw = repeat(_orient_data(transpose(Array(ds_comp["sw_flux_up"])), bot_at_1), 1, nrepeat)
    comp_flux_dn_sw = repeat(_orient_data(transpose(Array(ds_comp["sw_flux_dn"])), bot_at_1), 1, nrepeat)
    close(ds_comp)
    nlev, ncol_ds = size(comp_flux_up_lw)
    return comp_flux_up_lw[:, 1:ncol],
    comp_flux_dn_lw[:, 1:ncol],
    comp_flux_up_sw[:, 1:ncol],
    comp_flux_dn_sw[:, 1:ncol]
end
