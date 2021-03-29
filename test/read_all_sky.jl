
function setup_allsky_as(ds_lw_in, idx_gases, lookup, FT, DA, max_threads)

    deg2rad = FT(π) / FT(180)
    nlay = Int(ds_lw_in.dim["lay"])
    ncol = Int(ds_lw_in.dim["col"])
    nlev = nlay + 1
    ngas = lookup.n_gases

    #---no lat / long information for this
    lon = nothing
    lat = nothing

    p_lev = (ds_lw_in["p_lev"][:])
    lev_ind = p_lev[1, 1] > p_lev[1, end] ? (1:nlev) : (nlev:-1:1)
    lay_ind = p_lev[1, 1] > p_lev[1, end] ? (1:nlay) : (nlay:-1:1)

    p_lev = DA{FT,2}(transpose(p_lev[:, lev_ind]))
    p_lay = DA{FT,2}(transpose(ds_lw_in["p_lay"][:][:, lay_ind]))
    t_lev = DA{FT,2}(transpose(ds_lw_in["t_lev"][:][:, lev_ind]))
    t_lay = DA{FT,2}(transpose(ds_lw_in["t_lay"][:][:, lay_ind]))

    t_sfc = DA{FT,1}(ds_lw_in["t_sfc"][:])

    col_dry = DA{FT,2}(undef, nlay, ncol)
    # Reading volume mixing ratios 
    vmrat = zeros(FT, nlay, ncol, ngas)

    vmrat[:, :, idx_gases["h2o"]] .=
        Array{FT}(transpose(ds_lw_in["vmr_h2o"][:]))
    vmrat[:, :, idx_gases["o3"]] .= Array{FT}(transpose(ds_lw_in["vmr_o3"][:]))
    vmrat[:, :, idx_gases["co2"]] .=
        Array{FT}(transpose(ds_lw_in["vmr_co2"][:]))
    vmrat[:, :, idx_gases["n2o"]] .=
        Array{FT}(transpose(ds_lw_in["vmr_n2o"][:]))
    vmrat[:, :, idx_gases["co"]] .= Array{FT}(transpose(ds_lw_in["vmr_co"][:]))
    vmrat[:, :, idx_gases["ch4"]] .=
        Array{FT}(transpose(ds_lw_in["vmr_ch4"][:]))
    vmrat[:, :, idx_gases["o2"]] .= Array{FT}(transpose(ds_lw_in["vmr_o2"][:]))
    vmrat[:, :, idx_gases["n2"]] .= Array{FT}(transpose(ds_lw_in["vmr_n2"][:]))

    @show nlay
    @show ncol
    @show nlev
    @show ngas
    return nothing
end
