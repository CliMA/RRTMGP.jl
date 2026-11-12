using RRTMGP
using NCDatasets
using CairoMakie

const LT = RRTMGP.LookUpTables
const OP = RRTMGP.Optics
const AP = RRTMGP.ArtifactPaths

lkp_lw = NCDataset(AP.get_lookup_filename(:cloud, :lw), "r") do ds
    LT.LookUpCld(ds, Float64, Array)
end
lkp_sw = NCDataset(AP.get_lookup_filename(:cloud, :sw), "r") do ds
    LT.LookUpCld(ds, Float64, Array)
end
lkp_precip_lw = NCDataset(AP.get_lookup_filename(:precip, :lw), "r") do ds
    LT.LookUpPrecip(ds, Float64, Array)
end
lkp_precip_sw = NCDataset(AP.get_lookup_filename(:precip, :sw), "r") do ds
    LT.LookUpPrecip(ds, Float64, Array)
end

function eval_lut_curve(
    kernel,
    nsize,
    r_lwr,
    r_upr,
    ext_v,
    ssa_v,
    asy_v,
    r_grid,
)
    ext_out = similar(r_grid)
    ssa_out = similar(r_grid)
    asy_out = similar(r_grid)
    for (i, r) in enumerate(r_grid)
        τ, τ_ssa, τ_ssag =
            kernel(nsize, r_lwr, r_upr, ext_v, ssa_v, asy_v, r, 1.0)
        ext_out[i] = τ
        ssa_out[i] = τ > 0 ? τ_ssa / τ : 0.0
        asy_out[i] = τ_ssa > 0 ? τ_ssag / τ_ssa : 0.0
    end
    return ext_out, ssa_out, asy_out
end

# Representative bands: SW Band 10 (625–778 nm) and LW Band 6 (820–980 cm⁻¹, ~10.2–12.2 µm)
ib_sw, ib_lw = 10, 6

fig1 = Figure(size = (1050, 820), fontsize = 13)
Label(
    fig1[0, 1:2],
    "RRTMGP Cloud & Precipitation Lookup Tables vs. Effective Radius",
    fontsize = 16,
    font = :bold,
)

for (col, (lkp_cld, lkp_prc, ib, title_str)) in enumerate((
    (lkp_sw, lkp_precip_sw, ib_sw, "Shortwave Band $ib_sw (625–778 nm)"),
    (
        lkp_lw,
        lkp_precip_lw,
        ib_lw,
        "Longwave Band $ib_lw (820–980 cm⁻¹, 10.2–12.2 µm)",
    ),
))
    _, _, nsize_liq, nsize_ice, _ = lkp_cld.dims
    radliq_lwr, radliq_upr, radice_lwr, radice_upr = lkp_cld.bounds
    _, nsize_sno, nsize_rai = lkp_prc.dims
    radsno_lwr, radsno_upr, radrai_lwr, radrai_upr = lkp_prc.bounds

    r_liq = collect(range(radliq_lwr, radliq_upr, length = 100))
    r_ice = collect(range(radice_lwr, radice_upr, length = 100))
    r_sno = exp.(range(log(radsno_lwr), log(radsno_upr), length = 200))
    r_rai = exp.(range(log(radrai_lwr), log(radrai_upr), length = 200))

    ext_l, ssa_l, asy_l = eval_lut_curve(
        OP.compute_lookup_cld_liq_props,
        nsize_liq,
        radliq_lwr,
        radliq_upr,
        LT.getview_liqdata(lkp_cld, ib)...,
        r_liq,
    )
    ext_s, ssa_s, asy_s = eval_lut_curve(
        OP.compute_lookup_precip_props,
        nsize_sno,
        radsno_lwr,
        radsno_upr,
        LT.getview_snodata(lkp_prc, ib)...,
        r_sno,
    )
    ext_r, ssa_r, asy_r = eval_lut_curve(
        OP.compute_lookup_precip_props,
        nsize_rai,
        radrai_lwr,
        radrai_upr,
        LT.getview_raidata(lkp_prc, ib)...,
        r_rai,
    )

    ticks = (
        [2, 5, 10, 20, 50, 100, 200, 500, 1000, 5000],
        ["2", "5", "10", "20", "50", "100", "200", "500", "1000", "5000"],
    )
    ax_ext = Axis(
        fig1[1, col],
        title = title_str,
        ylabel = col == 1 ? "Extinction coeff. [m²/g]" : "",
        xscale = log10,
        yscale = log10,
        xticks = ticks,
    )
    ax_ssa = Axis(
        fig1[2, col],
        ylabel = col == 1 ? "Single-scattering albedo ω₀ [-]" : "",
        xscale = log10,
        xticks = ticks,
    )
    ax_asy = Axis(
        fig1[3, col],
        xlabel = "Effective radius r_eff [µm]",
        ylabel = col == 1 ? "Asymmetry parameter g [-]" : "",
        xscale = log10,
        xticks = ticks,
    )

    # Rain & Cloud liquid
    lines!(
        ax_ext,
        r_rai,
        ext_r,
        color = :dodgerblue3,
        linestyle = :dash,
        linewidth = 2.2,
        label = "Rain LUT (Mie, 5–5000 µm)",
    )
    lines!(
        ax_ssa,
        r_rai,
        ssa_r,
        color = :dodgerblue3,
        linestyle = :dash,
        linewidth = 2.2,
    )
    lines!(
        ax_asy,
        r_rai,
        asy_r,
        color = :dodgerblue3,
        linestyle = :dash,
        linewidth = 2.2,
    )

    lines!(
        ax_ext,
        r_liq,
        ext_l,
        color = :navy,
        linewidth = 2.8,
        label = "Cloud liquid LUT (2.5–21.5 µm)",
    )
    lines!(ax_ssa, r_liq, ssa_l, color = :navy, linewidth = 2.8)
    lines!(ax_asy, r_liq, asy_l, color = :navy, linewidth = 2.8)

    # Snow & Cloud ice (3 roughness levels)
    lines!(
        ax_ext,
        r_sno,
        ext_s,
        color = :crimson,
        linestyle = :dash,
        linewidth = 2.2,
        label = "Snow LUT (Yang 2013, 5–1100 µm)",
    )
    lines!(
        ax_ssa,
        r_sno,
        ssa_s,
        color = :crimson,
        linestyle = :dash,
        linewidth = 2.2,
    )
    lines!(
        ax_asy,
        r_sno,
        asy_s,
        color = :crimson,
        linestyle = :dash,
        linewidth = 2.2,
    )

    rgh_labels = ("smooth", "medium", "rough")
    rgh_colors = (:orange, :darkorange3, :firebrick4)
    for irgh in 1:3
        ext_i, ssa_i, asy_i = eval_lut_curve(
            OP.compute_lookup_cld_ice_props,
            nsize_ice,
            radice_lwr,
            radice_upr,
            LT.getview_icedata(lkp_cld, ib, irgh)...,
            r_ice,
        )
        lines!(
            ax_ext,
            r_ice,
            ext_i,
            color = rgh_colors[irgh],
            linewidth = 2.0,
            label = "Cloud ice LUT ($(rgh_labels[irgh]), 5–90 µm)",
        )
        lines!(ax_ssa, r_ice, ssa_i, color = rgh_colors[irgh], linewidth = 2.0)
        lines!(ax_asy, r_ice, asy_i, color = rgh_colors[irgh], linewidth = 2.0)
    end
end

Legend(
    fig1[4, 1:2],
    fig1.content[2],
    orientation = :horizontal,
    nbanks = 2,
    framevisible = true,
)

fig2 = Figure(size = (1050, 780), fontsize = 13)
Label(
    fig2[0, 1:2],
    "Existing RRTMGP Cloud Liquid & Ice Lookup Tables Across All Spectral Bands",
    fontsize = 16,
    font = :bold,
)

for (col, (lkp, wave_label)) in
    enumerate(((lkp_lw, "Longwave"), (lkp_sw, "Shortwave")))
    nband, _, nsize_liq, nsize_ice = lkp.dims
    radliq_lwr, radliq_upr, radice_lwr, radice_upr = lkp.bounds
    bands = 1:nband

    ax_e = Axis(
        fig2[1, col],
        title = "$wave_label Bands (1–$nband)",
        ylabel = col == 1 ? "Extinction [m²/g]" : "",
        yscale = log10,
    )
    ax_s = Axis(
        fig2[2, col],
        ylabel = col == 1 ? "Single-scattering albedo ω₀ [-]" : "",
    )
    ax_g = Axis(
        fig2[3, col],
        xlabel = "$wave_label band index",
        ylabel = col == 1 ? "Asymmetry parameter g [-]" : "",
    )

    for (r_l, c_l) in ((5.0, :royalblue1), (10.0, :dodgerblue3), (20.0, :navy))
        ext_b = [
            OP.compute_lookup_cld_liq_props(
                nsize_liq,
                radliq_lwr,
                radliq_upr,
                LT.getview_liqdata(lkp, b)...,
                r_l,
                1.0,
            )[1] for b in bands
        ]
        ssa_b = [
            begin
                τ, τs, _ = OP.compute_lookup_cld_liq_props(
                    nsize_liq,
                    radliq_lwr,
                    radliq_upr,
                    LT.getview_liqdata(lkp, b)...,
                    r_l,
                    1.0,
                )
                τs / τ
            end for b in bands
        ]
        asy_b = [
            begin
                _, τs, τsg = OP.compute_lookup_cld_liq_props(
                    nsize_liq,
                    radliq_lwr,
                    radliq_upr,
                    LT.getview_liqdata(lkp, b)...,
                    r_l,
                    1.0,
                )
                τsg / τs
            end for b in bands
        ]
        scatterlines!(
            ax_e,
            bands,
            ext_b,
            color = c_l,
            linewidth = 2,
            markersize = 6,
            label = "Liquid r_eff = $(Int(r_l)) µm",
        )
        scatterlines!(
            ax_s,
            bands,
            ssa_b,
            color = c_l,
            linewidth = 2,
            markersize = 6,
        )
        scatterlines!(
            ax_g,
            bands,
            asy_b,
            color = c_l,
            linewidth = 2,
            markersize = 6,
        )
    end

    for (r_i, c_i) in
        ((10.0, :gold3), (30.0, :darkorange2), (80.0, :firebrick3))
        ext_b = [
            OP.compute_lookup_cld_ice_props(
                nsize_ice,
                radice_lwr,
                radice_upr,
                LT.getview_icedata(lkp, b, 2)...,
                r_i,
                1.0,
            )[1] for b in bands
        ]
        ssa_b = [
            begin
                τ, τs, _ = OP.compute_lookup_cld_ice_props(
                    nsize_ice,
                    radice_lwr,
                    radice_upr,
                    LT.getview_icedata(lkp, b, 2)...,
                    r_i,
                    1.0,
                )
                τs / τ
            end for b in bands
        ]
        asy_b = [
            begin
                _, τs, τsg = OP.compute_lookup_cld_ice_props(
                    nsize_ice,
                    radice_lwr,
                    radice_upr,
                    LT.getview_icedata(lkp, b, 2)...,
                    r_i,
                    1.0,
                )
                τsg / τs
            end for b in bands
        ]
        scatterlines!(
            ax_e,
            bands,
            ext_b,
            color = c_i,
            linestyle = :dash,
            linewidth = 2,
            markersize = 6,
            label = "Ice (med. rgh) r_eff = $(Int(r_i)) µm",
        )
        scatterlines!(
            ax_s,
            bands,
            ssa_b,
            color = c_i,
            linestyle = :dash,
            linewidth = 2,
            markersize = 6,
        )
        scatterlines!(
            ax_g,
            bands,
            asy_b,
            color = c_i,
            linestyle = :dash,
            linewidth = 2,
            markersize = 6,
        )
    end
end

Legend(
    fig2[4, 1:2],
    fig2.content[2],
    orientation = :horizontal,
    nbanks = 2,
    framevisible = true,
)

fig3 = Figure(size = (1100, 820), fontsize = 12)
Label(
    fig3[0, 1:4],
    "RRTMGP Snow & Rain Optical Lookup Tables Across All Bands and Effective Radii",
    fontsize = 16,
    font = :bold,
)

datasets = (
    ("Snow LW", lkp_precip_lw, :sno),
    ("Snow SW", lkp_precip_sw, :sno),
    ("Rain LW", lkp_precip_lw, :rai),
    ("Rain SW", lkp_precip_sw, :rai),
)

for (col, (label_str, lkp_prc, kind)) in enumerate(datasets)
    nband, nsize_sno, nsize_rai = lkp_prc.dims
    radsno_lwr, radsno_upr, radrai_lwr, radrai_upr = lkp_prc.bounds
    nsize = kind === :sno ? nsize_sno : nsize_rai
    r_lwr = kind === :sno ? radsno_lwr : radrai_lwr
    r_upr = kind === :sno ? radsno_upr : radrai_upr
    r_grid = exp.(range(log(r_lwr), log(r_upr), length = nsize))

    ext_mat = zeros(nsize, nband)
    ssa_mat = zeros(nsize, nband)
    asy_mat = zeros(nsize, nband)
    for b in 1:nband
        ev, sv, gv =
            kind === :sno ? LT.getview_snodata(lkp_prc, b) :
            LT.getview_raidata(lkp_prc, b)
        ext_mat[:, b] .= ev
        ssa_mat[:, b] .= sv
        asy_mat[:, b] .= gv
    end

    ax1 = Axis(
        fig3[1, col],
        title = "$label_str: log₁₀(ext [m²/g])",
        ylabel = col == 1 ? "r_eff [µm]" : "",
        yscale = log10,
    )
    ax2 = Axis(
        fig3[2, col],
        title = "$label_str: ω₀ [-]",
        ylabel = col == 1 ? "r_eff [µm]" : "",
        yscale = log10,
    )
    ax3 = Axis(
        fig3[3, col],
        title = "$label_str: g [-]",
        xlabel = "Band index",
        ylabel = col == 1 ? "r_eff [µm]" : "",
        yscale = log10,
    )

    hm1 = heatmap!(
        ax1,
        1:nband,
        r_grid,
        log10.(max.(ext_mat', 1e-5)),
        colormap = :viridis,
        colorrange = (-3.5, 0.5),
    )
    hm2 = heatmap!(
        ax2,
        1:nband,
        r_grid,
        ssa_mat',
        colormap = :plasma,
        colorrange = (0.4, 1.0),
    )
    hm3 = heatmap!(
        ax3,
        1:nband,
        r_grid,
        asy_mat',
        colormap = :cividis,
        colorrange = (0.5, 1.0),
    )

    if col == 4
        Colorbar(fig3[1, 5], hm1, label = "log₁₀(ext [m²/g])")
        Colorbar(fig3[2, 5], hm2, label = "ω₀ [-]")
        Colorbar(fig3[3, 5], hm3, label = "g [-]")
    end
end
