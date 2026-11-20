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

# Representative bands: SW Band 10 (12850–16000 cm⁻¹, 625–778 nm) and
# LW Band 6 (820–980 cm⁻¹, 10.2–12.2 µm)
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

# Shared 2D (r_eff, wavenumber) lookup-table grid plotter reproducing the axes,
# merged LW+SW spectral grid, and `:jet` color scales of Figures 4 and 6 in
# Feng, Menzel, & Paynter (2025):
# - x-axis: effective radius r_eff [µm] on a shared [1, 5000] µm log₁₀ axis
#           (out-of-range sizes left blank/transparent)
# - y-axis: spectral wavenumber [cm⁻¹] on a single [200, 50000] cm⁻¹ log₁₀ axis
#           combining all 16 LW bands (< 3250 cm⁻¹; Band 1 clipped at 200 cm⁻¹
#           for display) and the 12 SW bands above 3250 cm⁻¹ (3250–50000 cm⁻¹)
# - Row 1: mass extinction coefficient [m²/g] (log₁₀ scale, `:jet`)
# - Row 2: single-scattering albedo ω₀ [-] on co-albedo log₁₀(1 - ω₀) scale
#          (`:jet`, ticks 0.9999, 0.999, 0.99, 0.9, 0 from bottom to top)
# - Row 3: asymmetry factor g [-] (linear 0–1 scale, `:jet`)
function build_lut_2d_figure(fig_title, specs; nr = 300)
    ncol = length(specs)
    fig = Figure(size = (840, 780), fontsize = 12)
    Label(fig[0, 1:ncol], fig_title, fontsize = 16, font = :bold)

    xlims = (1.0, 5000.0)
    ylims = (200.0, 50000.0)
    xticks = ([1, 10, 100, 1000], ["10⁰", "10¹", "10²", "10³"])
    yticks = (
        [400, 600, 1000, 10000, 20000],
        ["400", "600", "1000", "10000", "20000"],
    )

    base_r = exp.(range(log(xlims[1]), log(xlims[2]), length = nr + 1))
    hm_ext, hm_ssa, hm_asy = nothing, nothing, nothing

    for (col, spec) in enumerate(specs)
        (;
            col_title,
            lkp_lw_band,
            lkp_sw_band,
            r_lwr,
            r_upr,
            eval_lw,
            eval_sw,
        ) = spec

        nb_lw = size(lkp_lw_band, 2)
        sw_idx = [
            b for b in 1:size(lkp_sw_band, 2) if
            lkp_sw_band[1, b] >= 3250.0 - 1e-6
        ]
        sort!(sw_idx; by = b -> lkp_sw_band[1, b])

        nband = nb_lw + length(sw_idx)
        wn_edges = Vector{Float64}(undef, nband + 1)
        # Clamp far-IR Band 1 lower display edge (10 cm⁻¹) to 200 cm⁻¹ so
        # Band 1 does not dominate the log₁₀ y-axis, matching Figures 4 & 6 of
        # Feng et al. (2025).
        wn_edges[1] = max(Float64(lkp_lw_band[1, 1]), ylims[1])
        for b in 1:nb_lw
            wn_edges[b + 1] = Float64(lkp_lw_band[2, b])
        end
        for (k, b) in enumerate(sw_idx)
            wn_edges[nb_lw + k + 1] = Float64(lkp_sw_band[2, b])
        end

        # Insert exact [r_lwr, r_upr] boundaries so the valid domain edge is
        # crisp, with NaN (transparent/white) outside [r_lwr, r_upr].
        r_edges = filter(
            r ->
                abs(r - r_lwr) > 1e-7 * r_lwr &&
                    abs(r - r_upr) > 1e-7 * r_upr,
            base_r,
        )
        append!(r_edges, (Float64(r_lwr), Float64(r_upr)))
        sort!(r_edges)
        n_cells = length(r_edges) - 1
        r_centers = sqrt.(r_edges[1:(end - 1)] .* r_edges[2:end])
        in_range = (r_centers .>= r_lwr) .& (r_centers .<= r_upr)
        r_valid = r_centers[in_range]

        ext_mat = fill(NaN, n_cells, nband)
        ssa_mat = fill(NaN, n_cells, nband)
        asy_mat = fill(NaN, n_cells, nband)

        for b in 1:nb_lw
            e_b, s_b, g_b = eval_lw(b, r_valid)
            ext_mat[in_range, b] .= log10.(max.(e_b, 1e-4))
            ssa_mat[in_range, b] .= log10.(clamp.(1.0 .- s_b, 1e-4, 1.0))
            asy_mat[in_range, b] .= g_b
        end
        for (k, b) in enumerate(sw_idx)
            col_b = nb_lw + k
            e_b, s_b, g_b = eval_sw(b, r_valid)
            ext_mat[in_range, col_b] .= log10.(max.(e_b, 1e-4))
            ssa_mat[in_range, col_b] .= log10.(clamp.(1.0 .- s_b, 1e-4, 1.0))
            asy_mat[in_range, col_b] .= g_b
        end

        ax1 = Axis(
            fig[1, col],
            title = col_title,
            ylabel = col == 1 ? "Wavenumber [cm⁻¹]" : "",
            xscale = log10,
            yscale = log10,
            xticks = xticks,
            yticks = yticks,
            limits = (xlims, ylims),
        )
        ax2 = Axis(
            fig[2, col],
            ylabel = col == 1 ? "Wavenumber [cm⁻¹]" : "",
            xscale = log10,
            yscale = log10,
            xticks = xticks,
            yticks = yticks,
            limits = (xlims, ylims),
        )
        ax3 = Axis(
            fig[3, col],
            xlabel = "r_eff [µm]",
            ylabel = col == 1 ? "Wavenumber [cm⁻¹]" : "",
            xscale = log10,
            yscale = log10,
            xticks = xticks,
            yticks = yticks,
            limits = (xlims, ylims),
        )
        hidexdecorations!(ax1, grid = false, ticks = false)
        hidexdecorations!(ax2, grid = false, ticks = false)
        if col > 1
            hideydecorations!(ax1, grid = false, ticks = false)
            hideydecorations!(ax2, grid = false, ticks = false)
            hideydecorations!(ax3, grid = false, ticks = false)
        end

        hm_ext = heatmap!(
            ax1,
            r_edges,
            wn_edges,
            ext_mat,
            colormap = :jet,
            colorrange = (-3.2, 0.2),
            nan_color = :transparent,
        )
        hm_ssa = heatmap!(
            ax2,
            r_edges,
            wn_edges,
            ssa_mat,
            colormap = :jet,
            colorrange = (-4.0, 0.0),
            nan_color = :transparent,
        )
        hm_asy = heatmap!(
            ax3,
            r_edges,
            wn_edges,
            asy_mat,
            colormap = :jet,
            colorrange = (0.0, 1.0),
            nan_color = :transparent,
        )
    end

    Colorbar(
        fig[1, ncol + 1],
        hm_ext,
        label = "Extinction [m²/g]",
        ticks = ([-3, -2, -1, 0], ["10⁻³", "10⁻²", "10⁻¹", "10⁰"]),
    )
    Colorbar(
        fig[2, ncol + 1],
        hm_ssa,
        label = "Single-scattering albedo ω₀ [-]",
        ticks = ([-4, -3, -2, -1, 0], ["0.9999", "0.999", "0.99", "0.9", "0"]),
    )
    Colorbar(
        fig[3, ncol + 1],
        hm_asy,
        label = "Asymmetry factor g [-]",
        ticks = 0.0:0.2:1.0,
    )
    colgap!(fig.layout, 16)
    rowgap!(fig.layout, 12)
    return fig
end

function cld_liq_spec(col_title, lkp_lw, lkp_sw)
    _, _, nsize_liq_lw, _, _ = lkp_lw.dims
    _, _, nsize_liq_sw, _, _ = lkp_sw.dims
    r_lwr, r_upr = lkp_lw.bounds[1], lkp_lw.bounds[2]
    return (;
        col_title,
        lkp_lw_band = lkp_lw.bnd_lims_wn,
        lkp_sw_band = lkp_sw.bnd_lims_wn,
        r_lwr,
        r_upr,
        eval_lw = (b, r_vec) -> eval_lut_curve(
            OP.compute_lookup_cld_liq_props,
            nsize_liq_lw,
            r_lwr,
            r_upr,
            LT.getview_liqdata(lkp_lw, b)...,
            r_vec,
        ),
        eval_sw = (b, r_vec) -> eval_lut_curve(
            OP.compute_lookup_cld_liq_props,
            nsize_liq_sw,
            r_lwr,
            r_upr,
            LT.getview_liqdata(lkp_sw, b)...,
            r_vec,
        ),
    )
end

function cld_ice_spec(col_title, lkp_lw, lkp_sw; irgh = 2)
    _, _, _, nsize_ice_lw, _ = lkp_lw.dims
    _, _, _, nsize_ice_sw, _ = lkp_sw.dims
    r_lwr, r_upr = lkp_lw.bounds[3], lkp_lw.bounds[4]
    return (;
        col_title,
        lkp_lw_band = lkp_lw.bnd_lims_wn,
        lkp_sw_band = lkp_sw.bnd_lims_wn,
        r_lwr,
        r_upr,
        eval_lw = (b, r_vec) -> eval_lut_curve(
            OP.compute_lookup_cld_ice_props,
            nsize_ice_lw,
            r_lwr,
            r_upr,
            LT.getview_icedata(lkp_lw, b, irgh)...,
            r_vec,
        ),
        eval_sw = (b, r_vec) -> eval_lut_curve(
            OP.compute_lookup_cld_ice_props,
            nsize_ice_sw,
            r_lwr,
            r_upr,
            LT.getview_icedata(lkp_sw, b, irgh)...,
            r_vec,
        ),
    )
end

function precip_spec(col_title, lkp_lw, lkp_sw, kind)
    _, nsize_sno_lw, nsize_rai_lw = lkp_lw.dims
    _, nsize_sno_sw, nsize_rai_sw = lkp_sw.dims
    radsno_lwr, radsno_upr, radrai_lwr, radrai_upr = lkp_lw.bounds
    nsize_lw = kind === :sno ? nsize_sno_lw : nsize_rai_lw
    nsize_sw = kind === :sno ? nsize_sno_sw : nsize_rai_sw
    r_lwr = kind === :sno ? radsno_lwr : radrai_lwr
    r_upr = kind === :sno ? radsno_upr : radrai_upr
    get_view = kind === :sno ? LT.getview_snodata : LT.getview_raidata
    return (;
        col_title,
        lkp_lw_band = lkp_lw.bnd_lims_wn,
        lkp_sw_band = lkp_sw.bnd_lims_wn,
        r_lwr,
        r_upr,
        eval_lw = (b, r_vec) -> eval_lut_curve(
            OP.compute_lookup_precip_props,
            nsize_lw,
            r_lwr,
            r_upr,
            get_view(lkp_lw, b)...,
            r_vec,
        ),
        eval_sw = (b, r_vec) -> eval_lut_curve(
            OP.compute_lookup_precip_props,
            nsize_sw,
            r_lwr,
            r_upr,
            get_view(lkp_sw, b)...,
            r_vec,
        ),
    )
end

fig2 = build_lut_2d_figure(
    "Existing RRTMGP Cloud Liquid & Ice (Medium Roughness) Lookup Tables",
    (
        cld_liq_spec("Cloud Liquid (2.5–21.5 µm)", lkp_lw, lkp_sw),
        cld_ice_spec("Cloud Ice, med. rgh (5–90 µm)", lkp_lw, lkp_sw; irgh = 2),
    ),
)

fig3 = build_lut_2d_figure(
    "RRTMGP Rain & Snow Optical Lookup Tables (Feng et al., 2025)",
    (
        precip_spec("Rain (5–5000 µm)", lkp_precip_lw, lkp_precip_sw, :rai),
        precip_spec("Snow (5–1100 µm)", lkp_precip_lw, lkp_precip_sw, :sno),
    ),
)
