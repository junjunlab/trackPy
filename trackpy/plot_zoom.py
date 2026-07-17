"""Zoom panel extension for trackPy — top/bottom layout."""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from trackpy.plot import (
    _resolve_track_color, draw_signal_track, draw_gene_model,
    draw_strand_arrow, draw_ideogram, draw_isoform_row, add_highlights, IGV_COLORS
)


def _draw_zoom_trapezoid_vertical(ax_top, ax_bot, zoom_start, zoom_end,
                                   rs_full, re_full, trap_color_top="#E0E0E0",
                                   trap_color_bot="#404040", trap_smooth=200,
                                   flip=False):
    """Draw a vertical trapezoid between two panels.
    When flip=False: top=narrow(zoom region), bottom=wide(full panel).
    When flip=True:  top=wide(full panel), bottom=narrow(zoom region)."""
    bbox_t = ax_top.get_position()
    bbox_b = ax_bot.get_position()
    top_edge = bbox_t.y0
    bot_edge = bbox_b.y0 + bbox_b.height
    left_x = bbox_t.x0
    right_x = bbox_t.x0 + bbox_t.width

    full_span = re_full - rs_full if re_full > rs_full else 1
    zs_f = (zoom_start - rs_full) / full_span
    ze_f = (zoom_end - rs_full) / full_span
    narrow_l = left_x + (right_x - left_x) * zs_f
    narrow_r = left_x + (right_x - left_x) * ze_f
    wide_l = left_x
    wide_r = right_x

    if flip:
        top_l, top_r = wide_l, wide_r
        bot_l, bot_r = narrow_l, narrow_r
    else:
        top_l, top_r = narrow_l, narrow_r
        bot_l, bot_r = wide_l, wide_r

    r1, g1, b1 = int(trap_color_top[1:3], 16), int(trap_color_top[3:5], 16), int(trap_color_top[5:7], 16)
    r2, g2, b2 = int(trap_color_bot[1:3], 16), int(trap_color_bot[3:5], 16), int(trap_color_bot[5:7], 16)

    for i in range(trap_smooth):
        t = i / trap_smooth
        x_l = top_l + (bot_l - top_l) * t
        x_r = top_r + (bot_r - top_r) * t
        y0 = top_edge + (bot_edge - top_edge) * t
        y1 = top_edge + (bot_edge - top_edge) * (t + 1.0 / trap_smooth)
        r, g, b = int(r1 + (r2 - r1) * t), int(g1 + (g2 - g1) * t), int(b1 + (b2 - b1) * t)
        color = f"#{r:02x}{g:02x}{b:02x}"
        plt.gcf().patches.append(matplotlib.patches.Polygon(
            [[x_l, y1], [x_r, y1], [x_r, y0], [x_l, y0]],
            transform=plt.gcf().transFigure, facecolor=color, edgecolor=color,
            linewidth=0, alpha=1.0, closed=True))


def _draw_signal_block(gs, s0, col, track_labels, track_data, rs, re, ymax,
                       ymax_override, yscale, track_ymax, track_colors, colors,
                       ymax_pos, ymax_label_size, show_yticks, show_box, show_labels):
    """Draw a block of signal tracks."""
    axes = []
    for ti, lbl in enumerate(track_labels):
        ax = plt.gcf().add_subplot(gs[s0 + ti, col]); axes.append(ax)
        ym = ymax
        if yscale == "track" and ymax_override is None:
            ym = track_ymax.get(lbl, ymax)
        cl, cf = _resolve_track_color(str(lbl), track_colors, colors)
        draw_signal_track(ax, track_data.get(lbl, []), rs, re, ym, cl, cf,
                          str(lbl) if show_labels else None,
                          ymax_pos=ymax_pos, ymax_label_size=ymax_label_size,
                          show_yticks=show_yticks, show_box=show_box)
    return axes


def plot_faceted_zoom(genes, data, zoom_region, track_labels, ymax_values, colors, output,
                      title=None, gene_model_bottom=True, figsize=(12, 10),
                      show_coords=True, utr_ratio=0.5,
                      ymax_override=None, ymax_pos=None, ymax_label_size=8,
                      show_yticks=True, show_box=False, track_colors=None, yscale="gene",
                      gene_ratio=0.8, highlights=None, cytoband=None, cytoband_height=0.6,
                      trap_smooth=200, marker_size=0.01,
                      trap_color_top="#E0E0E0", trap_color_bot="#404040", trap_height=2.5,
                      zoom_position="bottom"):
    """Faceted plot with zoom-in panels stacked vertically.
    zoom_region: list of (start, end) tuples, one per gene.
    zoom_position: 'bottom' (full on top, zoom below) or 'top' (zoom on top, full below)."""
    gene_names = list(genes.keys())
    n_genes = len(gene_names)
    n_tracks = len(track_labels)
    # zoom_region is now a list of tuples, one per gene
    head = 1 if show_coords else 0
    ideo = 1 if cytoband else 0

    # Build layout: [header?] + [top-block signals] + [top-block gene] + [top-block cyto?]
    #                + [trap row] + [bot-block signals] + [bot-block gene] + [bot-block cyto?]
    sig_h = [2.0] * n_tracks
    hr = ([0.6] if show_coords else []) + sig_h + [gene_ratio] + ([cytoband_height] if cytoband else []) \
         + [0.8] + sig_h + [gene_ratio] + ([cytoband_height] if cytoband else [])

    # Row indices
    top_s0 = head
    top_gmr = top_s0 + n_tracks
    trap_row = top_gmr + 1 + (1 if cytoband else 0)
    bot_s0 = trap_row + 1
    bot_gmr = bot_s0 + n_tracks

    has_full_on_top = (zoom_position == "bottom")

    fig = plt.figure(figsize=figsize)
    tot = len(hr)
    gs = fig.add_gridspec(tot, n_genes, hspace=0.03, wspace=0.06,
                          height_ratios=hr,
                          left=0.07, right=0.94, top=0.95, bottom=0.05)

    for gi_idx, gn in enumerate(gene_names):
        gi = genes[gn]; gd = data[gn]
        chrom = gi["chr"]
        zoom_start, zoom_end = zoom_region[gi_idx]
        rs_full, re_full = gd["region"]
        ymax_full = ymax_override if ymax_override is not None else gd["ymax"]

        # Zoom region signal data
        zoom_tracks = {}
        for lbl in track_labels:
            vals = gd["tracks"].get(lbl, [])
            zoom_tracks[lbl] = [(s, e, v) for s, e, v in vals if e > zoom_start and s < zoom_end]
        all_zv = [v for lbl in track_labels for _, _, v in zoom_tracks[lbl]]
        zoom_ymax = float(np.percentile(all_zv, 99)) if all_zv else 1.0
        zoom_ymax = ymax_override if ymax_override is not None else zoom_ymax

        if show_coords:
            ax_h = fig.add_subplot(gs[0, gi_idx])
            ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1); ax_h.set_facecolor("white")
            zlabel = f"zoom: {chrom}:{zoom_start:,}-{zoom_end:,}"
            ax_h.text(0.5, 0.72, f"{gn} | {chrom}:{rs_full:,}-{re_full:,}",
                      transform=ax_h.transAxes, fontsize=7, fontweight="bold",
                      fontstyle="italic", color=colors["highlight"], ha="center", va="center")
            ax_h.text(0.5, 0.28, zlabel,
                      transform=ax_h.transAxes, fontsize=6, color="#888888", ha="center", va="center")
            ax_h.set_yticks([]); ax_h.set_xticks([])
            for sp in ax_h.spines.values(): sp.set_visible(False)

        tx = gi["transcripts"][0] if gi.get("transcripts") else {}

        # Track signal axes for trapezoid
        top_signal_axes = []
        bot_signal_axes = []

        # --- Top block ---
        if has_full_on_top:
            top_signal_axes = _draw_signal_block(gs, top_s0, gi_idx, track_labels, gd["tracks"],
                               rs_full, re_full, ymax_full, ymax_override, yscale,
                               gd.get("track_ymax", {}), track_colors, colors,
                               ymax_pos, ymax_label_size, show_yticks, show_box,
                               show_labels=(gi_idx == 0))
            ax_g_top = fig.add_subplot(gs[top_gmr, gi_idx])
            ax_g_top.set_xlim(rs_full, re_full); ax_g_top.set_ylim(0, 3)
            ax_g_top.set_facecolor("none")
            draw_gene_model(ax_g_top, tx, rs_full, re_full, colors, utr_ratio=utr_ratio)
            draw_strand_arrow(ax_g_top, gi["start"], gi["end"], gi.get("strand", "+"), colors)
            ax_g_top.xaxis.set_ticks([]); ax_g_top.set_yticks([])
            for sp in ["top", "right", "left"]: ax_g_top.spines[sp].set_visible(False)
            if cytoband:
                ax_cyto_top = fig.add_subplot(gs[top_gmr + 1, gi_idx])
                draw_ideogram(ax_cyto_top, chrom, gi["start"], gi["end"], cytoband, colors,
                             region_start=rs_full, region_end=re_full,
                             trap_smooth=trap_smooth, marker_size=marker_size,
                             trap_color_top=trap_color_top, trap_color_bot=trap_color_bot,
                             trap_height=trap_height)
        else:
            _draw_signal_block(gs, top_s0, gi_idx, track_labels, zoom_tracks,
                               zoom_start, zoom_end, zoom_ymax, ymax_override, yscale,
                               gd.get("track_ymax", {}), track_colors, colors,
                               ymax_pos, ymax_label_size, show_yticks, show_box,
                               show_labels=(gi_idx == 0))
            ax_g_top = fig.add_subplot(gs[top_gmr, gi_idx])
            ax_g_top.set_xlim(zoom_start, zoom_end); ax_g_top.set_ylim(0, 3)
            ax_g_top.set_facecolor("none")
            draw_gene_model(ax_g_top, tx, zoom_start, zoom_end, colors, utr_ratio=utr_ratio)
            draw_strand_arrow(ax_g_top, gi["start"], gi["end"], gi.get("strand", "+"), colors)
            ax_g_top.xaxis.set_ticks([]); ax_g_top.set_yticks([])
            for sp in ["top", "right", "left"]: ax_g_top.spines[sp].set_visible(False)
            if cytoband:
                ax_cyto_top = fig.add_subplot(gs[top_gmr + 1, gi_idx])
                draw_ideogram(ax_cyto_top, chrom, gi["start"], gi["end"], cytoband, colors,
                             region_start=zoom_start, region_end=zoom_end,
                             trap_smooth=trap_smooth, marker_size=marker_size,
                             trap_color_top=trap_color_top, trap_color_bot=trap_color_bot,
                             trap_height=trap_height)

        # --- Bottom block ---
        if has_full_on_top:
            bot_axes = _draw_signal_block(gs, bot_s0, gi_idx, track_labels, zoom_tracks,
                               zoom_start, zoom_end, zoom_ymax, ymax_override, yscale,
                               gd.get("track_ymax", {}), track_colors, colors,
                               ymax_pos, ymax_label_size, show_yticks, show_box,
                               show_labels=(gi_idx == 0))
            ax_g_bot = fig.add_subplot(gs[bot_gmr, gi_idx])
            ax_g_bot.set_xlim(zoom_start, zoom_end); ax_g_bot.set_ylim(0, 3)
            ax_g_bot.set_facecolor("none")
            draw_gene_model(ax_g_bot, tx, zoom_start, zoom_end, colors, utr_ratio=utr_ratio)
            draw_strand_arrow(ax_g_bot, gi["start"], gi["end"], gi.get("strand", "+"), colors)
            ax_g_bot.xaxis.set_ticks([]); ax_g_bot.set_yticks([])
            for sp in ["top", "right", "left"]: ax_g_bot.spines[sp].set_visible(False)
            if not cytoband:
                ax_g_top.set_xlabel(f"{chrom}:{rs_full:,}-{re_full:,}", fontsize=5, color="#888888")
                ax_g_bot.set_xlabel(f"{chrom}:{zoom_start:,}-{zoom_end:,}", fontsize=5, color="#888888")
            if cytoband:
                ax_cyto_bot = fig.add_subplot(gs[bot_gmr + 1, gi_idx])
                draw_ideogram(ax_cyto_bot, chrom, gi["start"], gi["end"], cytoband, colors,
                             region_start=zoom_start, region_end=zoom_end,
                             trap_smooth=trap_smooth, marker_size=marker_size,
                             trap_color_top=trap_color_top, trap_color_bot=trap_color_bot,
                             trap_height=trap_height)
            # Draw trapezoid: from full gene model bottom to zoom signal top
            _draw_zoom_trapezoid_vertical(ax_g_top, bot_axes[0] if bot_axes else ax_g_bot,
                                          zoom_start, zoom_end,
                                          rs_full, re_full, trap_color_top,
                                          trap_color_bot, trap_smooth)
        else:
            bot_axes2 = _draw_signal_block(gs, bot_s0, gi_idx, track_labels, gd["tracks"],
                               rs_full, re_full, ymax_full, ymax_override, yscale,
                               gd.get("track_ymax", {}), track_colors, colors,
                               ymax_pos, ymax_label_size, show_yticks, show_box,
                               show_labels=(gi_idx == 0))
            ax_g_bot = fig.add_subplot(gs[bot_gmr, gi_idx])
            ax_g_bot.set_xlim(rs_full, re_full); ax_g_bot.set_ylim(0, 3)
            ax_g_bot.set_facecolor("none")
            draw_gene_model(ax_g_bot, tx, rs_full, re_full, colors, utr_ratio=utr_ratio)
            draw_strand_arrow(ax_g_bot, gi["start"], gi["end"], gi.get("strand", "+"), colors)
            ax_g_bot.xaxis.set_ticks([]); ax_g_bot.set_yticks([])
            for sp in ["top", "right", "left"]: ax_g_bot.spines[sp].set_visible(False)
            if not cytoband:
                ax_g_top.set_xlabel(f"{chrom}:{zoom_start:,}-{zoom_end:,}", fontsize=5, color="#888888")
                ax_g_bot.set_xlabel(f"{chrom}:{rs_full:,}-{re_full:,}", fontsize=5, color="#888888")
            if cytoband:
                ax_cyto_bot = fig.add_subplot(gs[bot_gmr + 1, gi_idx])
                draw_ideogram(ax_cyto_bot, chrom, gi["start"], gi["end"], cytoband, colors,
                             region_start=rs_full, region_end=re_full,
                             trap_smooth=trap_smooth, marker_size=marker_size,
                             trap_color_top=trap_color_top, trap_color_bot=trap_color_bot,
                             trap_height=trap_height)
            # Trapezoid: top=wide (zoom panel), bottom=narrow (zoom region in full)
            _draw_zoom_trapezoid_vertical(ax_g_top, bot_axes2[0] if bot_axes2 else ax_g_bot,
                                          zoom_start, zoom_end,
                                          rs_full, re_full, trap_color_top,
                                          trap_color_bot, trap_smooth, flip=True)

    if title: fig.suptitle(str(title), fontsize=10, fontweight="bold", y=0.997, color="#1A1A1A")
    fig.text(0.07, 0.015, "CDS", fontsize=7, fontweight="bold", color=colors["cds"])
    fig.text(0.11, 0.015, "UTR", fontsize=7, fontweight="bold", color=colors["utr"])
    fig.savefig(str(output), dpi=300, facecolor="white", edgecolor="none")
    plt.close(fig)


def plot_isoforms_zoom(genes, data, zoom_region, track_labels, ymax_values, colors, output,
                       title=None, iso_h=0.35, figsize=(12, 12),
                       iso_label_pos="bottom", iso_label_size=3.8,
                       iso_align="top", wspace=0.06,
                       show_coords=True, show_isoform_label=True,
                       ymax_override=None, ymax_pos=None, ymax_label_size=8,
                       show_yticks=True, show_box=False, track_colors=None, yscale="gene",
                       highlights=None, cytoband=None, cytoband_height=0.6,
                       trap_smooth=200, marker_size=0.01,
                       trap_color_top="#E0E0E0", trap_color_bot="#404040", trap_height=2.5,
                       zoom_position="bottom"):
    """Isoforms plot with zoom-in panels stacked vertically.
    zoom_region: list of (start, end) tuples, one per gene."""
    gene_names = list(genes.keys())
    n_genes = len(gene_names)
    n_tracks = len(track_labels)
    max_iso = max(len(genes[gn]["transcripts"]) for gn in gene_names)
    head = 1 if show_coords else 0
    ideo = 1 if cytoband else 0

    # Layout: header rows + [top signals] + [top isoforms] + [trap] + [bot signals] + [bot isoforms] + xaxis + cyto
    hr = ([0.8] if show_coords else []) \
         + [2.0] * n_tracks + [iso_h] * max_iso \
         + [0.8] \
         + [2.0] * n_tracks + [iso_h] * max_iso \
         + [0.3] + ([cytoband_height] if cytoband else [])

    top_s0 = head
    top_iso_s = top_s0 + n_tracks
    trap_row = top_iso_s + max_iso
    bot_s0 = trap_row + 1
    bot_iso_s = bot_s0 + n_tracks
    xaxis_row = bot_iso_s + max_iso

    has_full_on_top = (zoom_position == "bottom")
    fig = plt.figure(figsize=figsize)
    tot = len(hr)
    gs = fig.add_gridspec(tot, n_genes, hspace=0.03, wspace=wspace,
                          height_ratios=hr,
                          left=0.06, right=0.88, top=0.95, bottom=0.05)

    for gi_idx, gn in enumerate(gene_names):
        gi = genes[gn]; gd = data[gn]
        chrom = gi["chr"]
        zoom_start, zoom_end = zoom_region[gi_idx]
        rs_full, re_full = gd["region"]
        ymax_full = ymax_override if ymax_override is not None else gd["ymax"]
        txs = gi["transcripts"]; niso = len(txs)
        top_last_iso_ax = None
        bot_first_sig_ax = None

        # Zoom signal data
        zoom_tracks = {}
        for lbl in track_labels:
            vals = gd["tracks"].get(lbl, [])
            zoom_tracks[lbl] = [(s, e, v) for s, e, v in vals if e > zoom_start and s < zoom_end]
        all_zv = [v for lbl in track_labels for _, _, v in zoom_tracks[lbl]]
        zoom_ymax = float(np.percentile(all_zv, 99)) if all_zv else 1.0
        zoom_ymax = ymax_override if ymax_override is not None else zoom_ymax

        if show_coords:
            ax_h = fig.add_subplot(gs[0, gi_idx])
            ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1); ax_h.set_facecolor("white")
            ax_h.text(0.5, 0.85, f"{chrom}:{rs_full:,}-{re_full:,}",
                      transform=ax_h.transAxes, fontsize=7, fontweight="bold",
                      color=colors["highlight"], ha="center", va="center")
            ax_h.text(0.5, 0.50, f"{gn} | {chrom} | {niso} iso",
                      transform=ax_h.transAxes, fontsize=7, fontweight="bold",
                      fontstyle="italic", color=colors["highlight"], ha="center", va="center")
            zlabel = f"zoom: {zoom_start:,}-{zoom_end:,}"
            ax_h.text(0.5, 0.15, zlabel, transform=ax_h.transAxes, fontsize=6,
                      color="#888888", ha="center", va="center")
            ax_h.set_yticks([]); ax_h.set_xticks([])
            for sp in ax_h.spines.values(): sp.set_visible(False)

        # --- Top block ---
        if has_full_on_top:
            _draw_signal_block(gs, top_s0, gi_idx, track_labels, gd["tracks"],
                               rs_full, re_full, ymax_full, ymax_override, yscale,
                               gd.get("track_ymax", {}), track_colors, colors,
                               ymax_pos, ymax_label_size, show_yticks, show_box,
                               show_labels=(gi_idx == 0))
            _, top_last_iso_ax = _draw_isoform_block(gs, top_iso_s, gi_idx, txs, max_iso,
                                rs_full, re_full, colors, iso_label_pos, iso_label_size,
                                show_isoform_label, iso_align)
        else:
            _draw_signal_block(gs, top_s0, gi_idx, track_labels, zoom_tracks,
                               zoom_start, zoom_end, zoom_ymax, ymax_override, yscale,
                               gd.get("track_ymax", {}), track_colors, colors,
                               ymax_pos, ymax_label_size, show_yticks, show_box,
                               show_labels=(gi_idx == 0))
            _, top_last_iso_ax = _draw_isoform_block(gs, top_iso_s, gi_idx, txs, max_iso,
                                zoom_start, zoom_end, colors, iso_label_pos, iso_label_size,
                                show_isoform_label, iso_align,
                                label_filter_region=(zoom_start, zoom_end))

        # --- Bottom block ---
        if has_full_on_top:
            bot_axes = _draw_signal_block(gs, bot_s0, gi_idx, track_labels, zoom_tracks,
                               zoom_start, zoom_end, zoom_ymax, ymax_override, yscale,
                               gd.get("track_ymax", {}), track_colors, colors,
                               ymax_pos, ymax_label_size, show_yticks, show_box,
                               show_labels=(gi_idx == 0))
            bot_first_sig_ax = bot_axes[0] if bot_axes else None
            _draw_isoform_block(gs, bot_iso_s, gi_idx, txs, max_iso, zoom_start, zoom_end,
                                colors, iso_label_pos, iso_label_size, show_isoform_label, iso_align,
                                label_filter_region=(zoom_start, zoom_end))
        else:
            bot_axes = _draw_signal_block(gs, bot_s0, gi_idx, track_labels, gd["tracks"],
                               rs_full, re_full, ymax_full, ymax_override, yscale,
                               gd.get("track_ymax", {}), track_colors, colors,
                               ymax_pos, ymax_label_size, show_yticks, show_box,
                               show_labels=(gi_idx == 0))
            bot_first_sig_ax = bot_axes[0] if bot_axes else None
            _draw_isoform_block(gs, bot_iso_s, gi_idx, txs, max_iso, rs_full, re_full,
                                colors, iso_label_pos, iso_label_size, show_isoform_label, iso_align)

        # x-axis at bottom
        ax_x = fig.add_subplot(gs[xaxis_row, gi_idx])
        ax_x.set_xlim(0, 1); ax_x.set_ylim(0, 1); ax_x.set_facecolor("none")
        ax_x.text(0.5, 0.5, f"{chrom}:{rs_full:,}-{re_full:,}",
                  transform=ax_x.transAxes, fontsize=6, color="#555555", ha="center", va="center")
        ax_x.set_yticks([]); ax_x.set_xticks([])
        for sp in ax_x.spines.values(): sp.set_visible(False)

        if cytoband:
            ax_cyto = fig.add_subplot(gs[xaxis_row + 1, gi_idx])
            draw_ideogram(ax_cyto, chrom, gi["start"], gi["end"], cytoband, colors,
                         region_start=rs_full, region_end=re_full,
                         trap_smooth=trap_smooth, marker_size=marker_size,
                         trap_color_top=trap_color_top, trap_color_bot=trap_color_bot,
                         trap_height=trap_height)

        # Trapezoid: from last top isoform row to first bottom signal row
        if top_last_iso_ax and bot_first_sig_ax:
            _draw_zoom_trapezoid_vertical(top_last_iso_ax, bot_first_sig_ax,
                                          zoom_start, zoom_end, rs_full, re_full,
                                          trap_color_top, trap_color_bot, trap_smooth,
                                          flip=not has_full_on_top)

    if title: fig.suptitle(str(title), fontsize=10, fontweight="bold", y=0.995, color="#1A1A1A")
    fig.text(0.06, 0.015, "CDS", fontsize=7, fontweight="bold", color=colors["cds"])
    fig.text(0.09, 0.015, "UTR", fontsize=7, fontweight="bold", color=colors["utr"])
    fig.savefig(str(output), dpi=300, facecolor="white", edgecolor="none")
    plt.close(fig)


def _draw_isoform_block(gs, iso_s, col, txs, max_iso, rs, re, colors,
                         iso_label_pos, iso_label_size, show_label, iso_align,
                         label_filter_region=None):
    """Draw a block of isoform rows. Returns (first_ax, last_ax).
    label_filter_region: (start, end) — only show labels for transcripts overlapping this region."""
    niso = len(txs)
    if iso_align == "top":
        off_top = 0
    elif iso_align == "center":
        off_top = (max_iso - niso) // 2
    else:
        off_top = max_iso - niso
    first_ax = None; last_ax = None
    for ti, tx in enumerate(txs):
        ax = plt.gcf().add_subplot(gs[iso_s + off_top + ti, col])
        if first_ax is None: first_ax = ax
        last_ax = ax
        ax.set_xlim(rs, re); ax.set_ylim(0, 2); ax.set_facecolor("none")
        sl = show_label
        if label_filter_region and sl:
            zs, ze = label_filter_region
            exons = tx.get("exons", [])
            if exons:
                tx_start = min(e[0] for e in exons)
                tx_end = max(e[1] for e in exons)
                if tx_end <= zs or tx_start >= ze:
                    sl = False
        draw_isoform_row(ax, tx, rs, re, colors, iso_label_pos, iso_label_size,
                         show_label=sl)
        ax.set_yticks([]); ax.set_xticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
    for er in range(max_iso - niso):
        if er < off_top or er >= off_top + niso:
            ax = plt.gcf().add_subplot(gs[iso_s + er, col])
            ax.set_facecolor("none"); ax.set_yticks([]); ax.set_xticks([])
            for sp in ax.spines.values(): sp.set_visible(False)
    return first_ax, last_ax
