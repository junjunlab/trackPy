"""Publication-quality IGV-style track plot engine."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

# ── IGV color scheme ────────────────────────────────────────────
IGV_COLORS = {
    "input_line": "#808080", "input_fill": "#D0D0D0",
    "ip_line": "#1A1A1A", "ip_fill": "#606060",
    "cds": "#1A1A1A", "utr": "#1A1A1A",
    "intron": "#1A1A1A", "noncoding": "#C0C0C0",
    "highlight": "#404040", "bracket": "#333333",
}

def _resolve_track_color(lbl, track_colors, colors):
    """Return (line_color, fill_color) for a track label."""
    if track_colors and lbl in track_colors:
        c = track_colors[lbl]
        return c, c
    is_ip = "m6A" in str(lbl) or "IP" in str(lbl) or "ip" in str(lbl)
    cl = colors.get("ip_line") if is_ip else colors.get("input_line")
    cf = colors.get("ip_fill") if is_ip else colors.get("input_fill")
    return cl, cf


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "axes.facecolor": "white", "figure.facecolor": "white",
})


# ── Drawing helpers ─────────────────────────────────────────────
def _has_overlap(s, e, regions):
    for rs, re in regions:
        if max(s, rs) < min(e, re): return True
    return False


def draw_gene_model(ax, gene_info, region_start, region_end, colors=None,
                    utr_ratio=0.8, gene_names=None):
    c = colors or IGV_COLORS
    exons = gene_info.get("exons", [])
    cds   = gene_info.get("CDS", [])
    utr5  = gene_info.get("UTR_5", [])
    utr3  = gene_info.get("UTR_3", [])
    if not exons: return

    cds_h = 0.8
    utr_h = cds_h * utr_ratio
    cds_y = 1.5 - cds_h / 2
    utr_y = 1.5 - utr_h / 2

    ax.plot([min(e[0] for e in exons), max(e[1] for e in exons)],
            [1.5, 1.5], linewidth=0.5, color=c["intron"])
    for es, ee in exons:
        in_cds = _has_overlap(es, ee, cds)
        if in_cds:
            for cs, ce in cds:
                o1, o2 = max(es, cs), min(ee, ce)
                if o1 < o2: ax.add_patch(Rectangle((o1, cds_y), o2 - o1, cds_h,
                                          facecolor=c["cds"], edgecolor="none"))
        for reg in (utr5, utr3):
            for us, ue in reg:
                o1, o2 = max(es, us), min(ee, ue)
                if o1 < o2: ax.add_patch(Rectangle((o1, utr_y), o2 - o1, utr_h,
                                          facecolor=c["utr"], edgecolor="none"))
        if not in_cds and not _has_overlap(es, ee, utr5) and not _has_overlap(es, ee, utr3):
            ax.add_patch(Rectangle((es, utr_y), ee - es, utr_h,
                         facecolor=c["utr"], edgecolor="none"))
    if gene_names:
        mid = (region_start + region_end) / 2
        ax.text(mid, 0.85, gene_names, fontsize=5, color="#555555", ha="center", va="top",
                style="italic", clip_on=True)


def draw_strand_arrow(ax, gs, ge, strand, colors=None, arrow_y=2.5):
    c = colors or IGV_COLORS
    if strand == "+":
        ax.annotate("", xy=(ge, arrow_y), xytext=(gs, arrow_y),
                    arrowprops=dict(arrowstyle="->", color=c["input_line"], lw=1.0))
    else:
        ax.annotate("", xy=(gs, arrow_y), xytext=(ge, arrow_y),
                    arrowprops=dict(arrowstyle="->", color=c["input_line"], lw=1.0))


def draw_isoform_row(ax, tx, rs, re, colors=None, label_pos="bottom", label_size=6,
                     show_label=True, strand=None):
    c = colors or IGV_COLORS
    exons = tx.get("exons", [])
    cds   = tx.get("CDS", [])
    utr5  = tx.get("UTR_5", [])
    utr3  = tx.get("UTR_3", [])
    is_coding = tx.get("is_coding", True)
    strand = strand or tx.get("strand")
    if exons:
        g_left = min(e[0] for e in exons)
        g_right = max(e[1] for e in exons)
        ax.plot([g_left, g_right], [1.0, 1.0], linewidth=0.3, color=c["intron"])
        # Draw strand arrows on intron line
        if strand:
            span = g_right - g_left if g_right > g_left else 1
            arrow_spacing = max(span / 8, 1)
            x = g_left + arrow_spacing * 0.5
            direction = 1 if strand == "+" else -1
            arrow_dx = span * 0.015 * direction
            while x < g_right - arrow_spacing * 0.5:
                ax.annotate("", xy=(x + arrow_dx, 1.0), xytext=(x - arrow_dx, 1.0),
                            arrowprops=dict(arrowstyle="->", color=c["intron"],
                            lw=0.4), clip_on=True)
                x += arrow_spacing
    for es, ee in exons:
        if is_coding:
            for cs, ce in cds:
                o1, o2 = max(es, cs), min(ee, ce)
                if o1 < o2: ax.add_patch(Rectangle((o1, 0.6), o2 - o1, 0.8,
                                          facecolor=c["cds"], edgecolor="none"))
            for reg in (utr5, utr3):
                for us, ue in reg:
                    o1, o2 = max(es, us), min(ee, ue)
                    if o1 < o2: ax.add_patch(Rectangle((o1, 0.75), o2 - o1, 0.5,
                                              facecolor=c["utr"], edgecolor="none"))
        else:
            ax.add_patch(Rectangle((es, 0.75), ee - es, 0.5,
                         facecolor=c["utr"], edgecolor="none"))
    if tx.get("id") and show_label:
        ex_starts = [e[0] for e in exons] if exons else [rs]
        ex_ends   = [e[1] for e in exons] if exons else [re]
        g_left  = min(ex_starts)
        g_right = max(ex_ends)
        span = g_right - g_left if g_right > g_left else 1

        if label_pos == "left":
            ax.text(g_left - span * 0.02, 1.0, tx["id"],
                    fontsize=label_size, color="#555555", va="center", ha="right",
                    family="monospace", clip_on=False)
        elif label_pos == "right":
            ax.text(g_right + span * 0.02, 1.0, tx["id"],
                    fontsize=label_size, color="#555555", va="center", ha="left",
                    family="monospace", clip_on=False)
        elif label_pos == "top":
            ax.text((g_left + g_right) / 2, 1.9, tx["id"],
                    fontsize=label_size, color="#555555", va="bottom", ha="center",
                    family="monospace", clip_on=False)
        elif label_pos == "bottom":
            ax.text((g_left + g_right) / 2, 0.1, tx["id"],
                    fontsize=label_size, color="#555555", va="top", ha="center",
                    family="monospace", clip_on=False)


def draw_signal_track(ax, data, rs, re, ymax, cl, cf, label=None,
                      ymax_pos=None, show_ymax_label=True,
                      ymax_label_size=8, show_yticks=True,
                      show_box=False):
    if data:
        xs, ys = [], []
        for ist, ied, val in data:
            s, e = max(ist, rs), min(ied, re)
            if s < e: xs.extend([s, e]); ys.extend([val, val])
        if xs:
            ax.fill_between(xs, ys, 0, color=cf, alpha=0.8, linewidth=0)
            ax.plot(xs, ys, color=cf, linewidth=0.3, alpha=1.0)
    ax.set_ylim(0, max(ymax * 1.25, 1.0))
    ax.set_xlim(rs, re); ax.set_facecolor("none")
    if label: ax.set_ylabel(label, fontsize=8, fontweight="bold", color="#1A1A1A")
    if show_yticks:
        ax.yaxis.set_major_locator(ticker.MaxNLocator(3, integer=True))
        ax.tick_params(axis="y", labelsize=6, pad=1)
    else:
        ax.set_yticks([])
    ax.tick_params(axis="x", labelbottom=False, length=0)
    if show_box:
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
        ax.spines["top"].set_linewidth(0.4); ax.spines["top"].set_color("#CCCCCC")
        ax.spines["right"].set_linewidth(0.4); ax.spines["right"].set_color("#CCCCCC")
    else:
        for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4); ax.spines["bottom"].set_color("#CCCCCC")
    ax.spines["left"].set_linewidth(0.4); ax.spines["left"].set_color("#CCCCCC")
    ax.axhline(y=0, color="#CCCCCC", linewidth=0.2)
    if show_ymax_label and ymax_pos:
        ax.text(ymax_pos[0], ymax_pos[1], f"[0-{ymax:.0f}]", transform=ax.transAxes,
                fontsize=ymax_label_size, color="#1A1A1A", ha="right", va="top")



def _load_cytoband(path):
    import gzip
    bands = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.strip().split("	")
            if len(parts) < 5: continue
            c = parts[0]
            if c not in bands: bands[c] = []
            bands[c].append((int(parts[1]), int(parts[2]), parts[3], parts[4]))
    return bands

def _gie_stain_color(stain):
    if stain == "gneg": return "#FFFFFF"
    if stain == "gpos33": return "#BFBFBF"
    if stain == "gpos50": return "#808080"
    if stain == "gpos66": return "#575757"
    if stain == "gpos75": return "#404040"
    if stain == "gpos100": return "#000000"
    if stain == "acen": return "#8B0000"
    if stain == "gvar": return "#000000"
    if stain == "stalk": return "#19A7CE"
    return "#E0E0E0"

def draw_ideogram(ax, chrom, gene_start, gene_end, cytoband_file, colors=None,
                  region_start=None, region_end=None, trap_color_top="#E0E0E0",
                  trap_color_bot="#404040", trap_height=2.5, trap_smooth=200,
                  marker_size=0.01):
    if not cytoband_file or not chrom: return
    bands = _load_cytoband(cytoband_file)
    chrom_bands = bands.get(chrom, [])
    if not chrom_bands: return
    chr_len = max(e for _, e, _, _ in chrom_bands)
    non_acen = [(s,e) for s,e,_,st in chrom_bands if st != "acen"]
    chr_start = min(s for s,_ in non_acen) if non_acen else 0
    chr_end = max(e for _,e in non_acen) if non_acen else chr_len
    ax.set_xlim(chr_start, chr_end)
    ax.set_ylim(0, 3.0)
    ax.set_facecolor("none")
    bar_y = 0.8; bar_h = 0.8
    for s, e, name, stain in chrom_bands:
        col = _gie_stain_color(stain)
        ax.add_patch(Rectangle((s, bar_y), e - s, bar_h, facecolor=col,
                     edgecolor="white", linewidth=0.1))
    ax.add_patch(Rectangle((chr_start, bar_y), chr_end - chr_start, bar_h,
                 facecolor="none", edgecolor="#333333", linewidth=0.5))
    marker_y = bar_y - 0.05; marker_h = bar_h + 0.1
    ax.add_patch(Rectangle((gene_start, marker_y), gene_end - gene_start, marker_h,
                 facecolor="#E74C3C", edgecolor="#C0392B", linewidth=0.5, alpha=0.8))
    # Red equilateral triangle at bottom of gene marker (figure coords)
    bbox = ax.get_position()
    gene_cx_frac = ((gene_start + gene_end) / 2 - chr_start) / (chr_end - chr_start)
    tri_cx_fig = bbox.x0 + gene_cx_frac * bbox.width
    tri_top_fig = bbox.y0 + (marker_y / 3.0) * bbox.height
    side = marker_size
    fig_w, fig_h = plt.gcf().get_size_inches()
    tri_h_fig = side * (3 ** 0.5) / 2 * (fig_w / fig_h)
    plt.gcf().patches.append(matplotlib.patches.Polygon(
        [[tri_cx_fig, tri_top_fig],
         [tri_cx_fig - side / 2, tri_top_fig - tri_h_fig],
         [tri_cx_fig + side / 2, tri_top_fig - tri_h_fig]],
        transform=plt.gcf().transFigure,
        facecolor="#E74C3C", edgecolor="#E74C3C", linewidth=0.5, alpha=0.8, closed=True))
    if region_start is not None and region_end is not None:
        bbox = ax.get_position()
        panel_top = bbox.y0 + bbox.height
        chr_top_f = (bar_y + bar_h) / 3.0
        chr_top_y = bbox.y0 + chr_top_f * bbox.height
        gs_f = (gene_start - chr_start) / (chr_end - chr_start)
        ge_f = (gene_end - chr_start) / (chr_end - chr_start)
        left_x = bbox.x0; right_x = bbox.x0 + bbox.width
        bot_x_l = bbox.x0 + gs_f * bbox.width
        bot_x_r = bbox.x0 + ge_f * bbox.width
        r1,g1,b1 = int(trap_color_top[1:3],16),int(trap_color_top[3:5],16),int(trap_color_top[5:7],16)
        r2,g2,b2 = int(trap_color_bot[1:3],16),int(trap_color_bot[3:5],16),int(trap_color_bot[5:7],16)
        for i in range(trap_smooth):
            t = i / trap_smooth
            x_l = left_x + (bot_x_l - left_x) * t
            x_r = right_x + (bot_x_r - right_x) * t
            y0 = panel_top + (chr_top_y - panel_top) * t
            y1 = panel_top + (chr_top_y - panel_top) * (t + 1.0/trap_smooth)
            r,g,b = int(r1+(r2-r1)*t),int(g1+(g2-g1)*t),int(b1+(b2-b1)*t)
            color = f"#{r:02x}{g:02x}{b:02x}"
            plt.gcf().patches.append(matplotlib.patches.Polygon(
                [[x_l,y1],[x_r,y1],[x_r,y0],[x_l,y0]],
                transform=plt.gcf().transFigure, facecolor=color, edgecolor=color,
                linewidth=0, alpha=0.6, closed=True))
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)


def add_highlights(axes, highlights, region_start, region_end, chrom=None):
    if not highlights: return
    for h in highlights:
        hc, hs, he, hcol = h
        if hc is not None and chrom is not None and hc != chrom: continue
        s, e = max(hs, region_start), min(he, region_end)
        if s < e:
            for ax in axes:
                ax.axvspan(s, e, color=hcol, alpha=0.2, linewidth=0, zorder=0)


def draw_xaxis(ax, rs, re, chrom=None):
    ax.set_xlim(rs, re); ax.set_ylim(0, 1); ax.set_facecolor("none")
    if chrom:
        ax.set_xlabel(f"{chrom}:{rs:,}-{re:,}", fontsize=6, color="#555555")
        ax.xaxis.set_ticks([])
    else:
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{x/1e6:.2f} Mb"))
        ax.tick_params(axis="x", labelsize=6, rotation=0, pad=2)
    ax.set_yticks([])
    for sp in ["top", "right", "left"]: ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4); ax.spines["bottom"].set_color("#CCCCCC")


# ── High-level plotting functions ───────────────────────────────

def plot_faceted(genes, data, track_labels, ymax_values, colors, output,
                 title=None, gene_model_bottom=True, figsize=(14, 6.5),
                 show_coords=True, wspace=0.04, utr_ratio=0.5,
                 ymax_override=None, ymax_pos=None, ymax_label_size=8,
                 show_yticks=True, show_box=False, track_colors=None, yscale="gene", gene_ratio=0.8, highlights=None, cytoband=None, cytoband_height=0.6, trap_smooth=200, marker_size=0.01,
                 trap_color_top="#E0E0E0", trap_color_bot="#404040", trap_height=2.5):
    """Create a faceted multi-gene track plot."""
    gene_names = list(genes.keys())
    n_genes = len(gene_names)
    n_tracks = len(track_labels)
    head = 1 if show_coords else 0
    ideo = 1 if cytoband else 0
    tot = n_tracks + 1 + head + ideo
    if gene_model_bottom:
        hr = ([0.35] if show_coords else []) + [2.0] * n_tracks + [gene_ratio] + ([cytoband_height] if cytoband else [])
        gmr = n_tracks + head
        s0 = head
    else:
        hr = ([0.35] if show_coords else []) + [gene_ratio] + [2.0] * n_tracks + ([cytoband_height] if cytoband else [])
        gmr = head
        s0 = head + 1
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(tot, n_genes, hspace=0.05, wspace=wspace,
                          height_ratios=hr,
                          left=0.07, right=0.96, top=0.91, bottom=0.11)

    for gi_idx, gn in enumerate(gene_names):
        gi = genes[gn]; gd = data[gn]; rs, re = gd["region"]; ymax = gd["ymax"]
        chrom = gi["chr"]

        # --- Header: chromosome coordinate row ---
        if show_coords:
            ax_h = fig.add_subplot(gs[0, gi_idx])
            ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1); ax_h.set_facecolor("white")
            ax_h.text(0.5, 0.7, f"{gn} | {chrom}:{rs:,}-{re:,}", transform=ax_h.transAxes,
                      fontsize=8, fontweight="bold", fontstyle="italic",
                      color=colors["highlight"], ha="center", va="center")
            ax_h.set_yticks([]); ax_h.set_xticks([])
            for sp in ax_h.spines.values(): sp.set_visible(False)

        # --- Signal tracks ---
        track_axes = []
        for ti, lbl in enumerate(track_labels):
            ax = fig.add_subplot(gs[s0 + ti, gi_idx]); track_axes.append(ax)
            ym = ymax_override if ymax_override is not None else (gd.get("track_ymax",{}).get(lbl,ymax) if yscale=="track" else ymax)
            cl, cf = _resolve_track_color(str(lbl), track_colors, colors)
            draw_signal_track(ax, gd["tracks"].get(lbl, []), rs, re, ym, cl, cf,
                              str(lbl) if gi_idx == 0 else None,
                              ymax_pos=ymax_pos, ymax_label_size=ymax_label_size,
                              show_yticks=show_yticks, show_box=show_box)

        # --- Gene model ---
        ax_g = fig.add_subplot(gs[gmr, gi_idx]); ax_g.set_xlim(rs, re); ax_g.set_ylim(0, 3); ax_g.set_facecolor("none")
        tx = gi["transcripts"][0] if gi.get("transcripts") else {}
        draw_gene_model(ax_g, tx, rs, re, colors, utr_ratio=utr_ratio,
                        gene_names=tx.get("gene_names"))
        draw_strand_arrow(ax_g, gi["start"], gi["end"], gi.get("strand", "+"), colors)
        if not cytoband: ax_g.set_xlabel(f"{chrom}:{rs:,}-{re:,}", fontsize=6, color="#888888")
        ax_g.xaxis.set_ticks([])
        ax_g.set_yticks([])
        for sp in ["top", "right", "left"]: ax_g.spines[sp].set_visible(False)

        if cytoband:
            ax_ideo = fig.add_subplot(gs[tot - 1, gi_idx])
            draw_ideogram(ax_ideo, chrom, gi["start"], gi["end"], cytoband, colors,
                         region_start=gd["region"][0], region_end=gd["region"][1],
                         trap_smooth=trap_smooth, marker_size=marker_size,
                         trap_color_top=trap_color_top, trap_color_bot=trap_color_bot,
                         trap_height=trap_height)
        if highlights: add_highlights(track_axes + [ax_g], highlights, rs, re, chrom=gi["chr"])

    if title: fig.suptitle(str(title), fontsize=10, fontweight="bold", y=0.997, color="#1A1A1A")
    fig.text(0.07, 0.028, "CDS", fontsize=7, fontweight="bold", color=colors["cds"])
    fig.text(0.12, 0.028, "UTR", fontsize=7, fontweight="bold", color=colors["utr"])
    fig.savefig(str(output), dpi=300, facecolor="white", edgecolor="none")
    plt.close(fig)


def plot_isoforms(genes, data, track_labels, ymax_values, colors, output,
                  title=None, iso_h=0.35, figsize=(15, 8),
                  iso_label_pos="bottom", iso_label_size=3.8,
                  iso_align="top", wspace=0.06,
                  ymax_override=None, ymax_pos=None, ymax_label_size=8,
                 show_yticks=True, show_box=False, track_colors=None, yscale="gene",
                 show_isoform_label=True, show_coords=True, highlights=None, gene_ratio=None,
                 cytoband=None, cytoband_height=0.6, trap_smooth=200, marker_size=0.01,
                 trap_color_top="#E0E0E0", trap_color_bot="#404040", trap_height=2.5):
    """Create track plot with isoform-level gene models."""
    gene_names = list(genes.keys())
    n_genes = len(gene_names)
    n_tracks = len(track_labels)
    max_iso = max(len(genes[gn]["transcripts"]) for gn in gene_names)
    head = 2 if show_coords else 0
    ideo = 1 if cytoband else 0
    tot = n_tracks + max_iso + 2 + head + ideo
    hr = ([0.35, 0.25] if show_coords else []) + [2.0] * n_tracks + [0.15] + [iso_h] * max_iso + [0.3] + ([cytoband_height] if cytoband else [])
    s0 = head

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(tot, n_genes, hspace=0.03, wspace=wspace,
                          height_ratios=hr,
                          left=0.06, right=0.88, top=0.95, bottom=0.05)

    for gi_idx, gn in enumerate(gene_names):
        gi = genes[gn]; gd = data[gn]; rs, re = gd["region"]; ymax = gd["ymax"]
        chrom = gi["chr"]
        txs = gi["transcripts"]; niso = len(txs)

        if show_coords:
            ax_h = fig.add_subplot(gs[0, gi_idx])
            ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1); ax_h.set_facecolor("white")
            ax_h.text(0.5, 0.7, f"{chrom}:{rs:,}-{re:,}", transform=ax_h.transAxes,
                      fontsize=8, fontweight="bold", color=colors["highlight"],
                      ha="center", va="center")
            ax_h.set_yticks([]); ax_h.set_xticks([])
            for sp in ax_h.spines.values(): sp.set_visible(False)
            ax_i = fig.add_subplot(gs[1, gi_idx])
            ax_i.set_xlim(0, 1); ax_i.set_ylim(0, 1); ax_i.set_facecolor("white")
            ax_i.text(0.5, 0.5, f"{gn} | {chrom} | {niso} isoform{'s' if niso>1 else ''}",
                      transform=ax_i.transAxes, fontsize=8, fontweight="bold",
                      fontstyle="italic", color=colors["highlight"],
                      ha="center", va="center")
            ax_i.set_yticks([]); ax_i.set_xticks([])
            for sp in ax_i.spines.values(): sp.set_visible(False)

        track_axes = []
        for ti, lbl in enumerate(track_labels):
            ax = fig.add_subplot(gs[s0 + ti, gi_idx]); track_axes.append(ax)
            ym = ymax_override if ymax_override is not None else (gd.get("track_ymax",{}).get(lbl,ymax) if yscale=="track" else ymax)
            cl, cf = _resolve_track_color(str(lbl), track_colors, colors)
            draw_signal_track(ax, gd["tracks"].get(lbl, []), rs, re, ym, cl, cf,
                              str(lbl) if gi_idx == 0 else None,
                              ymax_pos=ymax_pos, ymax_label_size=ymax_label_size,
                              show_yticks=show_yticks, show_box=show_box)

        # Strand direction arrow (between signals and isoforms)
        ax_strand = fig.add_subplot(gs[n_tracks + head, gi_idx])
        ax_strand.set_xlim(rs, re); ax_strand.set_ylim(0, 3)
        ax_strand.set_facecolor("none")
        draw_strand_arrow(ax_strand, gi["start"], gi["end"], gi.get("strand", "+"), colors, arrow_y=1.5)
        ax_strand.set_yticks([]); ax_strand.set_xticks([])
        for sp in ax_strand.spines.values(): sp.set_visible(False)

        iso_s = n_tracks + head + 1
        if iso_align == "top":
            off_top = 0
        elif iso_align == "center":
            off_top = (max_iso - niso) // 2
        else:  # bottom
            off_top = max_iso - niso
        iso_axes = []
        for ti, tx in enumerate(txs):
            ax = fig.add_subplot(gs[iso_s + off_top + ti, gi_idx])
            iso_axes.append(ax)
            ax.set_xlim(rs, re); ax.set_ylim(0, 2); ax.set_facecolor("none")
            draw_isoform_row(ax, tx, rs, re, colors, iso_label_pos, iso_label_size,
                             show_label=show_isoform_label)
            ax.set_yticks([]); ax.set_xticks([])
            for sp in ax.spines.values(): sp.set_visible(False)

        for er in range(max_iso - niso):
            if er < off_top or er >= off_top + niso:
                ax = fig.add_subplot(gs[iso_s + er, gi_idx])
                ax.set_facecolor("none"); ax.set_yticks([]); ax.set_xticks([])
                for sp in ax.spines.values(): sp.set_visible(False)

        ax_x = fig.add_subplot(gs[tot - 1 - ideo, gi_idx]); draw_xaxis(ax_x, rs, re, chrom=gi["chr"])
        if cytoband:
            ax_ideo = fig.add_subplot(gs[tot - 1, gi_idx])
            draw_ideogram(ax_ideo, gi["chr"], gi["start"], gi["end"], cytoband, colors,
                         region_start=gd["region"][0], region_end=gd["region"][1],
                         trap_smooth=trap_smooth, marker_size=marker_size,
                         trap_color_top=trap_color_top, trap_color_bot=trap_color_bot,
                         trap_height=trap_height)
        if highlights: add_highlights(track_axes + iso_axes + [ax_h, ax_i] if show_coords else track_axes + iso_axes, highlights, rs, re, chrom=gi["chr"])

    if title: fig.suptitle(str(title), fontsize=10, fontweight="bold", y=0.995, color="#1A1A1A")
    fig.text(0.06, 0.015, "CDS", fontsize=7, fontweight="bold", color=colors["cds"])
    fig.text(0.09, 0.015, "UTR", fontsize=7, fontweight="bold", color=colors["utr"])
    fig.text(0.13, 0.015, "Non-coding", fontsize=7, color=colors["noncoding"])
    fig.savefig(str(output), dpi=300, facecolor="white", edgecolor="none")
    plt.close(fig)


def _pack_transcripts_igv(txs):
    """Pack transcripts into rows IGV-style: non-overlapping share a row."""
    if not txs: return []
    sorted_txs = sorted(txs, key=lambda t: min(
        e[0] for e in t.get("exons", [[0, 0]])))
    rows = []
    for tx in sorted_txs:
        tx_exons = tx.get("exons", [])
        if not tx_exons:
            rows.append([tx])
            continue
        tx_start = min(e[0] for e in tx_exons)
        tx_end = max(e[1] for e in tx_exons)
        placed = False
        for row in rows:
            can_place = True
            for existing in row:
                ex = existing.get("exons", [])
                if not ex: continue
                ex_start = min(e[0] for e in ex)
                ex_end = max(e[1] for e in ex)
                if tx_start < ex_end and tx_end > ex_start:
                    can_place = False
                    break
            if can_place:
                row.append(tx)
                placed = True
                break
        if not placed:
            rows.append([tx])
    return rows


def plot_isoforms_regions(regions_data, data, track_labels, ymax_values, colors, output,
                          title=None, iso_h=0.35, figsize=(15, 8),
                          iso_label_pos="bottom", iso_label_size=3.8,
                          iso_align="top", wspace=0.06,
                          ymax_override=None, ymax_pos=None, ymax_label_size=8,
                          show_yticks=True, show_box=False, track_colors=None, yscale="gene",
                          show_isoform_label=True, show_coords=True, highlights=None,
                          cytoband=None, cytoband_height=0.6, trap_smooth=200, marker_size=0.01,
                          trap_color_top="#E0E0E0", trap_color_bot="#404040", trap_height=2.5):
    """Create isoform-level track plot for specified genomic regions.
    Each column = one region, showing all overlapping transcripts packed IGV-style."""
    region_labels = list(regions_data.keys())
    n_regions = len(region_labels)
    n_tracks = len(track_labels)

    # Pre-pack transcripts per region (IGV-style)
    packed_per_region = {}
    for rn in region_labels:
        packed_per_region[rn] = _pack_transcripts_igv(regions_data[rn]["transcripts"])
    max_rows = max(len(packed_per_region[rn]) for rn in region_labels)

    head = 2 if show_coords else 0
    ideo = 1 if cytoband else 0
    tot = n_tracks + max_rows + 1 + head + ideo
    hr = ([0.35, 0.25] if show_coords else []) + [2.0] * n_tracks + [iso_h] * max_rows + [0.3] + ([cytoband_height] if cytoband else [])
    s0 = head

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(tot, n_regions, hspace=0.03, wspace=wspace,
                          height_ratios=hr,
                          left=0.06, right=0.88, top=0.95, bottom=0.05)

    for gi_idx, rn in enumerate(region_labels):
        gi = regions_data[rn]; gd = data[rn]; rs, re = gd["region"]; ymax = gd["ymax"]
        chrom = gi["chr"]
        packed_rows = packed_per_region[rn]
        niso = len(regions_data[rn]["transcripts"])
        nrows = len(packed_rows)

        if show_coords:
            ax_h = fig.add_subplot(gs[0, gi_idx])
            ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1); ax_h.set_facecolor("white")
            ax_h.text(0.5, 0.7, f"{chrom}:{rs:,}-{re:,}", transform=ax_h.transAxes,
                      fontsize=8, fontweight="bold", color=colors["highlight"],
                      ha="center", va="center")
            ax_h.set_yticks([]); ax_h.set_xticks([])
            for sp in ax_h.spines.values(): sp.set_visible(False)
            ax_i = fig.add_subplot(gs[1, gi_idx])
            ax_i.set_xlim(0, 1); ax_i.set_ylim(0, 1); ax_i.set_facecolor("white")
            ax_i.text(0.5, 0.5, f"Region {gi_idx+1} | {chrom} | {niso} tx in {nrows} row{'s' if nrows>1 else ''}",
                      transform=ax_i.transAxes, fontsize=8, fontweight="bold",
                      fontstyle="italic", color=colors["highlight"],
                      ha="center", va="center")
            ax_i.set_yticks([]); ax_i.set_xticks([])
            for sp in ax_i.spines.values(): sp.set_visible(False)

        track_axes = []
        for ti, lbl in enumerate(track_labels):
            ax = fig.add_subplot(gs[s0 + ti, gi_idx]); track_axes.append(ax)
            ym = ymax_override if ymax_override is not None else (gd.get("track_ymax",{}).get(lbl,ymax) if yscale=="track" else ymax)
            cl, cf = _resolve_track_color(str(lbl), track_colors, colors)
            draw_signal_track(ax, gd["tracks"].get(lbl, []), rs, re, ym, cl, cf,
                              str(lbl) if gi_idx == 0 else None,
                              ymax_pos=ymax_pos, ymax_label_size=ymax_label_size,
                              show_yticks=show_yticks, show_box=show_box)

        iso_s = n_tracks + head
        if iso_align == "top":
            off_top = 0
        elif iso_align == "center":
            off_top = (max_rows - nrows) // 2
        else:
            off_top = max_rows - nrows
        iso_axes = []
        for ri, row_txs in enumerate(packed_rows):
            ax = fig.add_subplot(gs[iso_s + off_top + ri, gi_idx])
            iso_axes.append(ax)
            ax.set_xlim(rs, re); ax.set_ylim(0, 2); ax.set_facecolor("none")
            for tx in row_txs:
                tx_display = dict(tx)
                if tx.get("gene_name"):
                    tx_display["id"] = tx["gene_name"]
                draw_isoform_row(ax, tx_display, rs, re, colors, iso_label_pos, iso_label_size,
                                 show_label=show_isoform_label, strand=tx.get("strand"))
            ax.set_yticks([]); ax.set_xticks([])
            for sp in ax.spines.values(): sp.set_visible(False)

        for er in range(max_rows - nrows):
            if er < off_top or er >= off_top + nrows:
                ax = fig.add_subplot(gs[iso_s + er, gi_idx])
                ax.set_facecolor("none"); ax.set_yticks([]); ax.set_xticks([])
                for sp in ax.spines.values(): sp.set_visible(False)

        ax_x = fig.add_subplot(gs[tot - 1 - ideo, gi_idx]); draw_xaxis(ax_x, rs, re, chrom=gi["chr"])
        if cytoband:
            ax_ideo = fig.add_subplot(gs[tot - 1, gi_idx])
            draw_ideogram(ax_ideo, gi["chr"], gi["start"], gi["end"], cytoband, colors,
                         region_start=gd["region"][0], region_end=gd["region"][1],
                         trap_smooth=trap_smooth, marker_size=marker_size,
                         trap_color_top=trap_color_top, trap_color_bot=trap_color_bot,
                         trap_height=trap_height)
        if highlights: add_highlights(track_axes + iso_axes + [ax_h, ax_i] if show_coords else track_axes + iso_axes, highlights, rs, re, chrom=gi["chr"])

    if title: fig.suptitle(str(title), fontsize=10, fontweight="bold", y=0.995, color="#1A1A1A")
    fig.text(0.06, 0.015, "CDS", fontsize=7, fontweight="bold", color=colors["cds"])
    fig.text(0.09, 0.015, "UTR", fontsize=7, fontweight="bold", color=colors["utr"])
    fig.text(0.13, 0.015, "Non-coding", fontsize=7, color=colors["noncoding"])
    fig.savefig(str(output), dpi=300, facecolor="white", edgecolor="none")
    plt.close(fig)
