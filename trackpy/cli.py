#!/usr/bin/env python
"""trackPy CLI — publication-quality genomic track plots.

Usage:
    trackpy query  <file.bw> <region>                 Query bigWig values
    trackpy info   <file.bw>                          Show chromosome info
    trackpy plot   faceted <genes> -g <gtf> -b <bw...> [-l <labels>] Faceted plot
    trackpy plot   isoforms <genes> -g <gtf> -b <bw...> [-l <labels>] Isoform-level plot
    trackpy plot   regions <chr:start-end ...> -g <gtf> -b <bw...> [-l <labels>] Region-based isoform plot

Examples:
    trackpy info data/input.bw
    trackpy query data/input.bw chr7:10900000-10910000
    trackpy plot faceted Zscan4a Zscan4b -g genes.gtf -b a.bw b.bw c.bw d.bw -l Input1 Input2 IP1 IP2 -o out
    trackpy plot isoforms Myc Jun -g genes.gtf -b *.bw -l In1 In2 IP1 IP2 -o out
    trackpy plot regions chr7:10900000-11000000 chr7:14200000-14300000 -g genes.gtf -b in.bw ip.bw -l Input IP -o out
"""

import argparse
import sys
import os

# Ensure local trackpy is importable when running from source
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from trackpy.bigwig import BigWigReader
from trackpy.gtf import parse_annotations, load_gene_data, parse_regions, parse_faceted_regions
from trackpy.plot import (IGV_COLORS, plot_faceted, plot_isoforms, plot_isoforms_regions)
import numpy as np


def cmd_info(args):
    """Show chromosome info for a bigWig file."""
    with BigWigReader(args.bw) as bw:
        chroms = bw.chromosomes
    print(f"File: {args.bw}")
    print(f"Chromosomes: {len(chroms)}")
    for name in sorted(chroms.keys()):
        print(f"  {name}: {chroms[name]:,} bp")


def cmd_query(args):
    """Query a bigWig file for a genomic region."""
    chrom, _, coords = args.region.partition(":")
    start, _, end = coords.partition("-")
    start = int(start.replace(",", ""))
    end = int(end.replace(",", ""))

    with BigWigReader(args.bw) as bw:
        data = bw.query(chrom, start, end)

    if args.out:
        with open(args.out, "w") as f:
            f.write("start\tend\tvalue\n")
            for s, e, v in data:
                f.write(f"{s}\t{e}\t{v:.4f}\n")
        print(f"Wrote {len(data)} rows to {args.out}")
    else:
        for s, e, v in data:
            print(f"{s}\t{e}\t{v:.4f}")


def cmd_plot(args):
    """Create track plots (faceted, isoforms, or regions)."""
    # Organize bigWig files
    bw_paths = {}
    if args.labels:
        labels = args.labels
        if len(labels) != len(args.bw_files):
            print("ERROR: --labels count must match --bw-files count")
            sys.exit(1)
    else:
        labels = []
        for fp in args.bw_files:
            basename = os.path.splitext(os.path.basename(fp))[0]
            # Simplify common prefixes: GSM5746912_MS_2cell_IP_rep1 -> IP rep1
            basename = basename.replace("_", " ")
            labels.append(basename)

    for lbl, fp in zip(labels, args.bw_files):
        # Ensure unique keys
        key = lbl
        i = 1
        while key in bw_paths:
            i += 1
            key = f"{lbl} ({i})"
        bw_paths[key] = fp

    track_labels = list(bw_paths.keys())
    if not track_labels:
        print("ERROR: No bigWig files specified. Use -b.")
        sys.exit(1)

    colors = dict(IGV_COLORS)
    colors["cds"] = args.cds_color
    colors["utr"] = args.utr_color
    colors["intron"] = args.intron_color
    track_colors = None
    if args.track_colors:
        if len(args.track_colors) != len(args.bw_files):
            print(f"ERROR: --track-colors count ({len(args.track_colors)}) must match -b count ({len(args.bw_files)})")
            sys.exit(1)
        track_colors = dict(zip(track_labels, args.track_colors))

    out_base = args.output
    if out_base is None:
        out_base = "trackpy_output"
    out_dir = os.path.dirname(os.path.abspath(out_base)) or "."
    os.makedirs(out_dir, exist_ok=True)

    # Auto-compute wspace from y-axis labels if not set
    if args.wspace is not None:
        wspace = args.wspace
    else:
        if args.no_yticks:
            wspace = 0.02
        else:
            max_label_len = max((len(str(lbl)) for lbl in track_labels), default=5)
            wspace = max(0.06, 0.02 + 0.012 * max_label_len)

    # Parse highlights: "chr7:start-end" or "start-end"
    highlights = None
    if args.highlight:
        highlights = []
        for region_str, color in args.highlight:
            chrom = None
            coords = region_str
            if ":" in region_str:
                chrom, _, coords = region_str.partition(":")
            s, _, e = coords.partition("-")
            highlights.append((chrom, int(s.replace(",", "")), int(e.replace(",", "")), color))

    # Figure size
    default_sizes = {"faceted": (14, 6.5), "isoforms": (15, 8), "regions": (15, 8)}
    dw, dh = default_sizes[args.mode]
    figsize = (args.width or dw, args.height or dh)

    # Title
    title = f"{', '.join(args.genes)} | trackPy | mm10"

    if args.mode == "faceted":
        # Auto-detect: region (chr:start-end) or gene names
        is_region_input = any(":" in g for g in args.genes)
        if is_region_input:
            # Parse regions and discover genes
            parsed_regions = []
            for region_str in args.genes:
                chrom, _, coords = region_str.partition(":")
                start, _, end = coords.partition("-")
                if not chrom or not start or not end:
                    print(f"ERROR: invalid region format '{region_str}'")
                    sys.exit(1)
                chrom = chrom.strip()
                start = int(start.replace(",", ""))
                end = int(end.replace(",", ""))
                label = f"{chrom}:{start:,}-{end:,}"
                parsed_regions.append((chrom, start, end, label))
            print(f"  Searching for genes in {len(parsed_regions)} region(s)...")
            genes = parse_faceted_regions(args.gtf, parsed_regions)
            if not genes:
                print("ERROR: no genes found in specified regions")
                sys.exit(1)
            print(f"  Found {len(genes)} region(s) with transcripts")
        else:
            # Parse gene annotations by name
            genes = parse_annotations(args.gtf, args.genes)
            missing = [gn for gn in args.genes if not genes[gn].get("chr")]
            if missing:
                print(f"Warning: genes not found in annotation, skipped: {missing}")
                genes = {gn: gi for gn, gi in genes.items() if gi.get("chr")}
            if not genes:
                print("ERROR: no valid genes found in annotation")
                sys.exit(1)

        print(f"Loading data for {len(genes)} genes x {len(track_labels)} tracks...")
        data = load_gene_data(genes, bw_paths,
                              flank_up=args.flank_up, flank_down=args.flank_down)
        ymax_values = {gn: data[gn]["ymax"] for gn in genes}

        out = f"{out_base}.pdf"
        plot_faceted(genes, data, track_labels, ymax_values, colors, out,
                     title=title, gene_model_bottom=not args.gene_model_top,
                     show_coords=not args.no_coords, wspace=wspace,
                     utr_ratio=args.utr_ratio,
                     ymax_override=args.ymax,
                     ymax_pos=tuple(args.ymax_pos) if not args.no_range_label else None,
                     ymax_label_size=args.ymax_label_size,
                     show_yticks=not args.no_yticks,
                show_box=args.show_box,
                track_colors=track_colors,
                figsize=figsize,
                yscale=args.yscale, gene_ratio=args.gene_ratio, highlights=highlights, cytoband=args.cytoband, cytoband_height=args.cytoband_height,
                trap_smooth=args.trap_smooth, marker_size=args.marker_size,
                trap_color_top=args.trap_color[0], trap_color_bot=args.trap_color[1],
                trap_height=args.trap_height)
        print(f"  Saved: {out}")

    elif args.mode == "isoforms":
        genes = parse_annotations(args.gtf, args.genes)
        missing = [gn for gn in args.genes if not genes[gn].get("chr")]
        if missing:
            print(f"Warning: genes not found in annotation, skipped: {missing}")
            genes = {gn: gi for gn, gi in genes.items() if gi.get("chr")}
        if not genes:
            print("ERROR: no valid genes found in annotation")
            sys.exit(1)

        print(f"Loading data for {len(genes)} genes x {len(track_labels)} tracks...")
        data = load_gene_data(genes, bw_paths,
                              flank_up=args.flank_up, flank_down=args.flank_down)
        ymax_values = {gn: data[gn]["ymax"] for gn in genes}

        out = f"{out_base}.pdf"
        plot_isoforms(genes, data, track_labels, ymax_values, colors, out,
                      title=title, iso_h=args.isoform_height,
                      iso_label_pos=args.isoform_label_pos,
                      iso_label_size=args.isoform_label_size,
                      iso_align=args.isoform_align, wspace=wspace,
                      show_isoform_label=not args.no_isoform_label,
                      show_coords=not args.no_coords,
                      ymax_override=args.ymax,
                      ymax_pos=tuple(args.ymax_pos) if not args.no_range_label else None,
                      ymax_label_size=args.ymax_label_size,
                      show_yticks=not args.no_yticks,
                show_box=args.show_box,
                track_colors=track_colors,
                figsize=figsize,
                yscale=args.yscale, gene_ratio=args.gene_ratio, highlights=highlights, cytoband=args.cytoband, cytoband_height=args.cytoband_height,
                trap_smooth=args.trap_smooth, marker_size=args.marker_size,
                trap_color_top=args.trap_color[0], trap_color_bot=args.trap_color[1],
                trap_height=args.trap_height)
        print(f"  Saved: {out}")

    elif args.mode == "regions":
        # Parse region strings: chr:start-end
        parsed_regions = []
        for region_str in args.genes:
            chrom, _, coords = region_str.partition(":")
            start, _, end = coords.partition("-")
            if not chrom or not start or not end:
                print(f"ERROR: invalid region format '{region_str}', expected chr:start-end")
                sys.exit(1)
            chrom_clean = chrom.strip()
            start = int(start.replace(",", ""))
            end = int(end.replace(",", ""))
            label = f"{chrom_clean}:{start:,}-{end:,}"
            parsed_regions.append((chrom_clean, start, end, label))

        if not parsed_regions:
            print("ERROR: no valid regions specified")
            sys.exit(1)

        print(f"  Detected GTF format: {args.gtf}")
        print(f"  Scanning annotations for {len(parsed_regions)} regions...")
        regions_data = parse_regions(args.gtf, parsed_regions)

        for label, rd in regions_data.items():
            n_tx = len(rd["transcripts"])
            print(f"  {label}: {n_tx} transcript{'s' if n_tx != 1 else ''}")

        total_tx = sum(len(rd["transcripts"]) for rd in regions_data.values())
        if total_tx == 0:
            print("ERROR: no transcripts found in specified regions")
            sys.exit(1)

        print(f"Loading data for {len(regions_data)} regions x {len(track_labels)} tracks...")
        data = load_gene_data(regions_data, bw_paths,
                              flank_up=args.flank_up, flank_down=args.flank_down)
        ymax_values = {rn: data[rn]["ymax"] for rn in regions_data}

        out = f"{out_base}.pdf"
        plot_isoforms_regions(regions_data, data, track_labels, ymax_values, colors, out,
                              title=title, iso_h=args.isoform_height,
                              iso_label_pos=args.isoform_label_pos,
                              iso_label_size=args.isoform_label_size,
                              iso_align=args.isoform_align, wspace=wspace,
                              show_isoform_label=not args.no_isoform_label,
                              show_coords=not args.no_coords,
                              ymax_override=args.ymax,
                              ymax_pos=tuple(args.ymax_pos) if not args.no_range_label else None,
                              ymax_label_size=args.ymax_label_size,
                              show_yticks=not args.no_yticks,
                              show_box=args.show_box,
                              track_colors=track_colors,
                              figsize=figsize,
                              yscale=args.yscale, highlights=highlights, cytoband=args.cytoband, cytoband_height=args.cytoband_height,
                              trap_smooth=args.trap_smooth, marker_size=args.marker_size,
                              trap_color_top=args.trap_color[0], trap_color_bot=args.trap_color[1],
                              trap_height=args.trap_height)
        print(f"  Saved: {out}")


def main():
    parser = argparse.ArgumentParser(
        description="trackPy — publication-quality genomic track plots")
    sub = parser.add_subparsers(dest="command")

    # info
    p = sub.add_parser("info", help="Show bigWig chromosome info")
    p.add_argument("bw", help="bigWig file")
    p.set_defaults(func=cmd_info)

    # query
    p = sub.add_parser("query", help="Query bigWig values")
    p.add_argument("bw", help="bigWig file")
    p.add_argument("region", help="Region, e.g. chr7:10000000-11000000")
    p.add_argument("-o", "--out", help="Output file (default: stdout)")
    p.set_defaults(func=cmd_query)

    # plot
    p = sub.add_parser("plot", help="Create track plots")
    p.add_argument("mode", choices=["faceted", "isoforms", "regions"])
    p.add_argument("genes", nargs="+", help="Gene names or chr:start-end regions")
    p.add_argument("-g", "--gtf", required=True,
                   help="Annotation file (GTF or GFF3, .gz supported, auto-detected)")
    p.add_argument("-o", "--output", help="Output file base name")
    p.add_argument("-b", "--bw-files", nargs="+", required=True,
                   help="BigWig or bedGraph files (.bw, .bedgraph, .bedgraph.gz, auto-detected)")
    p.add_argument("-l", "--labels", nargs="+",
                   help="Display labels for each bigWig file (default: derived from filename)")
    p.add_argument("--gene-model-top", action="store_true",
                   help="Place gene model at top (default: bottom)")
    p.add_argument("--no-coords", action="store_true",
                   help="Hide chromosome coordinate header row")
    p.add_argument("--isoform-height", type=float, default=0.35,
                   help="Isoform row height, smaller = more compact (default: 0.35)")
    p.add_argument("--isoform-label-pos", choices=["left","right","top","bottom"],
                   default="bottom",
                   help="Position of transcript ID label (default: bottom)")
    p.add_argument("--isoform-label-size", type=float, default=6,
                   help="Font size for transcript ID label (default: 6)")
    p.add_argument("--no-isoform-label", action="store_true",
                   help="Hide transcript ID labels on isoform rows")
    p.add_argument("--isoform-align", choices=["top","center","bottom"],
                   default="top",
                   help="Vertical alignment of isoform rows within column (default: top)")
    p.add_argument("--flank-up", type=int, default=3000,
                   help="Flanking bp upstream of gene start (default: 3000)")
    p.add_argument("--flank-down", type=int, default=3000,
                   help="Flanking bp downstream of gene end (default: 3000)")
    p.add_argument("--wspace", type=float, default=None,
                   help="Horizontal spacing between columns (default: auto based on y-axis labels)")
    p.add_argument("--utr-ratio", type=float, default=0.5,
                   help="UTR/CDS height ratio in gene model (default: 0.5)")
    p.add_argument("--ymax", type=float, default=None,
                   help="Fixed y-axis maximum for all tracks (default: auto)")
    p.add_argument("--yscale", choices=["gene", "track"], default="gene",
                   help="Y-axis scale: gene (all tracks share scale per gene) or track (each track independent)")
    p.add_argument("--ymax-pos", type=float, nargs=2, default=[0.95, 0.95],
                   metavar=("X", "Y"),
                   help="Position of ymax label in axes coords (default: 0.95 0.95 top-right)")
    p.add_argument("--ymax-label-size", type=float, default=8,
                   help="Font size of ymax label (default: 8)")
    p.add_argument("--no-range-label", action="store_true",
                   help="Hide [0-xxx] y-axis range label")
    p.add_argument("--no-yticks", action="store_true",
                   help="Hide y-axis tick labels and values")
    p.add_argument("--show-box", action="store_true",
                   help="Show top and right border lines around each track")
    p.add_argument("--gene-ratio", type=float, default=0.8,
                   help="Gene model panel height relative to signal track (default: 0.8)")
    p.add_argument("--highlight", nargs=2, action="append", default=None,
                   metavar=("REGION", "COLOR"),
                   help="Highlight a genomic region, e.g. --highlight chr7:10901000-10902000 '#FF000020'. Repeat for multiple.")
    p.add_argument("--width", type=float, default=None,
                   help="Figure width in inches (default: 14 faceted, 15 isoforms)")
    p.add_argument("--height", type=float, default=None,
                   help="Figure height in inches (default: 6.5 faceted, 8 isoforms)")
    p.add_argument("--track-colors", nargs="+", default=None,
                   metavar="HEX",
                   help="Colors for each track (one per -b file), e.g. --track-colors #808080 #C0392B")
    p.add_argument("--cds-color", default="#1A1A1A",
                   help="CDS color (default: #1A1A1A)")
    p.add_argument("--utr-color", default="#1A1A1A",
                   help="UTR color (default: #1A1A1A)")
    p.add_argument("--intron-color", default="#1A1A1A",
                   help="Intron line color (default: #1A1A1A)")
    p.add_argument("--cytoband", default=None,
                   help="Path to cytoband file (.gz supported) for chromosome ideogram")
    p.add_argument("--trap-color", nargs=2, default=["#E0E0E0", "#404040"],
                   metavar=("TOP_COLOR", "BOTTOM_COLOR"),
                   help="Gradient colors for trapezoid (default: #E0E0E0 #404040)")
    p.add_argument("--trap-height", type=float, default=2.5,
                   help="Height of the trapezoid zoom indicator (default: 2.5)")
    p.add_argument("--trap-smooth", type=int, default=200,
                   help="Number of gradient slices for trapezoid, higher = smoother (default: 200)")
    p.add_argument("--marker-size", type=float, default=0.01,
                   help="Size of the red triangle marker below gene on cytoband (default: 0.01)")
    p.add_argument("--cytoband-height", type=float, default=0.6,
                   help="Height of the chromosome ideogram panel (default: 0.6)")
    p.set_defaults(func=cmd_plot)

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()
