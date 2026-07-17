# trackPy API Reference (v0.1.5.1)

## CLI Commands

### `trackpy info <file.bw>`
Print chromosome names and sizes.

### `trackpy query <file.bw> <region> [-o out.txt]`
Query values. Region: `chr:start-end`. Outputs tab-separated.

### `trackpy plot <mode> <genes...> [options]`
Three modes: `faceted`, `isoforms`, `regions`.
- `faceted`: gene names or `chr:start-end` regions (auto-detected). Regions collapse all genes per interval.
- `isoforms`: gene names only. Shows all transcripts per gene with strand arrows.
- `regions`: `chr:start-end` only. IGV-style packed transcript rows with intron strand arrows.

---

## Parameters

### Input / Output

| Param | Default | Description |
|-------|---------|-------------|
| `-g, --gtf` | required | GTF or GFF3 (`.gz` supported, auto-detected) |
| `-b, --bw-files` | required | BigWig (`.bw`) or bedGraph (`.bedgraph`, `.bedgraph.gz`), auto-detected |
| `-l, --labels` | filename | Display label per track. Same count as `-b`. |
| `-o, --output` | `trackpy_output` | Output base name. |

### Layout

| Param | Default | Description |
|-------|---------|-------------|
| `--flank-up` | `3000` | bp upstream of gene start |
| `--flank-down` | `3000` | bp downstream of gene end |
| `--wspace` | auto | Horizontal gap between columns. Auto based on label length. Override with float. |
| `--width` | 14 / 15 | Figure width in inches (faceted / isoforms) |
| `--height` | 6.5 / 8 | Figure height in inches |
| `--gene-model-top` | off | Place gene model above signal tracks |
| `--no-coords` | off | Hide coordinate header |
| `--show-box` | off | Show border on all 4 sides of each track |
| `--gene-ratio` | `0.8` | Gene model panel height relative to signal track |

### Gene Model

| Param | Default | Description |
|-------|---------|-------------|
| `--utr-ratio` | `0.5` | UTR height / CDS height (1.0 = equal) |
| `--cds-color` | `#1A1A1A` | CDS fill color. Also used for non-coding exons. |
| `--utr-color` | `#1A1A1A` | UTR fill color |
| `--intron-color` | `#1A1A1A` | Intron line color |

### Isoform

| Param | Default | Description |
|-------|---------|-------------|
| `--isoform-height` | `0.35` | Isoform row height |
| `--isoform-label-pos` | `bottom` | Transcript ID position: `left`, `right`, `top`, `bottom` |
| `--isoform-label-size` | `6` | Transcript ID font size |
| `--no-isoform-label` | off | Hide transcript ID labels |
| `--isoform-align` | `top` | Row alignment within column: `top`, `center`, `bottom` |

### Y-Axis

| Param | Default | Description |
|-------|---------|-------------|
| `--ymax` | auto (99th percentile) | Fixed y-axis ceiling for all tracks |
| `--yscale` | `gene` | `gene`: shared per gene. `track`: independent per track |
| `--ymax-pos` | `0.95 0.95` | Range label position in axes coords |
| `--ymax-label-size` | `8` | Range label font size |
| `--no-range-label` | off | Hide `[0-xxx]` label |
| `--no-yticks` | off | Hide y-axis ticks and values |

### Track Colors

| Param | Default | Description |
|-------|---------|-------------|
| `--track-colors` | auto | One HEX per `-b` file, same order. |

### Highlights

| Param | Default | Description |
|-------|---------|-------------|
| `--highlight` | — | `REGION COLOR`. Region: `chr:start-end` or `start-end`. Repeatable. |

### Chromosome Ideogram

| Param | Default | Description |
|-------|---------|-------------|
| `--cytoband` | — | Path to cytoband file (`.gz` supported). Enables chromosome ideogram below gene model. |
| `--trap-color` | `#E0E0E0 #404040` | Two HEX colors for trapezoid gradient: `TOP_COLOR BOTTOM_COLOR` |
| `--trap-height` | `2.5` | Trapezoid height |
| `--trap-smooth` | `200` | Number of gradient slices for trapezoid. Higher = smoother. |
| `--marker-size` | `0.01` | Red triangle marker size on cytoband (figure fraction). |
| `--cytoband-height` | `0.6` | Chromosome panel height |

When `--cytoband` is set:
- Full chromosome ideogram with IGV-standard cytoband coloring appears below each gene
- Red marker on the chromosome indicates the gene position
- Gray gradient trapezoid above the chromosome shows the zoom relationship (top=panel width, bottom=gene position)
- Gene model x-axis coordinate labels are hidden (redundant with ideogram)

### Zoom (faceted & isoforms, gene-name mode only)

| Param | Default | Description |
|-------|---------|-------------|
| `--zoom-region` | — | Sub-region(s) to magnify. `START-END` pairs, comma-separated. One per gene; single value applies to all genes. E.g. `--zoom-region 10904000-10905000,11006000-11009000` |
| `--zoom-position` | `bottom` | `bottom`: full gene on top, zoom below with trapezoid connecting narrow(zoom region)→wide(full zoom panel). `top`: reversed. |

When `--zoom-region` is set:
- Each gene column splits into full + zoom panels stacked vertically
- A gradient trapezoid connects the panels: top edge at zoom region position (narrow), bottom edge spanning full zoom panel width (wide)
- For `--zoom-position top`, the trapezoid is inverted (top=wide, bottom=narrow)
- Shared trapezoid colors via `--trap-color TOP BOTTOM`

---

## Python API

### `trackpy.BigWigReader(filepath)`
Pure Python bigWig reader. `query(chrom, start, end)` → `[(s,e,v),...]`. `chromosomes` property.

### `trackpy.BedGraphReader(filepath)`
BedGraph reader (supports `.gz`). Same interface.

### `trackpy.parse_annotations(filepath, gene_names)`
Parse GTF or GFF3 (auto-detected). Returns gene structures dict.

### `trackpy.load_gene_data(genes, bw_paths, flank_up=3000, flank_down=3000)`
Load bigWig/bedGraph data. Auto-syncs chr prefix. Returns `{gene: {region, chrom, tracks, ymax, track_ymax}}`.

### `trackpy.parse_regions(filepath, regions)`
Parse transcripts overlapping genomic intervals. `regions`: list of `(chrom, start, end, label)`.

### `trackpy.parse_faceted_regions(filepath, regions)`
Find all genes in regions, collapse transcripts per region.

### `trackpy.plot_faceted(..., cytoband=None)`
### `trackpy.plot_isoforms(...)`
### `trackpy.plot_isoforms_regions(..., cytoband=None)`
IGV-style packed transcript rows for genomic regions. Non-overlapping transcripts share a row.

### `trackpy.IGV_COLORS`

Default color dict. cytoband files provided in `demo/cytoband/` for mm10, mm39, hg19, hg38.
