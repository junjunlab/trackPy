# trackPy User Guide

> **v0.1.5** — Publication-quality genomic track plots from the command line.

---

## Table of Contents

1. [Installation](#1-installation)
2. [Quick Start](#2-quick-start)
3. [Core Concepts](#3-core-concepts)
4. [Commands Reference](#4-commands-reference)
5. [Plot Modes](#5-plot-modes)
6. [Parameters Reference](#6-parameters-reference)
7. [Recipes](#7-recipes)
8. [Python API](#8-python-api)

---

## 1. Installation

```bash
# From GitHub (recommended)
pip install git+https://github.com/junjunlab/trackPy.git

# Or clone and install locally
git clone https://github.com/junjunlab/trackPy.git
cd trackPy
pip install -e .
```

**Requirements:** Python >= 3.9, numpy >= 1.20, matplotlib >= 3.5

---

## 2. Quick Start

### 2.1 Check your data

```bash
# List chromosomes in a bigWig file
trackpy info demo/testdata/GSM5746912_MS_2cell_IP_rep1.bigWig

# Query signal values in a region
trackpy query demo/testdata/GSM5746912_MS_2cell_IP_rep1.bigWig chr7:10900000-10901000
```

**Output:**
```
File: demo/testdata/GSM5746912_MS_2cell_IP_rep1.bigWig
Chromosomes: 22
  chr1: 195,471,971 bp
  chr7: 145,441,459 bp
  ...
```

### 2.2 Your first plot

```bash
trackpy plot faceted Zscan4c Zscan4d Zscan4e Zscan4f \
  -g demo/testdata/Mus_musculus.GRCm38.102.gtf.gz \
  -b demo/testdata/GSM5746910_MS_2cell_Input_rep1.bigWig \
     demo/testdata/GSM5746911_MS_2cell_Input_rep2.bigWig \
     demo/testdata/GSM5746912_MS_2cell_IP_rep1.bigWig \
     demo/testdata/GSM5746913_MS_2cell_IP_rep2.bigWig \
  -l "Input rep1" "Input rep2" "IP rep1" "IP rep2" \
  --cytoband demo/cytoband/mm10_cytoBandIdeo.txt.gz \
  --show-box -o zscan4_family
```

![Faceted gene plot](demo/output/guide_faceted_genes.png)

### 2.3 Region-based plot (IGV-style isoform browser)

```bash
trackpy plot regions chr14:54835580-55001465 \
  -g demo/testdata/Mus_musculus.GRCm38.102.gtf.gz \
  -b demo/testdata/GSM5746910_MS_2cell_Input_rep1.bigWig \
     demo/testdata/GSM5746911_MS_2cell_Input_rep2.bigWig \
     demo/testdata/GSM5746912_MS_2cell_IP_rep1.bigWig \
     demo/testdata/GSM5746913_MS_2cell_IP_rep2.bigWig \
  -l "Input rep1" "Input rep2" "IP rep1" "IP rep2" \
  --show-box -o regions_demo
```

![Regions isoform plot](demo/output/guide_regions_iso.png)

### 2.4 Faceted region view (collapsed gene models per interval)

```bash
trackpy plot faceted chr14:54835580-55001465 chr7:73025897-76116527 \
  -g demo/testdata/Mus_musculus.GRCm38.102.gtf.gz \
  -b demo/testdata/GSM5746910_MS_2cell_Input_rep1.bigWig \
     demo/testdata/GSM5746911_MS_2cell_Input_rep2.bigWig \
     demo/testdata/GSM5746912_MS_2cell_IP_rep1.bigWig \
     demo/testdata/GSM5746913_MS_2cell_IP_rep2.bigWig \
  -l "Input rep1" "Input rep2" "IP rep1" "IP rep2" \
  --show-box -o faceted_regions
```

![Faceted region plot](demo/output/guide_faceted_regions.png)

---

## 3. Core Concepts

### 3.1 Input Files

| Type | Format | Extension | Description |
|------|--------|-----------|-------------|
| **Signal** | bigWig | `.bw`, `.bigwig` | Binary indexed signal track |
| **Signal** | bedGraph | `.bedgraph`, `.bedgraph.gz` | Text-based signal track |
| **Annotation** | GTF | `.gtf`, `.gtf.gz` | Gene/transcript annotations |
| **Annotation** | GFF3 | `.gff3`, `.gff3.gz` | Gene/transcript annotations |
| **Cytoband** | TSV | `.txt.gz` | Chromosome ideogram data |

> All formats are **auto-detected** by file extension. GTF vs GFF3 is also auto-detected.

### 3.2 Plot Modes

```
┌──────────────────────────────────────────────────────┐
│  faceted          isoforms           regions         │
│  ┌──┬──┬──┐      ┌──┬──┬──┐        ┌──┬──┬──┐      │
│  │G1│G2│G3│      │G1│G2│G3│        │R1│R2│R3│      │
│  ├──┼──┼──┤      ├──┼──┼──┤        ├──┼──┼──┤      │
│  │▄▄│▄▄│▄▄│      │▄▄│▄▄│▄▄│        │▄▄│▄▄│▄▄│      │
│  │▄▄│▄▄│▄▄│      │▄▄│▄▄│▄▄│        │▄▄│▄▄│▄▄│      │
│  ├──┼──┼──┤      ├──┼──┼──┤        ├──┼──┼──┤      │
│  │▬▬│▬▬│▬▬│      │▬▬│▬▬│▬▬│        │▬▬│▬▬│▬▬│ ← tx │
│  └──┴──┴──┘      │▬▬│▬▬│▬▬│        │▬▬│   │▬▬│      │
│  1 gene/column    │▬▬│   │▬▬│        │▬▬│   │▬▬│      │
│                   1 gene/column      └──┴──┴──┘      │
│                   all transcripts     1 region/column │
│                                       packed rows     │
└──────────────────────────────────────────────────────┘
```

### 3.3 Auto-Detection in Faceted Mode

The `faceted` mode intelligently detects input type:

| Input Pattern | Mode | Behavior |
|---------------|------|----------|
| `Myc Actb` (no `:`) | Gene-name | One column per gene, canonical transcript |
| `chr7:10Mb-11Mb` (has `:`) | Region | One column per region, all genes collapsed |

---

## 4. Commands Reference

### `trackpy info`
```
trackpy info <file.bw>
```
Print chromosome names and sizes from a bigWig file.

### `trackpy query`
```
trackpy query <file.bw> <region> [-o out.txt]
```
Dump signal values for a genomic region.
- `region`: format `chr:start-end` (e.g., `chr7:10900000-11000000`)
- `-o`: save to file instead of stdout

### `trackpy plot`
```
trackpy plot <mode> <items...> [options]
```
Generate publication-quality PDF track plots.
- `mode`: `faceted` | `isoforms` | `regions`
- `items`: gene names or `chr:start-end` regions

---

## 5. Plot Modes

### 5.1 `faceted` — Multi-Gene/Multi-Region Side-by-Side

Each column = one gene or one genomic region. Gene models shown as single collapsed row.

**Gene-name mode:**
```bash
trackpy plot faceted Myc Jun Actb \
  -g genes.gtf -b input.bw ip.bw -l Input IP -o out
```

**Region mode** (auto-detected when input contains `:`):
```bash
trackpy plot faceted chr14:54835580-55001465 chr7:73025897-76116527 \
  -g genes.gtf -b input.bw ip.bw -l Input IP -o out
```
> In region mode, all genes within each interval are collapsed into one gene model row. Gene names are displayed below the structure.

**Layout (default):**
```
┌──────────────────────────────┐
│  Gene A  │  Gene B  │ Gene C │  ← Header (gene + coordinates)
├──────────┼──────────┼────────┤
│  Signal  │  Signal  │ Signal │  ← Signal tracks
│  Track 1 │  Track 1 │ Track 1│
│  Signal  │  Signal  │ Signal │
│  Track 2 │  Track 2 │ Track 2│
├──────────┼──────────┼────────┤
│  ▬▬▬▬▬  │  ▬▬▬▬▬  │ ▬▬▬▬▬ │  ← Gene model + strand arrow
│  →       │    ←     │   →    │
├──────────┼──────────┼────────┤
│██████████│██████████│████████│  ← Chromosome ideogram (with --cytoband)
└──────────┴──────────┴────────┘
```

**With `--gene-model-top`:**
The gene model row moves above signal tracks:
```
Header → Gene Model → Signal Tracks → Cytoband
```

### 5.2 `isoforms` — All Transcripts per Gene

Each column = one gene. All transcripts shown as individual rows with IDs.

```bash
trackpy plot isoforms Myh6 Myh7 Bcl2l2 Pabpn1 \
  -g genes.gtf -b input.bw ip.bw -l Input IP -o out --show-box
```

![Isoforms gene plot](demo/output/guide_isoforms_genes.png)

### 5.3 `regions` — IGV-Style Region Browser

Each column = one genomic region. Transcripts are **packed IGV-style**: non-overlapping transcripts share the same row to maximize space efficiency. Strand arrows on intron lines indicate transcription direction.

```bash
trackpy plot regions chr14:54835580-55001465 chr19:5790000-5810000 \
  -g genes.gtf -b input.bw ip.bw -l Input IP -o out --show-box
```

**Packing example:**
```
Before (72 transcripts):    After IGV packing (12 rows):
┌────┐                      ┌────┐
│ tx1│                      │tx1 │  ← non-overlapping
│ tx2│                      │ tx3│
│ tx3│                      │tx2 │  ← overlaps tx1, tx3
│ tx4│                      │tx4 │  ← non-overlapping with tx5
│ .. │                      │ .. │
│tx72│                      │tx72│
└────┘                      └────┘
 72 rows → 12 rows (83% reduction)
```

**Label display:** Shows gene name only (not transcript ID).

---

## 6. Parameters Reference

### Input / Output
| Flag | Default | Description |
|------|---------|-------------|
| `-g, --gtf` | **required** | Annotation file (GTF/GFF3, `.gz` OK) |
| `-b, --bw-files` | **required** | Signal files (`.bw`, `.bedgraph`, `.bedgraph.gz`) |
| `-l, --labels` | filename | Display label per track |
| `-o, --output` | `trackpy_output` | Output PDF base name |

### Layout
| Flag | Default | Description |
|------|---------|-------------|
| `--flank-up` | `3000` | bp upstream padding |
| `--flank-down` | `3000` | bp downstream padding |
| `--width` | `14` / `15` | Figure width in inches (faceted / isoforms&regions) |
| `--height` | `6.5` / `8` | Figure height in inches |
| `--gene-model-top` | off | Place gene model above signal tracks |
| `--no-coords` | off | Hide coordinate header row |
| `--show-box` | off | Show border on all 4 sides of each track |
| `--gene-ratio` | `0.8` | Gene model row height relative to signal track |
| `--wspace` | auto | Horizontal gap between columns |

### Gene Model
| Flag | Default | Description |
|------|---------|-------------|
| `--utr-ratio` | `0.5` | UTR height / CDS height |
| `--cds-color` | `#1A1A1A` | CDS fill color |
| `--utr-color` | `#1A1A1A` | UTR fill color |
| `--intron-color` | `#1A1A1A` | Intron line color |

### Isoform
| Flag | Default | Description |
|------|---------|-------------|
| `--isoform-height` | `0.35` | Row height (smaller = more compact) |
| `--isoform-label-pos` | `bottom` | Transcript label: `left` / `right` / `top` / `bottom` |
| `--isoform-label-size` | `6` | Label font size |
| `--no-isoform-label` | off | Hide transcript labels |
| `--isoform-align` | `top` | Row alignment: `top` / `center` / `bottom` |

### Y-Axis
| Flag | Default | Description |
|------|---------|-------------|
| `--ymax` | auto (99%) | Fixed y-axis ceiling |
| `--yscale` | `gene` | `gene` (shared per gene) or `track` (independent) |
| `--ymax-pos` | `0.95 0.95` | Range label position in axes coords |
| `--ymax-label-size` | `8` | Range label font size |
| `--no-range-label` | off | Hide `[0-xxx]` labels |
| `--no-yticks` | off | Hide y-axis ticks |

### Cytoband (Chromosome Ideogram)
| Flag | Default | Description |
|------|---------|-------------|
| `--cytoband` | — | Path to cytoband file (`.gz` OK) |
| `--trap-color` | `#E0E0E0 #404040` | Trapezoid gradient: TOP BOTTOM |
| `--trap-height` | `2.5` | Trapezoid height |
| `--trap-smooth` | `200` | Gradient steps (higher = smoother) |
| `--marker-size` | `0.01` | Red triangle marker size |
| `--cytoband-height` | `0.6` | Chromosome panel height |

> When `--cytoband` is enabled:
> - Full chromosome ideogram with IGV-standard Giemsa staining appears below each gene
> - Red rectangle + triangle marks the gene position
> - Gray gradient trapezoid shows the zoom relationship (panel width → gene position)

### Other
| Flag | Default | Description |
|------|---------|-------------|
| `--track-colors` | auto | One HEX color per `-b` file |
| `--highlight` | — | `REGION COLOR` (repeatable) |

---

## 7. Recipes

### 7.1 Basic gene panel with ideogram
```bash
trackpy plot faceted Myc Jun Actb \
  -g genes.gtf -b input.bw ip.bw -l Input IP \
  --cytoband mm10_cytoBandIdeo.txt.gz --show-box -o panel
```

### 7.2 lncRNA isoform browser
```bash
trackpy plot isoforms Malat1 Xist H19 Airn Hotair \
  -g genes.gtf -b input.bw ip.bw -l Input IP \
  --show-box --isoform-label-pos left -o lncRNA
```

### 7.3 Explore a genomic interval (all transcripts)
```bash
trackpy plot regions chr14:54835580-55001465 \
  -g genes.gtf -b input.bw ip.bw -l Input IP \
  --show-box --cytoband mm10_cytoBandIdeo.txt.gz -o my_region
```

### 7.4 Faceted region comparison
```bash
trackpy plot faceted chr14:54835580-55001465 chr7:73025897-76116527 \
  -g genes.gtf -b input.bw ip.bw -l Input IP \
  --show-box --width 10 -o compare_regions
```

### 7.5 Custom track colors
```bash
trackpy plot faceted Myc Jun -g genes.gtf -b wt.bw ko.bw -l WT KO \
  --track-colors "#3498DB" "#E74C3C" -o custom_colors
```

### 7.6 Highlight specific regions
```bash
trackpy plot faceted Myc -g genes.gtf -b input.bw ip.bw -l Input IP \
  --highlight chr15:61900000-61901000 "#FF000020" \
  --highlight chr15:61902000-61903000 "#0000FF20" -o highlights
```

### 7.7 BedGraph (ATAC-seq) input
```bash
trackpy plot faceted Actb Myc -g genes.gff3.gz \
  -b wt.bedgraph.gz ko.bedgraph.gz -l WT KO \
  --track-colors "#3498DB" "#E74C3C" -o atac
```

### 7.8 Compact layout (no y-ticks, no range labels)
```bash
trackpy plot faceted Myc Jun Actb -g genes.gtf -b a.bw b.bw -l A B \
  --no-yticks --no-range-label --cytoband mm10_cytoBandIdeo.txt.gz -o compact
```

### 7.9 Gene model on top
```bash
trackpy plot faceted Myc Jun -g genes.gtf -b a.bw b.bw -l A B \
  --gene-model-top --cytoband mm10_cytoBandIdeo.txt.gz --show-box -o top
```

### 7.9 Gene model on top

```bash
trackpy plot faceted Zscan4b Zscan4c Zscan4d -g genes.gtf \
  -b input1.bw input2.bw ip1.bw ip2.bw -l In1 In2 IP1 IP2 \
  --gene-model-top --cytoband mm10_cytoBandIdeo.txt.gz --show-box -o top
```

![Gene model on top](demo/output/guide_gene_model_top.png)

### 7.10 lncRNA panel

```bash
trackpy plot faceted Malat1 Xist H19 Airn Hotair \
  -g genes.gtf -b in1.bw in2.bw ip1.bw ip2.bw -l In1 In2 IP1 IP2 \
  --cytoband mm10_cytoBandIdeo.txt.gz --show-box -o lncRNA
```

![lncRNA faceted](demo/output/guide_lncRNA_faceted.png)

---

## 8. Python API

```python
from trackpy import (
    BigWigReader, BedGraphReader,
    parse_annotations, load_gene_data,
    parse_regions, parse_faceted_regions,
    plot_faceted, plot_isoforms, plot_isoforms_regions,
    IGV_COLORS,
)

# ── Gene-name based ──────────────────────────────
genes = parse_annotations("genes.gtf.gz", ["Zscan4b", "Myc"])
data = load_gene_data(genes, {"Input": "in.bw", "IP": "ip.bw"})

plot_faceted(genes, data, ["Input", "IP"], data, IGV_COLORS,
             "faceted.pdf", cytoband="mm10_cytoBandIdeo.txt.gz", show_box=True)

plot_isoforms(genes, data, ["Input", "IP"], data, IGV_COLORS,
              "isoforms.pdf", show_box=True)

# ── Region-based ─────────────────────────────────
regions = [("7", 10900000, 11000000, "chr7:10.9-11.0Mb"),
           ("19", 5790000, 5810000, "chr19:5.79-5.81Mb")]

regions_data = parse_regions("genes.gtf.gz", regions)
plot_isoforms_regions(regions_data, data, ["Input", "IP"], {}, IGV_COLORS,
                      "regions.pdf", show_box=True)

faceted_data = parse_faceted_regions("genes.gtf.gz", regions)
plot_faceted(faceted_data, data, ["Input", "IP"], {}, IGV_COLORS,
             "faceted_regions.pdf", show_box=True)

# ── Low-level I/O ────────────────────────────────
with BigWigReader("signal.bw") as bw:
    print(bw.chromosomes)
    values = bw.query("chr7", 10900000, 11000000)

reader = BedGraphReader("signal.bedgraph.gz")
values = reader.query("chr7", 10900000, 11000000)
```

---

## License

MIT — [https://github.com/junjunlab/trackPy](https://github.com/junjunlab/trackPy)
