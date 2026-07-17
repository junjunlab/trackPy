"""trackPy -- Publication-quality genomic track plots via CLI."""

__version__ = "0.1.5.2"

from trackpy.bigwig import BigWigReader, BedGraphReader
from trackpy.plot import (
    IGV_COLORS,
    plot_faceted,
    plot_isoforms,
    plot_isoforms_regions,
)
from trackpy.plot_zoom import plot_faceted_zoom, plot_isoforms_zoom, plot_regions_zoom
from trackpy.gtf import parse_gtf, parse_gff3, parse_annotations, load_gene_data, parse_regions, parse_faceted_regions
