"""trackPy -- Publication-quality genomic track plots via CLI."""

__version__ = "0.1.4"

from trackpy.bigwig import BigWigReader, BedGraphReader
from trackpy.plot import (
    IGV_COLORS,
    plot_faceted,
    plot_isoforms,
)
from trackpy.gtf import parse_gtf, parse_gff3, parse_annotations, load_gene_data
