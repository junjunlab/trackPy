"""GTF and GFF3 gene structure parser."""

import re
import gzip
import os


def parse_gtf(gtf_file, gene_names):
    """Parse gene structures from GTF. Returns dict gene_name -> info."""
    target = set(gene_names)
    result = {gn: {"chr": None, "start": None, "end": None, "strand": None,
                   "transcripts": {}} for gn in target}

    opener = gzip.open if gtf_file.endswith(".gz") else open
    with opener(gtf_file, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attrs = {}
            for m in re.finditer(r'(\S+)\s+"([^"]*)"', parts[8]):
                attrs[m.group(1)] = m.group(2)
            gn = attrs.get("gene_name", "")
            if gn not in target:
                continue

            feat, chrom = parts[2], parts[0]
            start, end = int(parts[3]), int(parts[4])
            strand, tid = parts[6], attrs.get("transcript_id", "")

            if feat == "gene":
                result[gn].update(chr=chrom, start=start, end=end, strand=strand)
                continue
            if not tid:
                continue

            if tid not in result[gn]["transcripts"]:
                result[gn]["transcripts"][tid] = {
                    "id": tid,
                    "biotype": attrs.get("transcript_biotype", ""),
                    "exons": [], "CDS": [], "UTR_5": [], "UTR_3": [],
                    "is_coding": False,
                }
            tx = result[gn]["transcripts"][tid]
            if feat == "exon": tx["exons"].append([start, end])
            elif feat == "CDS": tx["CDS"].append([start, end]); tx["is_coding"] = True
            elif feat == "five_prime_utr": tx["UTR_5"].append([start, end])
            elif feat == "three_prime_utr": tx["UTR_3"].append([start, end])

    def _clean(iv):
        if not iv: return []
        u = sorted({(a, b) for a, b in iv})
        m = [list(u[0])]
        for s, e in u[1:]:
            if s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
            else: m.append([s, e])
        return m

    for gn in target:
        for tid, tx in result[gn]["transcripts"].items():
            for field in ["exons", "CDS", "UTR_5", "UTR_3"]:
                tx[field] = _clean(tx[field])
        def _tx_score(tx):
            # Prefer: coding, has both UTRs, more exons, longer CDS
            has_both_utr = (1 if (tx["UTR_5"] and tx["UTR_3"]) else 0)
            n_exons = len(tx["exons"])
            cds_len = sum(e - s for s, e in tx["CDS"]) if tx["CDS"] else 0
            return (not tx["is_coding"], -has_both_utr, -n_exons, -cds_len, tx["id"])

        result[gn]["transcripts"] = sorted(
            result[gn]["transcripts"].values(), key=_tx_score)
    return result


def parse_gff3(gff3_file, gene_names):
    """Parse gene structures from a GFF3 file (supports .gz compression).
    Uses two-pass: first find gene IDs by Name, then collect children by Parent."""
    target = set(gene_names)
    result = {gn: {"chr": None, "start": None, "end": None, "strand": None,
                   "transcripts": {}} for gn in target}

    opener = gzip.open if gff3_file.endswith(".gz") else open

    def _parse_attrs(attr_str):
        d = {}
        for part in attr_str.split(";"):
            part = part.strip()
            if "=" in part:
                k, v = part.split("=", 1)
                d[k] = v
        return d

    # Pass 1: find gene IDs matching target names
    gene_id_to_name = {}  # "gene:ENSMUSG..." -> "Myc"
    with opener(gff3_file, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            if parts[2] != "gene": continue
            attrs = _parse_attrs(parts[8])
            gn = attrs.get("Name", "")
            if gn in target:
                gid = attrs.get("ID", "")
                gene_id_to_name[gid] = gn
                result[gn].update(chr=parts[0], start=int(parts[3]),
                                  end=int(parts[4]), strand=parts[6])

    if not gene_id_to_name:
        return result

    # Pass 2: collect transcripts and features by Parent
    with opener(gff3_file, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            attrs = _parse_attrs(parts[8])
            parent = attrs.get("Parent", "")
            if not parent: continue

            feat = parts[2]
            start = int(parts[3])
            end = int(parts[4])

            if feat == "mRNA":
                # Parent is gene ID
                gn = gene_id_to_name.get(parent, "")
                if not gn: continue
                tid_raw = attrs.get("ID", "")
                tid = tid_raw.replace("transcript:", "") if tid_raw.startswith("transcript:") else tid_raw
                if not tid: continue
                result[gn]["transcripts"][tid] = {
                    "id": tid, "biotype": attrs.get("biotype", ""),
                    "exons": [], "CDS": [], "UTR_5": [], "UTR_3": [],
                    "is_coding": False,
                }
            elif feat in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
                # Parent is transcript ID; find which gene it belongs to
                tid = parent.replace("transcript:", "") if parent.startswith("transcript:") else parent
                # Find the gene that contains this transcript
                found_gn = None
                for gn in target:
                    if tid in result[gn]["transcripts"]:
                        found_gn = gn
                        break
                if not found_gn: continue
                tx = result[found_gn]["transcripts"][tid]
                if feat == "exon": tx["exons"].append([start, end])
                elif feat == "CDS": tx["CDS"].append([start, end]); tx["is_coding"] = True
                elif feat == "five_prime_UTR": tx["UTR_5"].append([start, end])
                elif feat == "three_prime_UTR": tx["UTR_3"].append([start, end])

    def _clean(iv):
        if not iv: return []
        u = sorted({(a, b) for a, b in iv})
        m = [list(u[0])]
        for s, e in u[1:]:
            if s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
            else: m.append([s, e])
        return m

    for gn in target:
        for tid, tx in result[gn]["transcripts"].items():
            for field in ["exons", "CDS", "UTR_5", "UTR_3"]:
                tx[field] = _clean(tx[field])
        def _tx_score(tx):
            has_both_utr = (1 if (tx["UTR_5"] and tx["UTR_3"]) else 0)
            n_exons = len(tx["exons"])
            cds_len = sum(e - s for s, e in tx["CDS"]) if tx["CDS"] else 0
            return (not tx["is_coding"], -has_both_utr, -n_exons, -cds_len, tx["id"])
        result[gn]["transcripts"] = sorted(
            result[gn]["transcripts"].values(), key=_tx_score)
    return result


def parse_annotations(filepath, gene_names):
    """Parse gene structures from GTF or GFF3 (auto-detected by extension).
    Supports .gz compressed files."""
    base = os.path.basename(filepath).lower()
    if "gff3" in base or "gff" in base:
        print(f"  Detected GFF3 format: {filepath}")
        return parse_gff3(filepath, gene_names)
    else:
        print(f"  Detected GTF format: {filepath}")
        return parse_gtf(filepath, gene_names)


def _open_reader(filepath):
    """Auto-detect file type and return appropriate reader."""
    base = filepath.lower()
    if base.endswith(".bw") or base.endswith(".bigwig") or base.endswith(".bigwig"):
        from trackpy.bigwig import BigWigReader
        return BigWigReader(filepath)
    elif base.endswith(".bedgraph") or base.endswith(".bedgraph.gz") or base.endswith(".bg.gz"):
        from trackpy.bigwig import BedGraphReader
        return BedGraphReader(filepath)
    else:
        # Try BigWig first, fallback to bedGraph
        from trackpy.bigwig import BigWigReader, BedGraphReader
        try:
            return BigWigReader(filepath)
        except ValueError:
            return BedGraphReader(filepath)


def load_gene_data(genes_dict, bw_paths, flank_up=3000, flank_down=3000):
    """Load bigWig/bedGraph data for all genes. Returns {gene: {region, tracks, ymax}}."""
    import numpy as np

    # Check chromosome prefix consistency
    bw_chroms = set()
    for fp in bw_paths.values():
        with _open_reader(fp) as rdr:
            bw_chroms.update(rdr.chromosomes.keys())
    bw_has_chr = any(c.startswith("chr") for c in bw_chroms)

    for gn, gi in genes_dict.items():
        gtf_chrom = gi.get("chr")
        if not gtf_chrom:
            print(f"  Warning: {gn} not found in annotation, skipping")
            continue
        gtf_has_chr = gtf_chrom.startswith("chr")
        if bw_has_chr and not gtf_has_chr:
            gi["chr"] = f"chr{gtf_chrom}"
        elif not bw_has_chr and gtf_has_chr:
            gi["chr"] = gtf_chrom[3:]

    track_labels = list(bw_paths.keys())
    data = {}
    for gn, gi in genes_dict.items():
        if not gi.get("chr"):
            continue
        rs = gi["start"] - flank_up
        re = gi["end"] + flank_down
        chrom = gi["chr"]
        data[gn] = {"region": (rs, re), "chrom": chrom, "tracks": {}}
        for lbl, fp in bw_paths.items():
            with _open_reader(fp) as rdr:
                data[gn]["tracks"][lbl] = rdr.query(chrom, rs, re)
        # Per-gene unified ymax (all tracks)
        all_vals = [v for lbl in track_labels for _, _, v in data[gn]["tracks"][lbl]]
        data[gn]["ymax"] = float(np.percentile(all_vals, 99)) if all_vals else 1.0
        # Per-track ymax
        data[gn]["track_ymax"] = {}
        for lbl in track_labels:
            vals = [v for _, _, v in data[gn]["tracks"][lbl]]
            data[gn]["track_ymax"][lbl] = float(np.percentile(vals, 99)) if vals else 1.0
    return data
