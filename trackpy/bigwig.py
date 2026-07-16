"""Pure Python bigWig file reader. No C dependencies required."""

import struct
import zlib


class BigWigReader:
    """Read and query bigWig files."""

    BIGWIG_MAGIC = 0x888FFC26

    def __init__(self, filepath):
        self._f = open(filepath, "rb")
        self._chrom_dict = {}
        self.chrom_to_id = {}
        self._read_header()
        self._read_chrom_tree()

    def _read_header(self):
        f = self._f
        f.seek(0)
        self.magic = struct.unpack("<I", f.read(4))[0]
        if self.magic != self.BIGWIG_MAGIC:
            raise ValueError(f"Not a valid bigWig file (magic={hex(self.magic)})")
        f.read(2)  # version
        self.zoom_levels = struct.unpack("<H", f.read(2))[0]
        self.chrom_tree_offset = struct.unpack("<Q", f.read(8))[0]
        self.unzoomed_data_count = struct.unpack("<Q", f.read(8))[0]
        self.unzoomed_index_offset = struct.unpack("<Q", f.read(8))[0]
        f.read(28)  # field_count..reserved

        self.zoom_headers = []
        for _ in range(self.zoom_levels):
            self.zoom_headers.append({
                "reduction_level": struct.unpack("<I", f.read(4))[0],
                "data_offset": struct.unpack("<Q", f.read(8))[0],
                "index_offset": struct.unpack("<Q", f.read(8))[0],
            })
            f.read(4)

    def _read_chrom_tree(self):
        f = self._f
        f.seek(self.chrom_tree_offset)
        struct.unpack("<I", f.read(4))  # magic
        block_size = struct.unpack("<I", f.read(4))[0]
        key_size   = struct.unpack("<I", f.read(4))[0]
        val_size   = struct.unpack("<I", f.read(4))[0]
        f.read(16)  # item_count + reserved
        self._visited = set()
        self._read_bptree(f.tell(), block_size, key_size, val_size)
        del self._visited
        self.chrom_to_id = {n: cid for n, (cid, _) in self._chrom_dict.items()}

    def _read_bptree(self, offset, bs, ks, vs):
        if offset in self._visited:
            return
        self._visited.add(offset)
        f = self._f
        f.seek(offset)
        is_leaf = struct.unpack("B", f.read(1))[0]
        f.read(1)
        count = struct.unpack("<H", f.read(2))[0]
        if count > 10000:
            return
        for _ in range(count):
            key = f.read(ks)
            if is_leaf:
                val = f.read(vs)
                cid = struct.unpack("<I", val[:4])[0]
                csz = struct.unpack("<I", val[4:8])[0]
                self._chrom_dict[key.rstrip(b"\x00").decode()] = (cid, csz)
            else:
                child = struct.unpack("<Q", f.read(8))[0]
                if child > 0:
                    self._read_bptree(child, bs, ks, vs)
        if not is_leaf:
            last = struct.unpack("<Q", f.read(8))[0]
            if last > 0:
                self._read_bptree(last, bs, ks, vs)

    def _read_rtree(self, offset):
        f = self._f
        f.seek(offset)
        f.read(48)  # header
        nodes, stack, seen = [], [f.tell()], set()
        while stack:
            off = stack.pop()
            if off in seen:
                continue
            seen.add(off)
            f.seek(off)
            is_leaf = struct.unpack("B", f.read(1))[0]
            f.read(1)
            count = struct.unpack("<H", f.read(2))[0]
            if count > 10000:
                continue
            for _ in range(count):
                sci = struct.unpack("<I", f.read(4))[0]
                sb  = struct.unpack("<I", f.read(4))[0]
                eci = struct.unpack("<I", f.read(4))[0]
                eb  = struct.unpack("<I", f.read(4))[0]
                if is_leaf:
                    doff = struct.unpack("<Q", f.read(8))[0]
                    dsz  = struct.unpack("<Q", f.read(8))[0]
                    nodes.append({"sci": sci, "sb": sb, "eci": eci, "eb": eb,
                                   "doff": doff, "dsz": dsz})
                else:
                    child = struct.unpack("<Q", f.read(8))[0]
                    stack.append(child)
        return nodes

    def resolve_chrom(self, chrom):
        if chrom in self.chrom_to_id:
            return self.chrom_to_id[chrom]
        if chrom.startswith("chr") and chrom[3:] in self.chrom_to_id:
            return self.chrom_to_id[chrom[3:]]
        if not chrom.startswith("chr") and f"chr{chrom}" in self.chrom_to_id:
            return self.chrom_to_id[f"chr{chrom}"]
        raise KeyError(f"Chromosome '{chrom}' not found in {list(self.chrom_to_id.keys())}")

    def query(self, chrom, start, end):
        """Return list of (start, end, value) for a genomic region."""
        cid = self.resolve_chrom(chrom)
        results = []
        for node in self._read_rtree(self.unzoomed_index_offset):
            if node["sci"] != cid and node["eci"] != cid:
                continue
            if node["eb"] < start or node["sb"] > end:
                continue
            self._f.seek(node["doff"])
            try:
                raw = zlib.decompress(self._f.read(node["dsz"]))
            except zlib.error:
                continue
            off = 0
            if len(raw) < 24:
                continue
            sec_start = struct.unpack_from("<I", raw, 4)[0]
            sec_end   = struct.unpack_from("<I", raw, 8)[0]
            sec_step  = struct.unpack_from("<I", raw, 12)[0]
            sec_span  = struct.unpack_from("<I", raw, 16)[0]
            sec_type  = struct.unpack_from("B", raw, 20)[0]
            sec_count = struct.unpack_from("<H", raw, 22)[0]
            off = 24

            if sec_type == 1:  # bedGraph
                for _ in range(sec_count):
                    if off + 12 > len(raw): break
                    i_s = struct.unpack_from("<I", raw, off)[0]; off += 4
                    i_e = struct.unpack_from("<I", raw, off)[0]; off += 4
                    val = struct.unpack_from("<f", raw, off)[0]; off += 4
                    if i_e > start and i_s < end:
                        results.append((i_s, i_e, val))
            elif sec_type in (2, 3):  # variableStep / fixedStep
                for i in range(sec_count):
                    if off + 4 > len(raw): break
                    val = struct.unpack_from("<f", raw, off)[0]; off += 4
                    i_s = sec_start + i * sec_step
                    i_e = i_s + sec_span
                    if i_e > start and i_s < end:
                        results.append((i_s, i_e, val))

        results.sort(key=lambda x: x[0])
        return results

    @property
    def chromosomes(self):
        return {n: sz for n, (_, sz) in self._chrom_dict.items()}

    def close(self):
        if self._f:
            self._f.close()
            self._f = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


import gzip as _gzip

class BedGraphReader:
    """Read bedGraph files (plain or .gz compressed).
    Format: chrom  start  end  value (tab-separated)."""

    def __init__(self, filepath):
        self.filepath = filepath
        self._is_gz = filepath.endswith(".gz")
        self._f = _gzip.open(filepath, "rt") if self._is_gz else open(filepath, "r")

    def query(self, chrom, start, end):
        """Return list of (start, end, value) in region."""
        results = []
        # Try with/without chr prefix
        alt_chrom = None
        if chrom.startswith("chr"):
            alt_chrom = chrom[3:]
        else:
            alt_chrom = f"chr{chrom}"

        self._f.seek(0)
        for line in self._f:
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            c = parts[0]
            if c != chrom and c != alt_chrom:
                continue
            s = int(parts[1])
            e = int(parts[2])
            if e <= start or s >= end:
                continue
            try:
                v = float(parts[3])
            except ValueError:
                continue
            results.append((s, e, v))
        results.sort(key=lambda x: x[0])
        return results

    @property
    def chromosomes(self):
        """Return {chrom: size} by scanning the file."""
        chroms = {}
        self._f.seek(0)
        for line in self._f:
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            c = parts[0]
            e = int(parts[2])
            if c not in chroms or e > chroms[c]:
                chroms[c] = e
        return chroms

    def close(self):
        if self._f:
            self._f.close()
            self._f = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
