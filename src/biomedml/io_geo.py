"""Parse a GEO series-matrix file (GSEnnn_series_matrix.txt.gz)."""
from __future__ import annotations
import gzip
import numpy as np


def parse_series_matrix(path):
    """Returns (X [samples x genes], gene_ids, meta dict of column->values)."""
    opener = gzip.open if path.endswith(".gz") else open
    meta = {}
    header = None
    rows = []
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_characteristics"):
                parts = line.rstrip("\n").split("\t")
                key = parts[1].split(":")[0].strip('" ')
                vals = [p.split(":", 1)[-1].strip('" ') for p in parts[1:]]
                meta.setdefault(key, vals)
            elif line.startswith("!Sample_geo_accession"):
                header = [p.strip('"') for p in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                cols = next(f).rstrip("\n").split("\t")
                genes = []
                data = []
                for l2 in f:
                    if l2.startswith("!series_matrix_table_end"):
                        break
                    p = l2.rstrip("\n").split("\t")
                    genes.append(p[0].strip('"'))
                    data.append([float(x) if x not in ("", "null", "NA") else float("nan")
                                 for x in p[1:]])
                X = np.array(data, dtype=np.float32).T  # samples x genes
                return X, genes, meta
    raise ValueError("no series matrix table found")
