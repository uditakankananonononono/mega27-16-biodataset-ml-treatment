"""Pull full METABRIC HT-12 expression matrix (all profile genes x mrna samples)
for the two-cohort selection arm (PREREG_TWOCOHORT.md section 2, discovery use only).
Gene-chunked, disk-checkpointed, resume-safe. Saves results/replication/metabric_matrix.npz
Usage: pull_metabric_matrix.py [n_chunks_per_run]
"""
import json, sys, time, urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results" / "replication"
CHUNK_DIR = OUT / "mb_matrix_chunks"
CHUNK_DIR.mkdir(parents=True, exist_ok=True)
API = "https://www.cbioportal.org/api"
GENES_PER_CHUNK = 200


def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Accept": "application/json",
                                                          "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)


def main(n_runs=99):
    idx_p = CHUNK_DIR / "index.json"
    if idx_p.exists():
        idx = json.load(open(idx_p))
        sample_ids, gene_ids, sym = idx["sample_ids"], idx["gene_ids"], idx["symbols"]
    else:
        sample_ids = get(f"{API}/sample-lists/brca_metabric_mrna/sample-ids")
        genes = get(f"{API}/genes?projection=SUMMARY&pageSize=100000")
        gene_ids = sorted(g["entrezGeneId"] for g in genes)
        sym = {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in genes}
        json.dump({"sample_ids": sample_ids, "gene_ids": gene_ids, "symbols": sym},
                  open(idx_p, "w"))
    print(f"samples={len(sample_ids)} genes={len(gene_ids)}", flush=True)
    done = 0
    for ci in range(0, len(gene_ids), GENES_PER_CHUNK):
        p = CHUNK_DIR / f"chunk_{ci}.npy"
        if p.exists():
            done += 1
            continue
        if done >= 0 and n_runs <= 0:
            break
        ids = gene_ids[ci:ci + GENES_PER_CHUNK]
        mat = np.full((len(ids), len(sample_ids)), np.nan, dtype=np.float32)
        col = {s: j for j, s in enumerate(sample_ids)}
        print(f"fetching chunk {ci}...", flush=True)
        body = json.dumps({"entrezGeneIds": ids, "sampleIds": sample_ids}).encode()
        rows = get(f"{API}/molecular-profiles/brca_metabric_mrna/molecular-data/fetch", body)
        rmap = {g: i for i, g in enumerate(ids)}
        for r in rows:
            mat[rmap[r["entrezGeneId"]], col[r["sampleId"]]] = r["value"]
        np.save(p, mat)
        n_runs -= 1
        print(f"chunk {ci} done ({mat.shape})", flush=True)
        time.sleep(0.5)
    print("pass complete", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 99)
