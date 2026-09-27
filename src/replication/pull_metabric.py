"""METABRIC replication pull + locked scoring (PREREG_REPLICATION.md, locked 2026-09-27 @e54d8bc).
Single scoring rule: score_i = sum_g sign(coef_discovery_g) * z(expr_ig)  (within-cohort z).
Primary test: directional AUC vs RFS_STATUS, one-sided 95% bootstrap CI (patient-level),
plus per-gene sign concordance vs discovery directions. No re-estimation on METABRIC data.
"""
import json, sys, time, urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results" / "replication"
OUT.mkdir(parents=True, exist_ok=True)
API = "https://www.cbioportal.org/api"

# Frozen panel (PREREG_REPLICATION.md section 1). GLTSCR1 is alias of BICRA (REPLICATION_ACCESS.md).
PANEL = [
    ("222013_x_at", "EEF2KMT", -1), ("219478_at", "WFDC1", 1), ("200876_s_at", "PSMB1", -1),
    ("210019_at", "CALML3", -1), ("221348_at", "NPPC", -1), ("206594_at", "PASK", 1),
    ("219445_at", "GLTSCR1", 1), ("32088_at", "BLZF1", 1), ("216862_s_at", "CMC4", -1),
    ("214277_at", "COX11", 1), ("204726_at", "CDH13", 1), ("214046_at", "FUT9", -1),
    ("202053_s_at", "ALDH3A2", -1), ("218787_x_at", "CWF19L1", -1), ("219756_s_at", "POF1B", 1),
]
ALIAS = {"GLTSCR1": "BICRA"}


def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Accept": "application/json",
                                                          "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def main():
    # 1. Entrez IDs
    entrez = {}
    for probe, gene, sign in PANEL:
        sym = ALIAS.get(gene, gene)
        if gene not in entrez:
            hits = [g for g in get(f"{API}/genes?keyword={sym}") if g["hugoGeneSymbol"] == sym]
            assert hits, f"gene not found: {sym}"
            entrez[gene] = hits[0]["entrezGeneId"]
        time.sleep(0.2)
    print("entrez:", entrez, flush=True)

    # 2. Sample list for the mrna profile
    sample_ids = get(f"{API}/sample-lists/brca_metabric_mrna/sample-ids")
    print("samples:", len(sample_ids), flush=True)

    # 3. mRNA values for panel genes (chunked sample fetch)
    expr = {g: {} for g in entrez}
    ids = list(entrez.values())
    B = 400
    for i in range(0, len(sample_ids), B):
        chunk = sample_ids[i:i + B]
        body = json.dumps({"entrezGeneIds": ids, "sampleIds": chunk}).encode()
        rows = get(f"{API}/molecular-profiles/brca_metabric_mrna/molecular-data/fetch", body)
        for r in rows:
            g = next(k for k, v in entrez.items() if v == r["entrezGeneId"])
            expr[g][r["sampleId"]] = r["value"]
        print(f"fetched {i + len(chunk)}/{len(sample_ids)}", flush=True)
        time.sleep(0.3)

    # 4. Clinical RFS (patient-level, chunked; sample ids == patient ids in this study)
    clin = []
    for i in range(0, len(sample_ids), 400):
        body = json.dumps({"attributeIds": ["RFS_STATUS", "RFS_MONTHS"],
                           "ids": sample_ids[i:i + 400]}).encode()
        clin += get(f"{API}/studies/brca_metabric/clinical-data/fetch?clinicalDataType=PATIENT&projection=SUMMARY", body)
        time.sleep(0.3)
    pat_rfs = {}
    for c in clin:
        pat_rfs.setdefault(c["patientId"], {})[c["clinicalAttributeId"]] = c["value"]
    samp2pat = {s["sampleId"]: s["patientId"]
                for s in get(f"{API}/studies/brca_metabric/samples?projection=SUMMARY&pageSize=100000")}
    np.save(OUT / "metabric_sample_ids.npy", np.array(sample_ids))
    json.dump({"entrez": entrez, "expr": expr, "pat_rfs": pat_rfs, "samp2pat": samp2pat},
              open(OUT / "metabric_raw.json", "w"))
    print("pulled:", {g: len(v) for g, v in expr.items()}, "patients with RFS:", len(pat_rfs), flush=True)


if __name__ == "__main__":
    main()
