"""SCAN-B (GSE96058) locked scoring - SECONDARY arm, PREREG_REPLICATION.md section 2/3.
GEO clinical annotation carries OS only -> declared OS-endpoint arm (prereg fallback),
NOT pooled with RFS arms. Same locked rule: score = sum sign_g * z(expr), no re-estimation.
"""
import csv, gzip, json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
DAT = ROOT / "data" / "scanb"
OUT = ROOT / "results" / "replication"
OUT.mkdir(parents=True, exist_ok=True)

# (matrix row name, discovery sign); EEF2KMT=FAM86A, BICRA=GLTSCR1 (REPLICATION_ACCESS.md)
PANEL = [("EEF2KMT", "FAM86A", -1), ("WFDC1", "WFDC1", 1), ("PSMB1", "PSMB1", -1),
         ("CALML3", "CALML3", -1), ("NPPC", "NPPC", -1), ("PASK", "PASK", 1),
         ("GLTSCR1", "GLTSCR1", 1), ("BLZF1", "BLZF1", 1), ("CMC4", "CMC4", -1),
         ("COX11", "COX11", 1), ("CDH13", "CDH13", 1), ("FUT9", "FUT9", -1),
         ("ALDH3A2", "ALDH3A2", -1), ("CWF19L1", "CWF19L1", -1), ("POF1B", "POF1B", 1)]


def parse_matrix(path):
    clin = {}
    with gzip.open(path, "rt") as f:
        titles = None
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
                for t in titles:
                    clin[t] = {}
            elif line.startswith("!Sample_characteristics_ch1") and titles:
                vals = line.rstrip("\n").split("\t")[1:]
                for t, v in zip(titles, vals):
                    v = v.strip('"')
                    if ": " in v:
                        k, val = v.split(": ", 1)
                        clin[t][k.lower()] = val
    return clin


def auc(scores, labels):
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores))
    ranks[order] = np.arange(1, len(scores) + 1)
    s_sorted = scores[order]
    i = 0
    while i < len(s_sorted):
        j = i
        while j + 1 < len(s_sorted) and s_sorted[j + 1] == s_sorted[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2
        i = j + 1
    pos = labels == 1
    return (ranks[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * (~pos).sum())


def main():
    clin = {}
    for f in ["GSE96058-GPL11154_series_matrix.txt.gz", "GSE96058-GPL18573_series_matrix.txt.gz"]:
        clin.update(parse_matrix(DAT / f))
    with open(DAT / "panel_rows.csv") as f:
        rows = list(csv.reader(f))
    cols = rows[0][1:]
    expr = {r[0].strip('"'): np.array([float(x) for x in r[1:]]) for r in rows[1:]}
    keep_cols, labels = [], []
    for i, c in enumerate(cols):
        rec = clin.get(c)
        if rec and rec.get("overall survival event", "NA") not in ("NA", ""):
            keep_cols.append(i)
            labels.append(int(float(rec["overall survival event"])))
    labels = np.array(labels)
    X = np.zeros((len(keep_cols), len(PANEL)))
    for j, (gene, rowname, sign) in enumerate(PANEL):
        v = expr[rowname][keep_cols]
        mu, sd = v.mean(), v.std()
        X[:, j] = (v - mu) / (sd if sd > 0 else 1.0) * sign
    score = X.sum(axis=1)
    a = auc(score, labels)
    rng = np.random.default_rng(2026)
    n = len(score)
    boots = []
    for _ in range(2000):
        idx = rng.integers(0, n, n)
        if labels[idx].sum() in (0, n):
            continue
        boots.append(auc(score[idx], labels[idx]))
    lo = float(np.percentile(boots, 5))
    conc = []
    for j, (gene, rowname, sign) in enumerate(PANEL):
        r = float(np.corrcoef(X[:, j] * sign, labels)[0, 1])
        conc.append({"gene": gene, "discovery_sign": sign, "scanb_corr": r,
                     "concordant": bool(np.sign(r) == sign)})
    n_match = sum(c["concordant"] for c in conc)
    res = {
        "cohort": "SCAN-B GSE96058 (GEO public, pulled 2026-09-27; RNA-seq, both platforms)",
        "endpoint": "OS (overall survival event) - DECLARED OS-endpoint arm; GEO annotation carries no relapse field (prereg section 2 fallback)",
        "n_samples": int(n), "n_events": int(labels.sum()),
        "genes_used": len(PANEL),
        "auc": float(a), "boot_lo_one_sided_95": lo,
        "sign_concordance": f"{n_match}/{len(PANEL)}",
        "per_gene": conc,
        "secondary_test_pass": bool(lo > 0.5 and n_match / len(PANEL) > 0.5),
        "gate": "PREREG_REPLICATION.md section 3 (SECONDARY 1)",
    }
    json.dump(res, open(OUT / "scanb_secondary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "per_gene"}, indent=1))


if __name__ == "__main__":
    main()
