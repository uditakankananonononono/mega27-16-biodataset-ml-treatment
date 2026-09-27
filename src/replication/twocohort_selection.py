"""Two-cohort cross-stability selection (PREREG_TWOCOHORT.md section 2, locked @b476288
BEFORE this ran). Discovery cohorts ONLY: GSE2034 (relapse) + METABRIC (RFS).
Recipe per cohort (identical to experiments/gene_stability.py): z-score, top-2000
symbols by |t|, 40 bootstrap L1 logistic fits (C=0.5, max_iter=500,
StratifiedKFold(5, shuffle, seed 0) rotating splits, bootstrap resample, rng seed 0),
stability = fraction of fits keeping majority sign. Panel = stability>=0.975 in BOTH
cohorts AND majority signs agree. No test-cohort data is touched here.
Usage: twocohort_selection.py [gse2034|metabric|panel]
"""
import json, sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results" / "replication"
OUT.mkdir(parents=True, exist_ok=True)


def tstats_chunked(X, y, step=2000):
    t = np.zeros(X.shape[1], dtype=np.float32)
    for i in range(0, X.shape[1], step):
        c = np.nan_to_num(X[:, i:i + step].astype(np.float32))
        sd = c.std(0) + 1e-8
        t[i:i + step] = (c[y == 1].mean(0) - c[y == 0].mean(0)) / sd
        del c
    return t


def stability_run(X, y, names):
    t = tstats_chunked(X, y)
    idx = np.argsort(-np.abs(t))[:2000]
    idx.sort()
    Xs = np.nan_to_num(X[:, idx].astype(np.float32))
    Xs = (Xs - Xs.mean(0)) / (Xs.std(0) + 1e-8)
    gsel = [names[i] for i in idx]
    rng = np.random.RandomState(0)
    skf = StratifiedKFold(5, shuffle=True, random_state=0)
    folds = list(skf.split(Xs, y))
    coefs = []
    for b in range(40):
        tr, _ = folds[b % 5]
        sub = rng.choice(tr, size=len(tr), replace=True)
        clf = LogisticRegression(max_iter=500, C=0.5).fit(Xs[sub], y[sub])
        coefs.append(clf.coef_[0])
        if (b + 1) % 10 == 0:
            print(f"  {b+1}/40", flush=True)
    C = np.array(coefs)
    sign = np.sign(C.mean(0))
    stability = (np.sign(C) == sign).mean(0)
    return {g: {"stability": float(s), "sign": int(sg), "mean_coef": float(c)}
            for g, s, sg, c in zip(gsel, stability, sign, C.mean(0))}


def gse2034():
    sys.path.insert(0, str(ROOT / "src"))
    from biomedml.io_geo import parse_series_matrix
    X, probes, meta = parse_series_matrix(str(ROOT / "data" / "GSE2034_matrix.txt.gz"))
    key = [k for k in meta if "relapse" in k.lower() or "metasta" in k.lower() or "event" in k.lower()]
    y = np.array([1 if v.strip().lower() in ("1", "yes", "true") else 0 for v in meta[key[0]]])
    p2s = json.load(open(ROOT / "results" / "u133a_probe_to_symbol.json"))
    X = X.astype(np.float32)
    t = np.abs(tstats_chunked(X, y))
    # max-|t| probe per symbol
    best = {}
    for i, p in enumerate(probes):
        s = p2s.get(p)
        if s and (s not in best or t[i] > best[s][1]):
            best[s] = (i, t[i])
    idxs = sorted(v[0] for v in best.values())
    names = [p2s[probes[i]] for i in idxs]
    print("symbols:", len(names), "samples:", X.shape[0], flush=True)
    res = stability_run(X[:, idxs], y, names)
    json.dump(res, open(OUT / "sel_gse2034.json", "w"))
    print("gse2034 arm done:", len(res), flush=True)


def metabric():
    z = np.load(OUT / "metabric_matrix.npz", allow_pickle=True)
    mat, gene_ids, sample_ids = z["mat"], list(z["gene_ids"]), list(z["sample_ids"])
    raw = json.load(open(OUT / "metabric_raw.json"))
    pat_rfs, samp2pat = raw["pat_rfs"], raw["samp2pat"]
    idxj = json.load(open(OUT / "mb_index_symbols.json")) if (OUT / "mb_index_symbols.json").exists() else None
    # symbols from index.json were saved in chunk dir (deleted); rebuild from pull index copy
    sym_map = json.load(open(OUT / "mb_symbols.json"))
    labels = np.array([1 if pat_rfs.get(samp2pat.get(s, ""), {}).get("RFS_STATUS", "").startswith("1")
                       else 0 if pat_rfs.get(samp2pat.get(s, ""), {}).get("RFS_STATUS", "").startswith("0")
                       else -1 for s in sample_ids])
    keep = labels >= 0
    X = mat[:, keep].T
    y = labels[keep]
    names = [sym_map.get(str(g), str(g)) for g in gene_ids]
    # dedupe symbols: keep max-|t| (chunked, no full-size copies)
    t = np.abs(tstats_chunked(X, y))
    best = {}
    for i, s in enumerate(names):
        if s and (s not in best or t[i] > best[s][1]):
            best[s] = (i, t[i])
    idxs = sorted(v[0] for v in best.values())
    names2 = [names[i] for i in idxs]
    print("symbols:", len(names2), "samples:", X.shape[0], "events:", int(y.sum()), flush=True)
    res = stability_run(X[:, idxs], y, names2)
    json.dump(res, open(OUT / "sel_metabric.json", "w"))
    print("metabric arm done:", len(res), flush=True)


def panel():
    a = json.load(open(OUT / "sel_gse2034.json"))
    b = json.load(open(OUT / "sel_metabric.json"))
    genes = []
    for g in set(a) & set(b):
        if a[g]["stability"] >= 0.975 and b[g]["stability"] >= 0.975 and a[g]["sign"] == b[g]["sign"]:
            genes.append({"gene": g, "sign": a[g]["sign"],
                          "stab_gse2034": round(a[g]["stability"], 3),
                          "stab_metabric": round(b[g]["stability"], 3)})
    genes.sort(key=lambda r: -min(r["stab_gse2034"], r["stab_metabric"]))
    if len(genes) > 30:
        genes = genes[:30]
    json.dump({"rule": "stability>=0.975 both cohorts + sign agreement (PREREG_TWOCOHORT.md)",
               "n_panel": len(genes), "panel": genes},
              open(OUT / "twocohort_panel.json", "w"), indent=1)
    print("PANEL:", len(genes))
    for g in genes:
        print(g)


if __name__ == "__main__":
    {"gse2034": gse2034, "metabric": metabric, "panel": panel}[sys.argv[1]]()
