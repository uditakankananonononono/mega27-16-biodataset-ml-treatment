"""Locked scoring of the frozen 15-gene panel on METABRIC (PREREG_REPLICATION.md @e54d8bc).
score_i = sum_g sign_g * z(expr_ig); primary test: directional AUC vs RFS_STATUS,
one-sided 95% bootstrap CI (patient-level, 2000 resamples); sign concordance per gene.
No re-estimation on METABRIC data. Run AFTER pull_metabric.py.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results" / "replication"

PANEL = [  # (probe, gene, discovery sign) - frozen, PREREG section 1
    ("222013_x_at", "EEF2KMT", -1), ("219478_at", "WFDC1", 1), ("200876_s_at", "PSMB1", -1),
    ("210019_at", "CALML3", -1), ("221348_at", "NPPC", -1), ("206594_at", "PASK", 1),
    ("219445_at", "GLTSCR1", 1), ("32088_at", "BLZF1", 1), ("216862_s_at", "CMC4", -1),
    ("214277_at", "COX11", 1), ("204726_at", "CDH13", 1), ("214046_at", "FUT9", -1),
    ("202053_s_at", "ALDH3A2", -1), ("218787_x_at", "CWF19L1", -1), ("219756_s_at", "POF1B", 1),
]


def auc(scores, labels):
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores))
    ranks[order] = np.arange(1, len(scores) + 1)
    # tie-average ranks
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
    n_pos, n_neg = pos.sum(), (~pos).sum()
    return (ranks[pos].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def main():
    raw = json.load(open(OUT / "metabric_raw.json"))
    expr, pat_rfs, samp2pat = raw["expr"], raw["pat_rfs"], raw["samp2pat"]
    sample_ids = list(np.load(OUT / "metabric_sample_ids.npy"))
    samples = [s for s in sample_ids if s in samp2pat]
    labels = np.array([1 if pat_rfs.get(samp2pat[s], {}).get("RFS_STATUS", "").startswith("1")
                       else 0 if pat_rfs.get(samp2pat[s], {}).get("RFS_STATUS", "").startswith("0")
                       else -1 for s in samples])
    keep = labels >= 0
    samples = [s for s, k in zip(samples, keep) if k]
    labels = labels[keep]
    X = np.zeros((len(samples), len(PANEL)))
    dropped = []
    for j, (probe, gene, sign) in enumerate(PANEL):
        v = np.array([expr[gene].get(s, np.nan) for s in samples], dtype=float)
        if np.isnan(v).all():
            dropped.append(gene)
            continue
        mu, sd = np.nanmean(v), np.nanstd(v)
        X[:, j] = np.where(np.isnan(v), 0.0, (v - mu) / (sd if sd > 0 else 1.0)) * sign
    score = X.sum(axis=1)
    a = auc(score, labels)
    rng = np.random.default_rng(2026)
    boots = np.array([auc(score[rng.integers(0, len(score), len(score))],
                          labels[rng.integers(0, len(labels), len(labels))]) for _ in range(2000)])
    # patient-level bootstrap: resample indices iid (patients == samples here, 1 tumour/patient)
    lo = np.percentile(boots, 5)  # one-sided 95%
    conc = []
    for j, (probe, gene, sign) in enumerate(PANEL):
        if gene in dropped:
            continue
        r = np.corrcoef(X[:, j] * sign, labels)[0, 1]  # raw-gene direction vs label
        conc.append((gene, sign, float(r), bool(np.sign(r) == sign)))
    n_match = sum(c[3] for c in conc)
    res = {
        "cohort": "METABRIC brca_metabric (cBioPortal public API, pulled 2026-09-27)",
        "n_samples": len(samples), "n_relapsed": int(labels.sum()),
        "genes_used": len(PANEL) - len(dropped), "genes_dropped_from_profile": dropped,
        "primary_auc": float(a), "boot_lo_one_sided_95": float(lo),
        "sign_concordance": f"{n_match}/{len(conc)}",
        "per_gene": [{"gene": g, "discovery_sign": s, "metabric_corr": r, "concordant": m}
                     for g, s, r, m in conc],
        "primary_test_pass": bool(lo > 0.5 and n_match / len(conc) > 0.5),
        "gate": "PREREG_REPLICATION.md section 3-4",
    }
    json.dump(res, open(OUT / "metabric_primary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "per_gene"}, indent=1))
    for c in conc:
        print(c)


if __name__ == "__main__":
    main()
