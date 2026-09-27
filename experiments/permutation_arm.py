"""Permutation-significance arm for the GSE2034 stability selection.

Null distribution: rerun the FULL locked selection recipe under label
permutations - prefilter (top-2000 |t|) recomputed on the permuted labels,
then 40 bootstrap L2 LR fits C=0.5, stability = sign-consistency fraction. For each permutation we
record (a) how many genes reach stability >= 0.975 (the locked panel
threshold) and (b) the max stability observed. The observed run had 268 genes
at stability >= 0.975 (sel_gse2034.json); p = (1 + #{perm count >= 268}) /
(1 + P). Checkpointed per permutation; resume-safe.
Usage: permutation_arm.py [n_perms]
"""
import json, os, sys, time
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
sys.path.insert(0, "src")
from biomedml.io_geo import parse_series_matrix
from biomedml.stats import de_genes

OUT = "results/permutation_arm.json"
CKPT = "results/permutation_arm_ckpt.json"
STAB_THRESH = 0.975
OBSERVED_N_STABLE = 268  # sel_gse2034.json, locked selection

def main():
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    X, genes, meta = parse_series_matrix("data/GSE2034_matrix.txt.gz")
    key = [k for k in meta if "relapse" in k.lower() or "metasta" in k.lower() or "event" in k.lower()]
    y = np.array([1 if v.strip().lower() in ("1", "yes", "true") else 0 for v in meta[key[0]]])
    X = np.nan_to_num(X)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    skf = StratifiedKFold(5, shuffle=True, random_state=0)
    folds = list(skf.split(X, y))
    state = {"done": [], "null_counts": [], "null_maxstab": []}
    if os.path.exists(CKPT):
        state = json.load(open(CKPT))
    for p in range(P):
        if p in state["done"]:
            continue
        rng = np.random.RandomState(1000 + p)
        yp = rng.permutation(y)
        tp, _ = de_genes(X, yp)  # full pipeline under permutation, incl. prefilter
        Xp = X[:, np.argsort(-np.abs(tp))[:2000]]
        brng = np.random.RandomState(p)
        coefs = []
        for b in range(40):
            tr, _ = folds[b % 5]
            sub = brng.choice(tr, size=len(tr), replace=True)
            clf = LogisticRegression(max_iter=500, C=0.5).fit(Xp[sub], yp[sub])
            coefs.append(clf.coef_[0])
        C = np.array(coefs)
        sign = np.sign(C.mean(0))
        stab = (np.sign(C) == sign).mean(0)
        state["done"].append(p)
        state["null_counts"].append(int((stab >= STAB_THRESH).sum()))
        state["null_maxstab"].append(float(stab.max()))
        json.dump(state, open(CKPT, "w"))
        print(f"perm {p}: n_stable={state['null_counts'][-1]} maxstab={state['null_maxstab'][-1]:.3f} ({len(state['done'])}/{P})", flush=True)
    nc = np.array(state["null_counts"]); ms = np.array(state["null_maxstab"])
    res = {"n_perms": len(nc), "threshold": STAB_THRESH,
           "observed_n_stable": OBSERVED_N_STABLE,
           "null_n_stable_median": float(np.median(nc)),
           "null_n_stable_max": int(nc.max()),
           "p_value_panel_size": (1 + int((nc >= OBSERVED_N_STABLE).sum())) / (1 + len(nc)),
           "null_maxstab_median": float(np.median(ms)),
           "null_counts": nc.tolist(), "null_maxstab": ms.tolist()}
    json.dump(res, open(OUT, "w"), indent=1)
    print("ARM DONE", json.dumps({k: v for k, v in res.items() if not isinstance(v, list)}), flush=True)

if __name__ == "__main__":
    main()
