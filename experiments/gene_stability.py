"""Gene-stability arm: which genes robustly carry the metastasis signal?

Bootstrap the 5-fold-CV logistic pipeline 40 times; a gene's stability is
the fraction of fits where its coefficient keeps the majority sign, times
mean |coef|. Outputs a ranked, falsifiable stable-panel candidate list.
"""
import json, os, sys
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
sys.path.insert(0, "src")
from biomedml.io_geo import parse_series_matrix

def main():
    X, genes, meta = parse_series_matrix("data/GSE2034_matrix.txt.gz")
    key = [k for k in meta if "relapse" in k.lower() or "metasta" in k.lower() or "event" in k.lower()]
    y = np.array([1 if v.strip().lower() in ("1", "yes", "true") else 0 for v in meta[key[0]]])
    X = np.nan_to_num(X)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    # prefilter to top 2000 by |t| to bound compute (same as pipeline)
    from biomedml.stats import de_genes
    t, _ = de_genes(X, y)
    idx = np.argsort(-np.abs(t))[:2000]
    Xs, gsel = X[:, idx], [genes[i] for i in idx]

    rng = np.random.RandomState(0)
    B, coefs = 40, []
    skf = StratifiedKFold(5, shuffle=True, random_state=0)
    folds = list(skf.split(Xs, y))
    for b in range(B):
        tr, _ = folds[b % 5]
        sub = rng.choice(tr, size=len(tr), replace=True)
        clf = LogisticRegression(max_iter=500, C=0.5).fit(Xs[sub], y[sub])
        coefs.append(clf.coef_[0])
    C = np.array(coefs)                      # B x G
    sign = np.sign(C.mean(0))
    stability = (np.sign(C) == sign).mean(0)
    score = stability * np.abs(C.mean(0))
    order = np.argsort(-score)[:30]
    panel = [dict(probe=gsel[i], stability=round(float(stability[i]), 3),
                  mean_coef=round(float(C.mean(0)[i]), 4),
                  direction="risk" if sign[i] > 0 else "protective") for i in order]
    os.makedirs("results", exist_ok=True)
    json.dump(dict(n_boot=40, n_genes_prefilter=2000, panel=panel),
              open("results/gene_stability.json", "w"), indent=1)
    for p in panel[:15]:
        print(p)

if __name__ == "__main__":
    main()
