"""Relapse prediction from tumor expression (GEO GSE2034, 286 samples):
DE biomarker panel + 1D-CNN vs logistic baseline, patient-level split."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from biomedml.io_geo import parse_series_matrix
from biomedml.stats import de_genes, bh_fdr, coexpression_graph
from biomedml.models import (ExprCNN, LogRegNP, train_torch, predict_torch,
                             auc_score)

OUT = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)
X, genes, meta = parse_series_matrix(os.path.join(
    os.path.dirname(__file__), "..", "data", "GSE2034_matrix.txt.gz"))
print("X", X.shape, "meta keys", list(meta)[:8], flush=True)
# GSE2034: characteristics include relapse status ("relapse (1=yes, 0=no)")
lab_key = None
for k in meta:
    vals = set(v.lower() for v in meta[k])
    if vals & {"yes", "no"} or "relapse" in k.lower():
        lab_key = k; break
print("label key:", lab_key, flush=True)
raw = meta[lab_key]
keep = [i for i, v in enumerate(raw) if v.strip().lower() in ("yes", "no", "1", "0")]
y = np.array([1 if raw[i].strip().lower() in ("yes", "1") else 0 for i in keep])
X = X[keep]
X = np.nan_to_num(np.log2(np.maximum(X, 1.0)))
X = (X - X.mean(0)) / (X.std(0) + 1e-6)
t, p = de_genes(X, y)
de_mask = bh_fdr(p, 0.05)
top_idx = np.argsort(p)[:500]  # feature selection INSIDE train below
rng = np.random.default_rng(5)
perm = rng.permutation(len(y))
n_tr = int(0.7 * len(y))
tr, te = perm[:n_tr], perm[n_tr:]
# feature selection on train only (no leakage)
t_tr, p_tr = de_genes(X[tr], y[tr])
sel = np.argsort(p_tr)[:400]
Xs_tr, Xs_te = X[tr][:, sel], X[te][:, sel]

cnn = ExprCNN(len(sel))
train_torch(cnn, Xs_tr, y[tr], np.arange(len(tr)), epochs=40, lr=1e-3, batch=32, seed=2)
pred_cnn = predict_torch(cnn, Xs_te)
acc_cnn = float((pred_cnn == y[te]).mean())

lr = LogRegNP(len(sel), lr=0.05, epochs=600, l2=1e-2)
lr.fit(Xs_tr, y[tr])
pred_lr = lr.predict(Xs_te)
acc_lr = float((pred_lr == y[te]).mean())
base = float(max((y[te] == 0).mean(), (y[te] == 1).mean()))
# AUC from logreg scores
scores = Xs_te @ (lr.W[:, 1] - lr.W[:, 0])
auc = auc_score(y[te], scores)

res = dict(dataset="GSE2034", n_samples=int(len(y)), n_genes=int(X.shape[1]),
           n_de_fdr05=int(de_mask.sum()),
           label_key=lab_key, prevalence=float(y.mean()),
           cnn_acc=acc_cnn, logreg_acc=acc_lr, logreg_auc=auc,
           majority_baseline=base,
           beats_baseline=bool(max(acc_cnn, acc_lr) > base + 0.03),
           top_de_genes=[genes[i] for i in np.argsort(p)[:25]])
with open(os.path.join(OUT, "results.json"), "w") as f:
    json.dump(res, f, indent=1)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
ax[0].bar(["1D-CNN", "LogReg", "Majority"], [acc_cnn, acc_lr, base],
          color=["#1f77b4", "#ff7f0e", "#7f7f7f"])
ax[0].set_ylabel("held-out accuracy"); ax[0].set_title("GSE2034 relapse prediction")
ax[1].scatter(-np.log10(p), t, s=2, alpha=0.4)
sig = de_mask
ax[1].scatter(-np.log10(p[sig]), t[sig], s=3, color="r", alpha=0.6)
ax[1].set_xlabel("-log10 p"); ax[1].set_ylabel("t statistic")
ax[1].set_title(f"Volcano: {sig.sum()} FDR<0.05 genes")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "figures", "gse2034.png"), dpi=150)
print("RESULT", {k: v for k, v in res.items() if k != "top_de_genes"}, flush=True)
