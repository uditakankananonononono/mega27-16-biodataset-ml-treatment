"""Pivot after the 70/30-split negative: 5-fold CV, top-100 DE panel per fold
(selection inside fold), logistic + CNN, report honest CV accuracy/AUC."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from biomedml.io_geo import parse_series_matrix
from biomedml.stats import de_genes
from biomedml.models import LogRegNP, ExprCNN, train_torch, predict_torch, auc_score

X, genes, meta = parse_series_matrix(os.path.join(
    os.path.dirname(__file__), "..", "data", "GSE2034_matrix.txt.gz"))
raw = meta["bone relapses (1=yes, 0=no)"]
keep = [i for i, v in enumerate(raw) if v.strip() in ("0", "1")]
y = np.array([int(raw[i]) for i in keep])
X = np.nan_to_num(np.log2(np.maximum(X[keep], 1.0)))
X = (X - X.mean(0)) / (X.std(0) + 1e-6)

rng = np.random.default_rng(9)
perm = rng.permutation(len(y))
folds = np.array_split(perm, 5)
accs_lr, aucs_lr, accs_cnn = [], [], []
base = []
for k in range(5):
    te = folds[k]
    tr = np.concatenate([folds[j] for j in range(5) if j != k])
    _, p_tr = de_genes(X[tr], y[tr])
    sel = np.argsort(p_tr)[:100]
    lr = LogRegNP(100, lr=0.2, epochs=800, l2=1e-2)
    lr.fit(X[tr][:, sel], y[tr])
    pred = lr.predict(X[te][:, sel])
    accs_lr.append(float((pred == y[te]).mean()))
    scores = X[te][:, sel] @ (lr.W[:, 1] - lr.W[:, 0])
    aucs_lr.append(auc_score(y[te], scores))
    cnn = ExprCNN(100)
    train_torch(cnn, X[tr][:, sel], y[tr], np.arange(len(tr)), epochs=30, lr=1e-3, seed=k)
    accs_cnn.append(float((predict_torch(cnn, X[te][:, sel]) == y[te]).mean()))
    base.append(float(max((y[te] == 0).mean(), (y[te] == 1).mean())))
    print("fold", k, "lr", round(accs_lr[-1], 3), round(aucs_lr[-1], 3),
          "cnn", round(accs_cnn[-1], 3), flush=True)
res = dict(cv=5, genes_per_fold=100,
           logreg_acc=[float(np.mean(accs_lr)), float(np.std(accs_lr))],
           logreg_auc=[float(np.mean(aucs_lr)), float(np.std(aucs_lr))],
           cnn_acc=[float(np.mean(accs_cnn)), float(np.std(accs_cnn))],
           majority=[float(np.mean(base))],
           verdict=("signal-above-chance" if np.mean(aucs_lr) > 0.6
                    else "honest-negative: expression panel does not beat chance robustly"))
OUT = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "cv_results.json"), "w") as f:
    json.dump(res, f, indent=1)
print("CV RESULT", res, flush=True)
