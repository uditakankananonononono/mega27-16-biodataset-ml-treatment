"""Graph arm: smooth expression over the train-fold co-expression graph
(GCN-style propagation over genes), then classify. Compare CV AUC vs the
plain top-100 panel from pivot_cv.py."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from biomedml.io_geo import parse_series_matrix
from biomedml.stats import de_genes, coexpression_graph
from biomedml.models import LogRegNP, auc_score

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
auc_plain, auc_graph = [], []
for k in range(5):
    te = folds[k]
    tr = np.concatenate([folds[j] for j in range(5) if j != k])
    _, p_tr = de_genes(X[tr], y[tr])
    sel = np.argsort(p_tr)[:100]
    # plain
    lr = LogRegNP(100, lr=0.2, epochs=800, l2=1e-2, seed=k)
    lr.fit(X[tr][:, sel], y[tr])
    auc_plain.append(auc_score(y[te], X[te][:, sel] @ (lr.W[:, 1] - lr.W[:, 0])))
    # graph-smoothed: build co-expression graph on train fold only
    A, _ = coexpression_graph(X[tr], sel, thresh=0.5)
    A_hat = A + np.eye(len(sel), dtype=np.float32)
    d = A_hat.sum(1)
    A_hat = A_hat / np.sqrt(np.outer(d, d) + 1e-8)
    alpha = 0.5  # half original + half neighborhood signal
    Xg_tr = (1 - alpha) * X[tr][:, sel] + alpha * (X[tr][:, sel] @ A_hat)
    Xg_te = (1 - alpha) * X[te][:, sel] + alpha * (X[te][:, sel] @ A_hat)
    lr2 = LogRegNP(100, lr=0.2, epochs=800, l2=1e-2, seed=k)
    lr2.fit(Xg_tr, y[tr])
    auc_graph.append(auc_score(y[te], Xg_te @ (lr2.W[:, 1] - lr2.W[:, 0])))
    print("fold", k, "plain", round(auc_plain[-1], 3), "graph", round(auc_graph[-1], 3), flush=True)

res = dict(arm="coexpression-graph-smoothed (alpha=0.5, train-fold graph, |r|>=0.5)",
           auc_plain=[float(np.mean(auc_plain)), float(np.std(auc_plain))],
           auc_graph=[float(np.mean(auc_graph)), float(np.std(auc_graph))],
           graph_helps=bool(np.mean(auc_graph) > np.mean(auc_plain) + 0.01))
with open(os.path.join(os.path.dirname(__file__), "..", "results",
                       "graph_arm_results.json"), "w") as f:
    json.dump(res, f, indent=1)
print("GRAPH ARM", res, flush=True)
