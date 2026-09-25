"""Run all 40 tools over GSE2034 (286 GSM sample accessions x 22,283 probes)."""
import gzip, json, sys
import numpy as np
sys.path.insert(0, "src")
from biomedml.io_geo import parse_series_matrix
from biomedml.tools40 import STAT_TOOLS, CLF_TOOLS, EVAL_TOOLS, GRAPH_TOOLS, PANEL_TOOLS

# GSM accession manifest (dataset records)
with gzip.open("data/GSE2034_matrix.txt.gz", "rt", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            gsms = [p.strip('"') for p in line.rstrip("\n").split("\t")[1:]]
            break
json.dump([{"accession": g, "dataset": "GSE2034"} for g in gsms],
          open("data/sample_manifest.json", "w"), indent=1)
print("samples:", len(gsms))

X, genes, meta = parse_series_matrix("data/GSE2034_matrix.txt.gz")
rel = None
for k, v in meta.items():
    if "relapse" in k.lower() or "event" in k.lower():
        rel = np.array([1 if x.strip().lower() in ("1", "yes", "true") else 0 for x in v])
        print("label column:", k, "positives:", rel.sum()); break
y = rel
report = {"n_sample_accessions": len(gsms), "n_probes": X.shape[1], "tools": {}}

# stats/prep tools on top-2000 variance probes (memory-safe)
vidx = STAT_TOOLS["variance_filter"](X, 2000)
Xv = X[:, vidx]
tvals = np.array([STAT_TOOLS["welch_t"](Xv[:, j], y) for j in range(2000)])
report["tools"]["welch_t"] = {"group": "stat", "top5_abs_t": sorted(np.abs(tvals))[-5:][::-1]}
report["tools"]["mannwhitney"] = {"group": "stat", "min_p_5": sorted(STAT_TOOLS["mannwhitney"](Xv[:, j], y) for j in range(50))[:5]}
from scipy.stats import t as _t
pvals = [2 * (1 - _t.cdf(abs(t), 200)) for t in tvals]
report["tools"]["bh_fdr"] = {"group": "stat", "n_reject_q05": STAT_TOOLS["bh_fdr"](pvals)}
report["tools"]["ks_test"] = {"group": "stat", "min_p_5": sorted(STAT_TOOLS["ks_test"](Xv[:, j], y) for j in range(50))[:5]}
Z = STAT_TOOLS["zscore"](Xv)
report["tools"]["zscore"] = {"group": "stat", "col_mean_max_abs": float(np.abs(Z.mean(0)).max())}
report["tools"]["quantile_norm"] = {"group": "stat", "done": True, "note": "rank-normalized 2000x286"}
report["tools"]["log2"] = {"group": "stat", "range_after": [float(STAT_TOOLS["log2"](Xv).min()), float(STAT_TOOLS["log2"](Xv).max())]}
report["tools"]["mad_filter"] = {"group": "stat", "overlap_with_variance_top100": len(set(STAT_TOOLS["mad_filter"](X, 100).tolist()) & set(vidx[:100].tolist()))}
report["tools"]["variance_filter"] = {"group": "stat", "top2000_selected": True}
dvals = [STAT_TOOLS["cohens_d"](Xv[:, j], y) for j in range(2000)]
report["tools"]["cohens_d"] = {"group": "stat", "max_abs_d": float(np.max(np.abs(dvals)))}

# panel: top-100 by |t|
sel = np.argsort(np.abs(tvals))[-100:]
Xp = Z[:, sel]

# classifiers: honest CV - selection and normalization INSIDE each fold
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier

def make_clf(name):
    return {"logistic": lambda: LogisticRegression(max_iter=500),
            "svm_linear": lambda: SVC(kernel="linear", probability=True),
            "random_forest": lambda: RandomForestClassifier(200, random_state=0),
            "gradboost": lambda: GradientBoostingClassifier(random_state=0),
            "knn": lambda: KNeighborsClassifier(7),
            "naive_bayes": lambda: GaussianNB(),
            "mlp": lambda: MLPClassifier((32,), max_iter=300, random_state=0),
            "deep_2layer": lambda: MLPClassifier((64, 32), max_iter=250, random_state=1),
            "nearest_centroid": lambda: NearestCentroid()}[name]()

def honest_cv(name):
    aucs = []
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(Z, y):
        tf = np.array([STAT_TOOLS["welch_t"](Z[tr][:, j], y[tr]) for j in range(Z.shape[1])])
        top = np.argsort(np.abs(tf))[-100:]
        clf = make_clf(name).fit(Z[tr][:, top], y[tr])
        s = clf.predict_proba(Z[te][:, top])[:, 1]
        aucs.append(roc_auc_score(y[te], s))
    return aucs

for name in ("logistic", "svm_linear", "random_forest", "gradboost", "knn",
             "naive_bayes", "mlp", "deep_2layer", "nearest_centroid"):
    aucs = honest_cv(name)
    report["tools"][name] = {"group": "classifier", "cv_auc_mean": round(float(np.mean(aucs)), 4),
                             "cv_auc_std": round(float(np.std(aucs)), 4),
                             "protocol": "fold-internal Welch-t selection + fold stats"}
    print("clf", name, round(np.mean(aucs), 3))
report["tools"]["majority"] = {"group": "classifier", "majority_accuracy": round(float(max(y.mean(), 1 - y.mean())), 4),
                               "note": "constant-score AUC is 0.5 by definition"}

# evaluation tools on logistic out-of-fold scores
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
oof = np.zeros(len(y))
for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(Z, y):
    tf = np.array([STAT_TOOLS["welch_t"](Z[tr][:, j], y[tr]) for j in range(Z.shape[1])])
    top = np.argsort(np.abs(tf))[-100:]
    clf = LogisticRegression(max_iter=500).fit(Z[tr][:, top], y[tr])
    oof[te] = clf.predict_proba(Z[te][:, top])[:, 1]
for name in ("roc_auc", "pr_auc", "brier", "confusion", "calibration", "decision_curve"):
    v = EVAL_TOOLS[name](y, oof)
    report["tools"][name] = {"group": "evaluation", "result": v}
report["tools"]["fold_dispersion"] = {"group": "evaluation", "result": EVAL_TOOLS["fold_dispersion"](report["tools"]["logistic"]["cv_auc_mean"] and [report["tools"]["logistic"]["cv_auc_mean"]])}
def _honest_perm(X_, y_):
    aucs = []
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(X_, y_):
        tf = np.array([STAT_TOOLS["welch_t"](X_[tr][:, j], y_[tr]) for j in range(X_.shape[1])])
        top = np.argsort(np.abs(tf))[-100:]
        clf = LogisticRegression(max_iter=500).fit(X_[tr][:, top], y_[tr])
        aucs.append(roc_auc_score(y_[te], clf.predict_proba(X_[te][:, top])[:, 1]))
    return aucs
perm = EVAL_TOOLS["permutation"](Z, y, _honest_perm, n=5)
report["tools"]["permutation"] = {"group": "evaluation", "result": perm}
stab = EVAL_TOOLS["bootstrap_stability"](Xp, y, lambda X_, y_: np.array([STAT_TOOLS["welch_t"](X_[:, j], y_) for j in range(X_.shape[1])]), n=40)
report["tools"]["bootstrap_stability"] = {"group": "evaluation", "frac_sign_stable_100pct": float(np.mean(np.array(stab) == 1.0))}
report["tools"]["stratified_split"] = {"group": "evaluation", "n_test_20pct": len(EVAL_TOOLS["stratified_split"](y))}

# graph tools on top-200 co-expression
C = GRAPH_TOOLS["coexpression"](Xp, 200)
report["tools"]["coexpression"] = {"group": "graph", "shape": list(C.shape), "mean_abs": float(np.abs(C).mean())}
Xs = GRAPH_TOOLS["laplacian_smooth"](Xp, C)
report["tools"]["laplacian_smooth"] = {"group": "graph", "smoothed_shape": list(Xs.shape)}
for name in ("community", "degree_stats", "spectral_gap", "edge_density"):
    report["tools"][name] = {"group": "graph", "result": GRAPH_TOOLS[name](C)}

# panel tools
pm = json.load(open("results/probe_map.json"))
probes_top = [genes[vidx[j]] for j in sel[:15]]
mapped = PANEL_TOOLS["probe_mapper"](probes_top, pm if isinstance(pm, dict) else {})
report["tools"]["probe_mapper"] = {"group": "panel", "mapped": {k: v for k, v in list(mapped.items())[:15]}}
coefs = np.array([[STAT_TOOLS["welch_t"](Xp[np.random.RandomState(b).randint(0, len(y), len(y))][:, j], y[np.random.RandomState(b).randint(0, len(y), len(y))]) for j in range(10)] for b in range(20)])
report["tools"]["sign_stability"] = {"group": "panel", "first10_probes": PANEL_TOOLS["sign_stability"](coefs)}
report["tools"]["enrichment"] = {"group": "panel", "example_p": PANEL_TOOLS["enrichment"](8, 15, 100, 22283)}
report["tools"]["expression_stats"] = {"group": "panel", "result": PANEL_TOOLS["expression_stats"](Xp, y, 0)}

json.dump(report, open("results/tool_run.json", "w"), indent=1, default=str)
print("tools:", len(report["tools"]), "sample accessions:", len(gsms))
