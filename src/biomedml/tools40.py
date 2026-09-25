"""40 named tools for item 16 (GSE2034 relapse ML). Hermetic, no network.
Groups: 10 stats/prep, 10 classifiers, 10 evaluation, 6 graph, 4 panel."""
import math
import numpy as np

# ---- stats / preprocessing (10) ----
def t_welch(x, ybin):
    a, b = x[ybin == 1], x[ybin == 0]
    t = (a.mean() - b.mean()) / np.sqrt(a.var(ddof=1)/len(a) + b.var(ddof=1)/len(b))
    return float(t)

def t_mannwhitney(x, ybin):
    from scipy.stats import mannwhitneyu
    return float(mannwhitneyu(x[ybin == 1], x[ybin == 0]).pvalue)

def t_bh_fdr(pvals, q=0.05):
    p = np.asarray(pvals); m = len(p); order = np.argsort(p)
    thresh = (np.arange(1, m + 1) / m) * q
    passed = p[order] <= thresh
    k = passed.sum() if not passed.any() else np.max(np.where(passed)) + 1
    return int(k)

def t_ks_test(x, ybin):
    from scipy.stats import ks_2samp
    return float(ks_2samp(x[ybin == 1], x[ybin == 0]).pvalue)

def t_zscore(X):
    return (X - X.mean(0)) / (X.std(0) + 1e-12)

def t_quantile_norm(X):
    from scipy.stats import rankdata
    R = np.apply_along_axis(rankdata, 0, X)
    sm = np.sort(X, axis=0).mean(1)
    out = np.empty_like(X)
    for j in range(X.shape[1]):
        out[:, j] = sm[(R[:, j] - 1).astype(int)]
    return out

def t_log2(X): return np.log2(np.maximum(X, 1e-6))

def t_mad_filter(X, top=1000):
    mad = np.median(np.abs(X - np.median(X, 0)), 0)
    return np.argsort(mad)[-top:]

def t_variance_filter(X, top=1000):
    return np.argsort(X.var(0))[-top:]

def t_cohens_d(x, ybin):
    a, b = x[ybin == 1], x[ybin == 0]
    s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / (s + 1e-12))

STAT_TOOLS = {"welch_t": t_welch, "mannwhitney": t_mannwhitney, "bh_fdr": t_bh_fdr,
 "ks_test": t_ks_test, "zscore": t_zscore, "quantile_norm": t_quantile_norm,
 "log2": t_log2, "mad_filter": t_mad_filter, "variance_filter": t_variance_filter,
 "cohens_d": t_cohens_d}

# ---- classifiers (10) ----
def _cv(clf, X, y, k=5):
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    return cross_val_score(clf, X, y, cv=StratifiedKFold(k, shuffle=True, random_state=0),
                           scoring="roc_auc").tolist()

def c_logistic(X, y):
    from sklearn.linear_model import LogisticRegression
    return _cv(LogisticRegression(max_iter=500), X, y)

def c_svm(X, y):
    from sklearn.svm import SVC
    return _cv(SVC(kernel="linear", probability=True), X, y)

def c_random_forest(X, y):
    from sklearn.ensemble import RandomForestClassifier
    return _cv(RandomForestClassifier(200, random_state=0), X, y)

def c_gradboost(X, y):
    from sklearn.ensemble import GradientBoostingClassifier
    return _cv(GradientBoostingClassifier(random_state=0), X, y)

def c_knn(X, y):
    from sklearn.neighbors import KNeighborsClassifier
    return _cv(KNeighborsClassifier(7), X, y)

def c_naive_bayes(X, y):
    from sklearn.naive_bayes import GaussianNB
    return _cv(GaussianNB(), X, y)

def c_mlp(X, y):
    from sklearn.neural_network import MLPClassifier
    return _cv(MLPClassifier((32,), max_iter=300, random_state=0), X, y)

def c_cnn1d(X, y):
    from sklearn.neural_network import MLPClassifier
    return _cv(MLPClassifier((64, 32), max_iter=250, random_state=1), X, y)  # deep proxy

def c_majority(X, y):
    p = max(y.mean(), 1 - y.mean())
    return [float(p)] * 5

def c_centroid(X, y):
    from sklearn.neighbors import NearestCentroid
    return _cv(NearestCentroid(), X, y)

CLF_TOOLS = {"logistic": c_logistic, "svm_linear": c_svm, "random_forest": c_random_forest,
 "gradboost": c_gradboost, "knn": c_knn, "naive_bayes": c_naive_bayes, "mlp": c_mlp,
 "deep_2layer": c_cnn1d, "majority": c_majority, "nearest_centroid": c_centroid}

# ---- evaluation (10) ----
def e_roc_auc(y, s):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s))

def e_pr_auc(y, s):
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(y, s))

def e_brier(y, s):
    from sklearn.metrics import brier_score_loss
    return float(brier_score_loss(y, s))

def e_confusion(y, s, t=0.5):
    from sklearn.metrics import confusion_matrix
    return confusion_matrix(y, (s > t).astype(int)).tolist()

def e_calibration(y, s, bins=6):
    edges = np.linspace(0, 1, bins + 1); out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (s >= a) & (s < b)
        out.append(float(y[m].mean()) if m.sum() else None)
    return out

def e_fold_dispersion(aucs):
    return {"mean": float(np.mean(aucs)), "std": float(np.std(aucs))}

def e_permutation(X, y, clf_fn, n=15):
    rng = np.random.RandomState(0)
    obs = np.mean(clf_fn(X, y))
    cnt = sum(np.mean(clf_fn(X, rng.permutation(y))) >= obs for _ in range(n))
    return {"observed": float(obs), "p": (cnt + 1) / (n + 1)}

def e_bootstrap_stability(X, y, select_fn, n=40):
    rng = np.random.RandomState(0)
    signs = []
    for _ in range(n):
        idx = rng.randint(0, len(y), len(y))
        signs.append(np.sign(select_fn(X[idx], y[idx])))
    signs = np.array(signs)
    med = np.sign(np.median(signs, 0))
    return ((signs == med).mean(0)).tolist()

def e_decision_curve(y, s, ts=(0.1, 0.2, 0.3)):
    out = {}
    for t in ts:
        pred = s > t
        tp = ((pred == 1) & (y == 1)).sum(); fp = ((pred == 1) & (y == 0)).sum()
        out[str(t)] = float(tp / len(y) - fp / len(y) * (t / (1 - t)))
    return out

def e_stratified_split(y, frac=0.2, seed=0):
    rng = np.random.RandomState(seed)
    te = []
    for c in np.unique(y):
        idx = np.where(y == c)[0]; rng.shuffle(idx)
        te += list(idx[: int(len(idx) * frac)])
    return sorted(te)

EVAL_TOOLS = {"roc_auc": e_roc_auc, "pr_auc": e_pr_auc, "brier": e_brier,
 "confusion": e_confusion, "calibration": e_calibration,
 "fold_dispersion": e_fold_dispersion, "permutation": e_permutation,
 "bootstrap_stability": e_bootstrap_stability, "decision_curve": e_decision_curve,
 "stratified_split": e_stratified_split}

# ---- graph (6) ----
def g_coexpression(X, top=200):
    C = np.corrcoef(X[:, :top].T)
    return C

def g_laplacian_smooth(X, C, alpha=0.3):
    W = np.clip(C, 0, None) ** 8
    D = np.diag(W.sum(1) + 1e-12)
    L = np.linalg.inv(D) @ W
    return (1 - alpha) * X + alpha * X @ L.T

def g_community(C, thresh=0.3):
    A = (np.abs(C) > thresh).astype(int); np.fill_diagonal(A, 0)
    n = A.shape[0]; seen = np.zeros(n, bool); sizes = []
    for s in range(n):
        if seen[s]: continue
        stack, seen[s], c = [s], True, 0
        while stack:
            u = stack.pop(); c += 1
            for w in np.where(A[u])[0]:
                if not seen[w]: seen[w] = True; stack.append(w)
        sizes.append(c)
    return {"n_components": len(sizes), "largest": max(sizes)}

def g_degree_stats(C, thresh=0.3):
    A = (np.abs(C) > thresh).astype(int); np.fill_diagonal(A, 0)
    d = A.sum(1)
    return {"mean_degree": float(d.mean()), "max_degree": int(d.max())}

def g_spectral_gap(C):
    W = np.clip(C, 0, None)
    D = np.diag(W.sum(1) + 1e-12)
    L = np.eye(len(W)) - np.linalg.inv(np.sqrt(D)) @ W @ np.linalg.inv(np.sqrt(D))
    ev = np.linalg.eigvalsh(L)
    return float(ev[1] - ev[0])

def g_edge_density(C, thresh=0.3):
    A = (np.abs(C) > thresh); np.fill_diagonal(A, False)
    n = len(A)
    return float(A.sum() / (n * (n - 1)))

GRAPH_TOOLS = {"coexpression": g_coexpression, "laplacian_smooth": g_laplacian_smooth,
 "community": g_community, "degree_stats": g_degree_stats,
 "spectral_gap": g_spectral_gap, "edge_density": g_edge_density}

# ---- panel / mapping (4) ----
def p_probe_mapper(probes, table):
    return {p: table.get(p) for p in probes}

def p_sign_stability(coefs):
    med = np.sign(np.median(coefs, 0))
    return ((np.sign(coefs) == med).mean(0)).tolist()

def p_enrichment(overlap, panel, pool, universe):
    from scipy.stats import hypergeom
    return float(1 - hypergeom.cdf(overlap - 1, universe, panel, pool))

def p_expression_stats(X, ybin, idx):
    return {"mean_relapse": float(X[ybin == 1][:, idx].mean()),
            "mean_nonrelapse": float(X[ybin == 0][:, idx].mean())}

PANEL_TOOLS = {"probe_mapper": p_probe_mapper, "sign_stability": p_sign_stability,
 "enrichment": p_enrichment, "expression_stats": p_expression_stats}

ALL_TOOLS = {**STAT_TOOLS, **CLF_TOOLS, **EVAL_TOOLS, **GRAPH_TOOLS, **PANEL_TOOLS}
assert len(ALL_TOOLS) == 40, len(ALL_TOOLS)
