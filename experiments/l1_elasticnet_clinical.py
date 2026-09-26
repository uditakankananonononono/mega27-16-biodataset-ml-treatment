"""L1 rung (PREREG_IMPROVEMENT.md, exact protocol locked 2026-09-26):
fold-contained elastic-net logreg on all probes + clinical covariates,
5 folds x 5 seeds, paired vs locked Wang76 comparator (gate G1).
Resumable: per-seed results saved; existing seeds skipped."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from biomedml.io_geo import parse_series_matrix
from biomedml.models import auc_score
from sklearn.linear_model import LogisticRegression
import gzip, csv

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUTP = os.path.join(ROOT, "results", "l1_elasticnet.json")
SEEDS = [9, 19, 29, 39, 49]

def load_all():
    X, genes, meta = parse_series_matrix(os.path.join(ROOT, "data", "GSE2034_matrix.txt.gz"))
    with gzip.open(os.path.join(ROOT, "data", "GSE2034_matrix.txt.gz"), "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                gsms = [p.strip('"') for p in line.rstrip("\n").split("\t")[1:]]
                break
    raw = meta["bone relapses (1=yes, 0=no)"]
    keep = [i for i, v in enumerate(raw) if v.strip() in ("0", "1")]
    y = np.array([int(raw[i]) for i in keep])
    kept = [gsms[i] for i in keep]
    Xl = np.log2(np.maximum(X[keep], 0.1))
    clin = {r["accession"]: r for r in json.load(open(os.path.join(ROOT, "data", "gse2034_clinical_vdx.json")))}
    def num(v):
        try: return float(v)
        except (TypeError, ValueError): return np.nan
    age = np.array([num(clin[g].get("age")) for g in kept])
    er = np.array([num(clin[g].get("er")) for g in kept])
    grade = np.array([num(clin[g].get("grade")) for g in kept])
    return Xl, y, np.vstack([age, er, grade]).T, genes, kept

def wang_scores(Xl, y, genes, kept):
    sig = list(csv.DictReader(open(os.path.join(ROOT, "results", "sig_gene76.csv"))))
    probes = [r["probe"] for r in sig]
    coeff = np.array([float(r["std.cox.coefficient"]) for r in sig])
    ersig = np.array([int(r["er"]) for r in sig])
    gset = {g: i for i, g in enumerate(genes)}
    sel = [gset[p] for p in probes]
    clin = {r["accession"]: r for r in json.load(open(os.path.join(ROOT, "data", "gse2034_clinical_vdx.json")))}
    er = np.array([int(clin[g]["er"]) for g in kept])
    W = np.zeros(len(y))
    for arm, const in ((1, 313.5), (0, 280.0)):
        m = er == arm
        am = ersig == arm
        W[m] = const + Xl[m][:, [sel[i] for i in np.where(am)[0]]] @ coeff[am]
    return W

def run_seed(seed, Xl, y, C):
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(y))
    folds = np.array_split(perm, 5)
    aucs, aucs_clin, maj = [], [], []
    for k in range(5):
        te = folds[k]
        tr = np.concatenate([folds[j] for j in range(5) if j != k])
        mu, sd = Xl[tr].mean(0), Xl[tr].std(0) + 1e-6
        Xtr = (Xl[tr] - mu) / sd
        Xte = (Xl[te] - mu) / sd
        cmu = np.nanmean(C[tr], axis=0); csd = np.nanstd(C[tr], axis=0) + 1e-6
        cmed = np.nanmedian(C[tr], axis=0)
        Ctr = np.where(np.isnan(C[tr]), cmed, C[tr]); Cte = np.where(np.isnan(C[te]), cmed, C[te])
        keepc = np.nanstd(Ctr, axis=0) > 0  # drop constant clinical cols (train fold)
        Ctr = (Ctr[:, keepc] - cmu[keepc]) / csd[keepc]
        Cte = (Cte[:, keepc] - cmu[keepc]) / csd[keepc]
        m = LogisticRegression(penalty="elasticnet", solver="saga", C=0.05,
                               l1_ratio=0.8, max_iter=3000, tol=1e-3, random_state=seed)
        m.fit(np.hstack([Xtr, Ctr]), y[tr])
        aucs.append(auc_score(y[te], m.predict_proba(np.hstack([Xte, Cte]))[:, 1]))
        mc = LogisticRegression(C=1.0, max_iter=1000, random_state=seed)
        mc.fit(Ctr, y[tr])
        aucs_clin.append(auc_score(y[te], mc.predict_proba(Cte)[:, 1]))
        maj.append(float(max((y[te] == 0).mean(), (y[te] == 1).mean())))
    return aucs, aucs_clin, maj

def main():
    res = json.load(open(OUTP)) if os.path.exists(OUTP) else {"per_seed": {}, "wang_fold_aucs": None}
    Xl, y, C, genes, kept = load_all()
    W = wang_scores(Xl, y, genes, kept)
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        perm = rng.permutation(len(y))
        folds = np.array_split(perm, 5)
        w_aucs = [auc_score(y[te], W[te]) for te in folds]
        if res["wang_fold_aucs"] is None:
            res["wang_fold_aucs"] = {}
        res["wang_fold_aucs"][str(seed)] = [float(a) for a in w_aucs]
        if str(seed) in res["per_seed"]:
            print("seed", seed, "cached", flush=True)
            continue
        aucs, aucs_clin, maj = run_seed(seed, Xl, y, C)
        res["per_seed"][str(seed)] = {"l1": [float(a) for a in aucs],
                                      "clinical_only": [float(a) for a in aucs_clin],
                                      "majority": [float(a) for a in maj]}
        json.dump(res, open(OUTP, "w"), indent=1)
        print("seed", seed, "l1", [round(a, 3) for a in aucs], flush=True)
    json.dump(res, open(OUTP, "w"), indent=1)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
