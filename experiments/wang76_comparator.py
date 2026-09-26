"""Comparator A (preregistered PREREG_IMPROVEMENT.md): Wang et al. 2005 Lancet
76-gene relapse score, re-implemented exactly per genefu::gene76 scoring
(constants 313.5 ER+ / 280 ER-, std.cox.coefficients from sig.gene76),
evaluated against the GEO bone-relapse label on the lane's identical folds
(seed 9 permutation, 5 folds, same as pivot_cv.py).

Inputs:
- data/GSE2034_matrix.txt.gz (bundled; verified = log2(clip(X,0.1)) of the
  same samples redistributed in breastCancerVDX, expression identity r=1.0)
- data/gse2034_clinical_vdx.json (ER etc. per GSM; mapping in
  results/vdx_gsm_map.json, r=1.0 identity on 22,283 common probes)
- results/sig_gene76.csv (probe, std.cox.coefficient, er arm) from genefu data
Outputs: results/wang76_comparator.json
"""
import json, os, sys, gzip
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from biomedml.io_geo import parse_series_matrix
from biomedml.models import auc_score

ROOT = os.path.join(os.path.dirname(__file__), "..")
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

import csv
sig = list(csv.DictReader(open(os.path.join(ROOT, "results", "sig_gene76.csv"))))
probes = [r["probe"] for r in sig]
coeff = np.array([float(r["std.cox.coefficient"]) for r in sig])
ersig = np.array([int(r["er"]) for r in sig])
gset = {g: i for i, g in enumerate(genes)}
missing = [p for p in probes if p not in gset]
assert not missing, f"Wang probes missing: {missing}"

clin = {r["accession"]: r for r in json.load(open(os.path.join(ROOT, "data", "gse2034_clinical_vdx.json")))}
er = np.array([int(clin[g]["er"]) for g in kept])
Xl = np.log2(np.maximum(X[keep], 0.1))
sel = [gset[p] for p in probes]
W = np.zeros(len(y))
for arm, const in ((1, 313.5), (0, 280.0)):
    m = er == arm
    am = ersig == arm
    W[m] = const + Xl[m][:, [sel[i] for i in np.where(am)[0]]] @ coeff[am]

rng = np.random.default_rng(9)
perm = rng.permutation(len(y))
folds = np.array_split(perm, 5)
fold_aucs = [auc_score(y[te], W[te]) for te in folds]
res = {
    "comparator": "Wang2005 gene76 relapse score (genefu sig.gene76 coefficients, ER-conditional 313.5/280)",
    "n": int(len(y)), "er_plus": int((er == 1).sum()), "probes_used": len(probes),
    "fold_aucs": [float(a) for a in fold_aucs],
    "pooled_auc": float(auc_score(y, W)),
    "fold_protocol": "seed9 permutation, 5 folds (identical to pivot_cv.py)",
    "er_source": "breastCancerVDX vdx clinical via r=1.0 expression-identity map",
}
json.dump(res, open(os.path.join(ROOT, "results", "wang76_comparator.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
