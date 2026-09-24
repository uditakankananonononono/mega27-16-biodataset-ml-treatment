"""Per-gene statistics for the stability panel + top-100 DE table."""
import json, sys
import numpy as np
sys.path.insert(0, "src")
from biomedml.io_geo import parse_series_matrix
from biomedml.stats import de_genes

X, genes, meta = parse_series_matrix("data/GSE2034_matrix.txt.gz")
key = [k for k in meta if "relapse" in k.lower()][0]
y = np.array([1 if v.strip() == "1" else 0 for v in meta[key]])
X = np.nan_to_num(X)
t, p = de_genes(X, y)
order = np.argsort(-np.abs(t))[:100]
top = [dict(probe=genes[i], t=round(float(t[i]), 3),
            mean_relapse=round(float(X[y == 1, i].mean()), 3),
            mean_control=round(float(X[y == 0, i].mean()), 3)) for i in order]
json.dump(top, open("results/top100_de.json", "w"), indent=1)
panel = json.load(open("results/gene_stability.json"))["panel"]
rows = []
for g in panel:
    i = genes.index(g["probe"])
    rows.append(dict(probe=g["probe"], t=round(float(t[i]), 3),
                     mean_relapse=round(float(X[y == 1, i].mean()), 3),
                     mean_control=round(float(X[y == 0, i].mean()), 3),
                     rank_by_t=int((np.abs(t) > abs(t[i])).sum()) + 1))
json.dump(rows, open("results/panel_stats.json", "w"), indent=1)
print("panel t-stats:", [(r["probe"], r["t"], r["rank_by_t"]) for r in rows[:5]])
