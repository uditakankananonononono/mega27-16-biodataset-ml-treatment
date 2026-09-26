# MEGA27-16 - GSE2034 bone-relapse - IMPROVEMENT PREREGISTRATION (locked before outcomes)
Current state (verified results/*.json): prevalence 0.2413; pass-1 70/30 CNN 0.640/logreg 0.535 vs majority
0.767 (honest negative, PRESERVED); pivot 5-fold CV logreg AUC 0.632+/-0.055, CNN acc 0.724+/-0.072 vs
majority 0.759 - "signal-above-chance" but NOT a benchmark beat.

## Reframed benchmark (the majority baseline is a floor, not a benchmark)
- COMPARATOR A (primary): Wang et al. 2005 (Lancet) 76-gene signature score - the published same-dataset
  same-task model, re-implemented and evaluated on OUR identical folds. Europe PMC 15721472.
- COMPARATOR B (sanity): majority-class, logreg-on-DE-100 (current best), clinical-only (age/ER).
- LOCKED QUESTION Q1: does the improved model beat Comparator A on bone-relapse AUC under identical folds?
- LOCKED GATE G1: mean AUC diff > 0, paired bootstrap 95% CI excludes 0, 5 folds x 5 seeds.

## Improvement ladder (in order; stop at first gate pass, keep going for discovery)
L1 fold-contained elastic-net + clinical covariates; L2 SVM-RFE / GBM (XGBoost);
L3 pathway-level features (MSigDB Hallmark) to kill probe-level noise; L4 calibrated ensemble.

## NEW DISCOVERY gate (locked)
D1: a novel bone-relapse gene panel that (i) passes G1 AND (ii) replicates directionally (AUC>0.5, same-sign
effect) on an INDEPENDENT cohort with bone-metastasis labels (candidates: GSE14020, GSE2603, GSE5327 -
verify label fields before locking final choice). External replication = the discovery. If all fail:
RULE 6 - ask ChatGPT for redirection, pivot on strongest, verbatim log.
