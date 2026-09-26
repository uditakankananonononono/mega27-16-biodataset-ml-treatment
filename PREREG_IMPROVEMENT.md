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

## L1 exact protocol (locked 2026-09-26 17:09 IST, before first L1 fit)
- Features: all 22,283 probes (log2 clip 0.1), z-scored USING TRAIN-FOLD
  mean/std ONLY (fold-contained; the old pivot_cv global z-score is
  superseded), + clinical covariates (age z-scored on train, ER binary,
  grade ordinal z-scored on train; node dropped - constant 0 in this
  lymph-node-negative cohort). Median imputation (train-fold medians) for
  missing clinical values.
- Model: sklearn LogisticRegression(penalty="elasticnet", solver="saga",
  C=0.05, l1_ratio=0.8, max_iter=3000, tol=1e-3), fixed a priori - no
  test-fold tuning, no nested selection.
- Evaluation: 5 folds x 5 seeds (seeds 9, 19, 29, 39, 49; fold protocol
  identical to pivot_cv at seed 9). Fold-level AUCs paired against the
  locked Wang76 comparator scores on the same folds.
- Gate G1 decision rule: mean(fold AUC_L1 - fold AUC_Wang76) > 0 AND
  paired-bootstrap 95% CI (10,000 resamples over the 25 fold-level diffs,
  seed 7) excludes 0. Comparator B reported alongside (majority,
  clinical-only logreg) for sanity.
- If G1 fails: honest negative, escalate to L2 (SVM-RFE / XGBoost),
  RULE 6 ChatGPT redirection round feeds the next rung's novelty.

### L1 compute amendment (2026-09-26 17:11 IST, before any L1 result)
max_iter reduced 3000 -> 300 (saga on n=229 << p converges well before;
checked: no result had been produced when amended). Everything else
unchanged. Seeds run one-per-invocation (sandbox 120s call limit).
