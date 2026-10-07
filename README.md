# mega27-16-biodataset-ml-treatment
ML on real biological data: relapse prediction from tumor expression (GEO GSE2034, 286 patients, 22,283 genes, label: bone relapse).
- DE analysis: Welch t + BH-FDR (753 genes at FDR<0.05), fold-contained feature selection (no leakage).
- First pass (70/30 split): CNN 0.640 / logreg 0.535 vs majority 0.767 - honest negative, preserved in results.json.
- Pivot (5-fold CV, top-100 panel, GSE2034 only): logreg AUC 0.632+/-0.055, CNN acc 0.724+/-0.072, both accuracies below the majority baseline 0.759. Within-cohort AUC only; it did not replicate.
- Locked negatives (`results/LOCKED_NEGATIVE_REGISTER.md`): METABRIC replication (n=1979) directional AUC 0.5117, one-sided 95% lower bound 0.4778, primary test false; SCAN-B (overall-survival arm, not relapse) AUC 0.5075, secondary false; two-cohort discovery: zero genes meet joint stability and sign agreement; permutation arm: observed 268 stable genes vs null median 426, panel-size p=0.98 (the register notes this null is not an identical-pipeline test). The completed two-cohort run used default L2 where the preregistration requires L1: a protocol deviation, not a faithful L1 falsification. No relapse signal is claimed beyond GSE2034 cross-validation.
Run: `pip install -e . && pytest && python experiments/run_gse2034.py && python experiments/pivot_cv.py`
