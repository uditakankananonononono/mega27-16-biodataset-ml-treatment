# mega27-16-biodataset-ml-treatment
ML on real biological data: relapse prediction from tumor expression (GEO GSE2034, 286 patients, 22,283 genes, label: bone relapse).
- DE analysis: Welch t + BH-FDR (753 genes at FDR<0.05), fold-contained feature selection (no leakage).
- First pass (70/30 split): CNN 0.640 / logreg 0.535 vs majority 0.767 - honest negative, preserved in results.json.
- Pivot (5-fold CV, top-100 panel): logreg AUC 0.632+/-0.055, CNN acc 0.724+/-0.072 vs majority 0.759 - weak but real signal above chance; verdict and both passes documented.
Run: `pip install -e . && pytest && python experiments/run_gse2034.py && python experiments/pivot_cv.py`
