# Source-based inventory audit, October 1
286 GSM sample accessions are samples of one cohort, not independent datasets;15 genes are records, not datasets. 301 records cannot certify120 datasets.
experiments/tool_inventory.py: quantile_norm only reports done; no call. honest_cv uses globally variance-filtered and globally z-scored Z. Welch-t selection is fold-internal, normalization is not. fold_dispersion receives one mean value, not recorded folds. Do not certify those outputs as fully fold-internal, valid fold dispersion or40 executed functions. Preserved old tool_run.json unchanged. No new fits.
External tool_run inventory needs independent execution/provenance review; installed or status=ok does not alone certify scientific use.
