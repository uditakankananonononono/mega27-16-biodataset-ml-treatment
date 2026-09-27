# PREREG_TWOCOHORT — locked two-cohort discovery design (lane-16 pivot, redirection candidate 1)

Locked: 2026-09-27, BEFORE any cross-cohort selection is run. This is the
pivot ordered after the locked-gate falsification of the single-cohort panel
(METABRIC @97ecc5d, SCAN-B @e1bbdc7). Direction chosen by owner via main:
candidate 1 (two-cohort discovery). Candidate 2 (subtype-stratified panel) is
documented in the paper as the considered alternative, NOT executed now.

Lesson carried (locked evidence): sign-stability inside ONE cohort is not a
generalization certificate; the new rule requires cross-cohort stability
BEFORE any test cohort is touched.

## 1. Cohorts and roles (frozen)

- DISCOVERY (selection only): GSE2034 (Affymetrix U133A, relapse endpoint,
  as in the original pipeline) + METABRIC brca_metabric (Illumina HT-12 v3,
  RFS_STATUS endpoint, cBioPortal pull of 2026-09-27).
- TEST (never used for any selection decision): PRIMARY = SCAN-B GSE96058
  (OS-endpoint arm, declared fallback per PREREG_REPLICATION.md section 2 +
  REPLICATION_ACCESS.md dated note); SECONDARY = GSE2603 bone-metastasis arm.
- DECLARED PRIOR USE: SCAN-B and GSE2603 were previously scored with the OLD
  panel (committed results). No selection decision for the NEW panel uses any
  test-cohort data; their role as untouched test arms for the new panel is
  preserved by construction and this declaration.

## 2. Frozen selection rule (runs on discovery cohorts only)

Per cohort, independently:
1. Map to gene symbols (GSE2034: probe -> symbol via the existing probe_map;
   multi-probe symbols keep the max-|t| probe; METABRIC: symbols native).
2. Within-cohort z-score; prefilter top 2000 symbols by |t| vs the cohort
   endpoint (same rule as the original pipeline).
3. 40 bootstrap L1 logistic fits (C=0.5, max_iter=500, StratifiedKFold(5,
   shuffle, seed 0) rotating splits, bootstrap resample of the training fold,
   rng seed 0) - identical recipe to experiments/gene_stability.py.
4. Gene stability = fraction of 40 fits whose coefficient keeps the
   majority sign.

A gene JOINS THE PANEL iff: stability >= 0.975 in BOTH cohorts AND the
majority signs agree across cohorts. If >30 genes pass, keep the top 30 by
min(stability_GSE2034, stability_METABRIC) * |mean coef sign product|; if 0
pass, the design is declared failed as designed (no threshold loosening -
a loosened threshold would be a NEW prereg).

## 3. Frozen scoring and tests

Score rule (unchanged shape): sum over panel genes of sign * z(expression),
signs = cross-cohort majority sign; within-cohort z on the test cohort;
NO coefficient or membership estimation on test data.

- PRIMARY TEST: AUC vs SCAN-B OS event; success = one-sided 95% bootstrap CI
  (2000 patient-level resamples, seed 2026) excludes 0.5 AND per-gene sign
  concordance > 50%.
- SECONDARY: same rule vs GSE2603 bone-metastasis event.

## 4. Falsification gate

Primary fails => the two-cohort panel is declared NOT validated, with full
numbers in abstract/results/limitations. No re-selection on test cohorts, no
threshold loosening, no endpoint switching. A further redesign requires a new
owner decision and a new prereg.

## 5. Honesty notes

- Both falsification results of the single-cohort panel stay prominent
  (abstract/results/limitations) per owner directive 12:56 PM.
- Compute: shared 2-core box; selection runs interleave with locked grind
  work; runtime does not affect the locked rule.
