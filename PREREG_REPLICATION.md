# PREREG_REPLICATION — Locked external-replication rung (lane-16)

Locked: 2026-09-27, BEFORE any METABRIC / SCAN-B / additional replication
analysis on the panel below. Amendment queue A2. Dated commit precedes all
replication computation; any post-lock change is declared in the paper as a
protocol deviation, not silently absorbed.

## 1. Frozen panel (selected on GSE2034 ONLY)

Bootstrap sign-stability selection (n_boot=40, L1 logistic, relapse-vs-control,
GSE2034 discovery cohort). 30 probes passed stability >= 0.975; the 15 that
map to a unique gene symbol form the locked 15-gene panel. The 15 unmapped
probes are carried as secondary probe-level features and are NOT part of the
primary replication test.

| probe | gene | direction (discovery) | stability | mean coef |
|---|---|---|---|---|
| 222013_x_at | EEF2KMT | protective | 1.000 | -0.1151 |
| 219478_at | WFDC1 | risk | 1.000 | 0.0989 |
| 200876_s_at | PSMB1 | protective | 1.000 | -0.0967 |
| 210019_at | CALML3 | protective | 1.000 | -0.0906 |
| 221348_at | NPPC | protective | 1.000 | -0.0892 |
| 206594_at | PASK | risk | 1.000 | 0.0869 |
| 219445_at | GLTSCR1 | risk | 1.000 | 0.0855 |
| 32088_at | BLZF1 | risk | 1.000 | 0.0835 |
| 216862_s_at | CMC4 | protective | 1.000 | -0.0824 |
| 214277_at | COX11 | risk | 0.975 | 0.0818 |
| 204726_at | CDH13 | risk | 1.000 | 0.0788 |
| 214046_at | FUT9 | protective | 1.000 | -0.0780 |
| 202053_s_at | ALDH3A2 | protective | 1.000 | -0.0769 |
| 218787_x_at | CWF19L1 | protective | 1.000 | -0.0768 |
| 219756_s_at | POF1B | risk | 1.000 | 0.0767 |

Panel coefficients are REFIT inside each replication cohort's training-free
scoring rule: locked score = sum over genes of sign(coef_discovery) *
z(expression). No coefficient magnitude, threshold, or gene membership is
re-estimated on replication data. This is the single locked scoring rule.

## 2. Frozen cohorts and endpoints

- R1 (primary): METABRIC (via cBioPortal public API if freely accessible;
  access attempt documented either way). Endpoint: relapse / distant
  metastasis-free survival as coded by the METABRIC clinical annotation
  (DSS/relapse status field chosen and named in the analysis script header
  before the script is run).
- R2 (secondary): SCAN-B if freely accessible without application; if gated,
  document the gating honestly and substitute the strongest free cohort
  found, named here-in-spirit by an appended dated note (substitution is a
  declared protocol event, not a silent swap).
- R3 (fold-in, PRIOR-PLANNED): GSE2603 bone-metastasis arm (already run
  pre-verdict @6d2fffc; re-scored under this locked rule for comparability;
  not claimed as verdict novelty).

Harmonization: gene-symbol mapping per cohort platform; samples restricted to
primary breast tumours; batch/scale handled by within-cohort z-scoring only.

## 3. Locked statistical tests (one primary, two secondary)

- PRIMARY TEST: directional AUC of the locked score on METABRIC relapse
  endpoint. Success = AUC > 0.5 with one-sided 95% bootstrap CI (cluster-level
  where samples cluster, e.g. by hospital/site) excluding 0.5, AND mean
  per-gene sign concordance with discovery directions > 50%.
- SECONDARY 1: same locked test on SCAN-B (or declared substitute).
- SECONDARY 2: GSE2603 arm under the locked rule.
- No endpoint, threshold, gene, or direction changes after this lock.

## 4. Falsification gate

If the PRIMARY TEST fails (AUC CI includes 0.5 or sign concordance <= 50%):
the panel is declared NOT externally validated. The paper states this plainly
in abstract, results, and limitations; the stability-selection method claim
reframes to "discovery-cohort-internal stability only"; and the negative
result is preserved with full numbers. A failure does NOT trigger panel
re-selection on replication data (that would be double-dipping and is
explicitly out of scope).

## 5. Honesty notes

- Expert adjudication of labels was considered and WAIVED by the owner
  (2026-09-27, recorded in amendment queue); residual label-noise risk is
  stated in both papers.
- AUC ~0.632 discovery-cohort band is the method's current ceiling and is
  reported as such, not as clinical utility.
