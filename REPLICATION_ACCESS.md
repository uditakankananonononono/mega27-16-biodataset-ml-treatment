# REPLICATION_ACCESS — access investigation (queue step 2), dated 2026-09-27

## METABRIC (R1, primary) — FREE ACCESS CONFIRMED
- cBioPortal public API, study `brca_metabric`, publicStudy=true
  (Nature 2012 / Nat Commun 2016; 2,509 primary breast tumours).
- mRNA expression: `brca_metabric_mrna` (Illumina HT-12 v3 microarray,
  1,980 samples). Z-score profile also available.
- Locked endpoint available verbatim: `RFS_STATUS` / `RFS_MONTHS`
  (Relapse Free Status) — matches the frozen prereg endpoint without
  re-definition. OS_MONTHS/OS_STATUS also present for secondary arms.
- Panel gene availability (cBioPortal gene index, Entrez):
  14/15 resolve directly; GLTSCR1 is an alias of BICRA (current symbol).
  Actual in-profile presence verified at pull time; any gene absent from
  the HT-12 profile is dropped from the score WITH a dated declaration
  (score uses remaining genes; no substitution).
- API spot checks (2026-09-27): study metadata, molecular profiles,
  clinical attributes, and gene endpoints all returned live.

## SCAN-B (R2, secondary) — FREE ACCESS CONFIRMED (no application)
- GEO `GSE96058`, public since 2018-03-12: SCAN-B (NCT02306096)
  population-based cohort, 3,273 primary breast cancers, RNA-seq,
  median follow-up 52 months. Downloadable directly from GEO FTP/HTTPS.
- Endpoint: OS / RFI fields in the associated clinical annotation
  (exact field names fixed in the analysis script header before run,
  per prereg section 2; if no relapse field exists in the GEO clinical
  file, SCAN-B is scored on OS and declared an OS-endpoint arm, not
  silently compared to RFS arms).

## GSE2603 (R3, PRIOR-PLANNED fold-in)
- Already run pre-verdict (@6d2fffc); re-score under locked rule only.

## Next
Build `src/replication/pull_metabric.py` (API pull, gene-symbol ->
Entrez, within-cohort z-score, locked score, directional AUC + bootstrap
cluster CI, sign concordance). Prereg lock @e54d8bc precedes ALL of this.
