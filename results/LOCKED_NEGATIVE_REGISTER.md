# Locked negative register, October 1, 2026

These are existing measured results, not new evaluations.

- METABRIC replication: n=1979, events803,15 genes, directional AUC0.5117236807549919, one-sided95% lower bound0.47780695557920955; primary test false.
- SCAN-B: n=3409, events353,15 genes, AUC0.5075289589605921, lower bound0.4790065356391322; secondary test false. This is explicitly an overall-survival arm, not relapse.
- Two-cohort discovery: zero genes meet locked joint stability and sign agreement. No panel rescued on replication outcomes.
- Permutation arm:100 permutations, observed stable genes268, null median426,max630; panel-size p0.9801980198019802. Stable-panel count is not evidence above this null.

Sources: results/replication/metabric_primary.json,scanb_secondary.json,twocohort_panel.json and results/permutation_arm.json.

A different pretreatment therapy-response question is a candidate for research, not an approved protocol. It requires genuine response labels, clinical baseline, external cohort and matched prior art. Existing relapse/OS outcomes cannot be renamed treatment benefit. No thresholds were loosened, no replication/test reselection was performed and no new training was started for this register. Paper revision and steering consultation remain open.

## Matched-null scope correction
Source inspection October1 found the completed permutation implementation uses probe-level features and de_genes filtering, whereas the observed268 count is symbol-level, with max-|t| probe collapse and tstats_chunked filtering. The permutation folds are formed from original labels before label permutation; observed folds use actual labels. The numeric p0.9801980198 remains the recorded implementation output, but is not a validated identical-pipeline significance test for the268 symbol-level count. Do not use it to certify a general null-calibrated tool or overwrite it silently. A separately preregistered matched positive/null-control study would be needed; none has started.

The observed selector describes L1 fits in its docstring, but code omits penalty and installed LogisticRegression default is L2. The methods description must report the effective L2 implementation, and any conflict with the preregistration remains explicit. No completed result was recomputed.

## Protocol deviation, not faithful L1 falsification
PREREG_TWOCOHORT.md section2 step3 explicitly requires L1; completed code used defaultL2. Panel0 is the measured L2 result only. The intended L1 experiment is unverified. Do not label this failed-as-designed under the L1 preregistration. No preregistration was amended retroactively, and no new L1 execution was started.
