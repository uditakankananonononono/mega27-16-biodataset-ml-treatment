# External-tool provenance audit, October 7, 2026

## Scope and verdict

Static inspection of the current tracked source and historical JSON only. No models, API inventory, or locked experiments were rerun. The 40-tool scientific-use gate remains uncertified. This is an audit milestone, not a novelty result or project completion.

The registry and output each contain 48 names in matching order. The historical output reports 45 ok and three failed. These are recorded execution statuses, not 45 independently verified, valid cohort analyses. Historical output is preserved unchanged.

## Material issues

- `bctpy` records `mean_clustering: Infinity`. This is not a finite scientific estimate, and Infinity is not strict JSON. The source does not clear the correlation matrix diagonal before passing it to the graph functions. That is a candidate input-validity problem, not an established explanation for the output.
- `scikit-optimize` selects its top 200 probes using all labels and normalizes all samples before BayesSearchCV. Its recorded best CV AUC 0.7658 is a model-selection output with upstream leakage, not independent predictive performance. No replacement score was calculated.
- `statsmodels` selects a best probe from the same cohort before its logit inference. The recorded top-probe p-value does not account for that selection. Its 2000 FDR rejections also need raw p-values and input-scale checks before interpretation.
- The shared context chooses top probes by whole-cohort label association. Descriptive plots, SHAP, enrichment, and feature-selection summaries on that context are exploratory, not held-out discoveries. `_cv_auc` and the SMOTE function separately select on training folds, so the shared-context issue should not be generalized to those functions.
- `mygene` reports a cached context result. The builder swallows annotation exceptions and can replace symbols with a hard-coded list. It stores no raw query-response mapping, fallback flag, or mapping version. The reported 129 mappings are compatible with multiple hits for 100 requested probes, not proof of 129 independent genes or successful one-to-one mapping. The function's text says top 25, while the builder requests top 100.
- `CrossRef` takes the first fuzzy search hit without identity checks. Its recorded DOI `10.1016/s1043-321x(06)80430-3` is not the primary Wang 2005 DOI. PubMed PMID 15721472 and the Erasmus institutional record identify `10.1016/S0140-6736(05)17947-1`. The historical wrong hit is not silently overwritten.
- KEGG records zero pathways without checking HTTP status. Reactome can turn any dictionary response into an empty list via `.get('pathways', [])`. Neither recorded zero establishes biological absence. QuickGO records five null term names despite status ok. Raw responses and schema checks are missing.
- gseapy only lists libraries. Sympy checks a fixed formula. WordNet gives a definition. PubMed, Europe PMC, NCBI eutils, and CrossRef are search/provenance checks. They are not additional cohort experiments. NLTK counts symbol tokens, not clinical language or biological validation.
- Several entries reuse the same provider or underlying data. A provider/client pair and multiple searches do not automatically establish distinct scientific tools under the program's counting rule. No arbitrary revised certified count is substituted here.
- The runner's success criterion is that a function returned without an exception. It records no code revision, data checksum, package versions, raw response snapshots, semantic checks, or fold-level prediction artifacts for this inventory. Existing source and outputs establish plausibility, not complete execution provenance.

## Primary-source checks

Fetched October 7, 2026:

- https://pubmed.ncbi.nlm.nih.gov/15721472/ : Wang et al., Lancet 2005, DOI 10.1016/S0140-6736(05)17947-1; 286 untreated lymph-node-negative patients.
- https://repub.eur.nl/pub/57509 : institutional record independently identifies the same primary DOI and article.
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034 : series record connects GSE2034 to citation 15721472. The cohort's records do not become independent datasets.

## What remains open

Agree and apply a distinct-tool counting rubric; validate response identity and schema; recover raw responses if they exist; capture data/code/package provenance; inspect original fold outputs and preprocessing for each predictive claim; audit native paper claims and regenerate and visually inspect it after corrections. No threshold change, historical protocol repair, test reselection, or training restart is authorized by this audit.
