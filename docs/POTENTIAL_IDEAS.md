# GenoPred V2 — Potential Ideas

This is a living backlog of scientifically plausible additions. Inclusion here does not mean that
a dataset is approved for model training. Every idea must pass provenance, biological relevance,
sample-independence, genome-build and leakage checks before promotion into the feature schema.

Last reviewed: **2026-08-26**

## Evidence tiers

1. **Primary driver evidence:** PAAD somatic mutation and copy-number recurrence in independent
   patient cohorts.
2. **Tumour-context evidence:** PAAD expression, chromatin and regulatory features.
3. **Functional evidence:** pancreatic-lineage dependency or perturbation experiments.
4. **Metastatic-state evidence:** expression associated with CTCs, metastases or plasticity.
5. **Cross-cancer or cross-species evidence:** validation only unless a specific transfer-learning
   experiment demonstrates benefit without negative transfer.

## Ideas from Mohit Kumar Jolly's work

### GSE114704 — human PDAC PDX single-cell RNA-seq

- **Priority:** High supporting dataset.
- **Status:** Candidate; not downloaded yet.
- **Source:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114704>
- **Contents:** 70 quality-controlled human cells from the Panc-265 metastatic PDX: 37 primary
  tumour cells, 23 liver-metastasis cells and 10 circulating tumour cells.
- **Potential features:** primary-to-metastasis differential expression, primary-to-CTC differential
  expression, metastatic-state specificity, fraction of CTC/metastatic cells expressing a gene, and
  agreement with TCGA-PAAD bulk expression.
- **Proposed role:** Orthogonal metastatic-expression validation and candidate prioritization.
- **Hard restriction:** Do not turn differential expression into driver/non-driver labels. This is a
  single PDX and does not establish population-level mutation causality.
- **Relevant paper:** <https://pubmed.ncbi.nlm.nih.gov/32499301/>

### GSE51372 — pancreatic circulating tumour cells

- **Priority:** Medium/optional.
- **Status:** Candidate; not downloaded yet.
- **Source:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51372>
- **Contents:** Single-cell RNA-seq of pancreatic CTCs, matched primary tumour cells and controls
  from a mouse pancreatic-cancer model.
- **Potential use:** Test whether metastatic programs found in GSE114704 reproduce independently.
- **Hard restriction:** Validation only. Mouse-to-human orthologue mapping and older sequencing
  technology make this unsuitable as a primary training cohort.

### EMT and epithelial–mesenchymal plasticity scores

- **Priority:** High method addition; no extra primary cohort required.
- **Status:** Proposed for TCGA-PAAD RNA-seq.
- **Method source:**
  <https://www.frontiersin.org/journals/bioengineering-and-biotechnology/articles/10.3389/fbioe.2020.00220/full>
- **Potential features:** 76GS, KS and multinomial-logistic-regression EMT scores; gene-to-EMT
  correlation; hybrid epithelial/mesenchymal association; plasticity-stratified expression.
- **Proposed role:** Tumour-state covariates and biological interpretation.
- **Hard restriction:** EMT association is not proof of driver status. Scores must be calculated only
  within training folds whenever they influence feature selection or normalization.

## Recommended evaluation experiment

Build the PAAD-only baseline first. Then add the Jolly-derived evidence in controlled ablations:

1. TCGA-PAAD mutation/CNV baseline.
2. Add TCGA-PAAD expression.
3. Add EMT/plasticity-derived features.
4. Add GSE114704 metastatic-expression evidence.
5. Use GSE51372 only as an external cross-species consistency check.

Report the change in held-out performance, calibration and ranked-gene stability at every step. A
feature group that does not improve untouched PAAD validation should not be retained simply because
it is biologically interesting.

## Other future ideas

- Add ICGC APGI-AU/PACA-CA as independent genomic validation when access is resolved.
- Add CPTAC PDAC protein abundance as orthogonal evidence rather than a substitute for genomics.
- Add DepMap pancreatic-lineage CRISPR dependency as functional validation, not driver labels.
- Evaluate tumour-purity and stromal-expression confounding before using bulk RNA features.
- Compare PAAD-from-scratch against pan-cancer pretraining and retain transfer only if it improves
  patient-disjoint PAAD evaluation.
