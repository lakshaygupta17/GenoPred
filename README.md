# GenoPred

GenoPred V2 is a reproducible research pipeline for ranking pancreatic ductal adenocarcinoma
(PDAC/PAAD) driver-gene candidates from multi-omic and graph data.

The project replaces the stateful GenoPred notebooks with a provenance-aware pipeline. Its first
goal is not to reproduce a headline score; it is to construct a biologically coherent dataset and
an evaluation design that can support trustworthy conclusions.

## Current stage

**Stage 6 complete; initial Stage 8 functional assessment completed (weak external support; Stage 7 pending)**

Last updated: **2026-11-04**

| Stage | Status | Deliverable |
|---|---|---|
| 0. Audit V1 | Complete | Identified label leakage, fixed-test reuse, coordinate and notebook defects |
| 1. Foundation | Complete | Package structure, source manifest, feature schema and validation gates |
| 2. Acquisition | Complete | Download manifests, checksums and immutable raw files |
| 3. Harmonization | Complete | GRCh38 coordinates, gene identity map and assay QC |
| 4. Labels | Complete | Leakage-controlled PAAD positive/unlabelled evidence table |
| 5. Dataset | Complete | Model-ready gene feature table and graph |
| 6. Baselines | Complete | Non-graph and PAAD-from-scratch reference models |
| 7. Transfer learning | Not started | Pan-cancer and MCF-7 transfer with negative-transfer controls |
| 8. External validation | Initial functional assessment complete; broader validation pending | Frozen shortlist tested against Sanger Project Score; weak associations and no top-100 enrichment |

## Built so far

- A `src/` Python package with a minimal mandatory runtime dependency set.
- A machine-readable data-source manifest covering GEO, TCGA, ICGC, CPTAC, DepMap and OncoKB.
- A proposed feature schema with explicit source, transformation and missing-data policy.
- A hard validation gate preventing silent hg19/GRCh37 and GRCh38 mixing.
- Feature-matrix validation for duplicate genes, wrong dimensions, NaN and infinity.
- A small CLI for project status and configuration validation.
- Unit tests for coordinate and configuration safeguards.
- Data-directory rules that keep raw data immutable and out of Git.
- Reproducible GEO and GDC inventory generation with explicit core, supporting and validation tiers.
- A resumable parallel downloader with disk-space checks and SHA-256/MD5 verification.
- Verified core snapshots for GSE149103, GSE124229 and TCGA-PAAD.
- Stage 3 file-schema and sample-identity inventories covering every core file.
- Exact GDC file-to-case, sample-type and aliquot metadata with a checksummed source snapshot.
- A living potential-ideas backlog with scientific evidence tiers and leakage restrictions.
- Checksum-recorded UCSC liftOver and GENCODE v36 reference assets.
- A GRCh38 patient ATAC atlas with interval- and signal-retention QC.
- A canonical Ensembl/HGNC gene identity table and an unambiguous tumour-only TCGA cohort.
- A mentor-ready technical brief documenting V1 defects, V2 corrections, transfer-learning risks,
  completed evidence and the remaining Stage 4-8 work at the time of the audit.
- A plain-language glossary and notation guide covering the biology, datasets, file formats,
  genome coordinates, graph tensors, modelling terms, metrics and transfer-learning vocabulary.
- A checksum-pinned IntOGen label-evidence snapshot and a deterministic 19,962-gene PAAD label
  table with 17 externally replicated positives and 19,945 unlabelled genes.
- Explicit positive-unlabelled modelling policy; no false reliable negatives manufactured from
  absence of catalog evidence.
- A frozen 19,962-by-20 finite model matrix containing 10 PAAD biological features and 10 explicit
  missing-value flags, aligned exactly to the Stage 4 gene order.
- A STRING v12.0 high-confidence physical graph with 170,342 directed edges and auditable mappings.
- A fixed-split nnPU logistic baseline and a STRING random-walk baseline, including the complete
  validation-selection trail and exploratory evaluation on a previously examined test split.

The first baselines have now been trained. Raw sources remain immutable, preprocessing was fitted
inside the approved training partitions, and unlabelled genes were not asserted to be true
negatives.

## Acquired core data

| Source | Verified files | Size | Contents |
|---|---:|---:|---|
| GSE149103 | 18/18 | 5.63 GiB | PANC-1/HPNE ATAC, H3K4me3, H3K27ac, CTCF and Hi-C |
| GSE124229 | 1/1 | 11.7 MiB | Patient ATAC atlas count matrix |
| TCGA-PAAD | 542/542 | 1.30 GiB | STAR RNA counts, masked somatic mutations and ASCAT3 gene-level CNV |

The 80 GSE124229 patient bigWigs (50.9 GiB) are retained in the validation inventory but were
intentionally deferred. The compact atlas count matrix provides the correct first-pass input without
spending roughly 51 GiB before its sample annotations and genome build pass QC.

## Stage 3 QC snapshot

- 561/561 acquired files inspected with zero schema failures.
- 542 TCGA files are GRCh38-compatible; all 19 hg19 GEO containers now have an enforced
  harmonization route and cannot be directly mixed with GRCh38.
- TCGA resolves to 185 cases and 375 distinct sample barcodes across file associations.
- TCGA schemas are internally uniform: 183 STAR-count files, 179 ASCAT3 CNV files and 180 MAFs.
- GSE124229 contains 80 assay columns in one coordinate-by-count matrix; 119,299/119,466
  intervals (99.860%) and 98.442% of signal mass were retained after liftover to GRCh38.
- All 16 GSE149103 bigWigs are valid hg19 tracks with consistent shared chromosome sizes; both
  Hi-C archives are valid and contain 24 chromosome files apiece plus macOS metadata entries.
- The GENCODE v36 map contains 60,660 assembly gene records and resolves 100% of identifiers in
  representative TCGA RNA (60,660) and CNV (60,623) files.
- The canonical tumour-only cohort contains 165 cases with RNA, CNV and mutation data.
- Detailed outputs are in `manifests/stage3_schema_report.json`,
  `manifests/stage3_file_inventory.csv` and `manifests/stage3_sample_inventory.csv`.
- Future dataset and modelling candidates are tracked in `docs/POTENTIAL_IDEAS.md`.
- Stage 3 transformation and cohort rules are documented in `docs/STAGE3_QC.md`.

## Stage 4 label snapshot

- Gene universe: 19,962 GENCODE v36 protein-coding genes.
- Positive: 17 genes replicated across at least two qualifying external cohort families.
- Unlabelled: 19,945 genes; these are unknown, not presumed non-drivers.
- Reliable negative: 0; absence from a catalog was not converted into a biological negative.
- Qualifying families: CPTAC, ICGC PACA-AU and ICGC PACA-CA.
- TCGA-PAAD evidence is audit-only because TCGA mutation prevalence will be a model feature.
- PCAWG evidence is supporting-only because its pancreatic cohort overlaps ICGC donors.
- Since CPTAC and ICGC contribute to label evidence, they are not untouched Stage 8 validation
  sources for this label version; later validation must be orthogonal or genuinely held out.
- Full rules and limitations are documented in `docs/STAGE4_LABELS.md`.

## Stage 5 dataset snapshot

- Matrix: 19,962 genes by 20 columns: 10 biological features plus 10 missing flags.
- Modalities: TCGA-PAAD expression, damaging mutation and CNV; patient PDAC ATAC; PANC-1 versus
  HPNE ATAC, CTCF, H3K4me3 and H3K27ac.
- Current dataset is `stage5-core-v2`. Mutation parsing is corrected, ATAC replicates are combined
  within 40 patients, and distorted promoter mappings are rejected.
- CNV uses a per-case autosomal modal baseline. Sex-chromosome CNV is missing pending a validated
  sex-aware calculation; 953 genes have unavailable CNV values. Missingness remains explicit.
- Graph: STRING v12.0 physical links at score >=700, yielding 85,171 undirected interactions
  stored as 170,342 directed edges across 10,599 connected genes.
- `LRP1B` is the only positive isolated by that strict graph; Stage 6 must retain feature-only
  prediction for isolated genes and report sensitivity rather than weakening the threshold.
- Hi-C features are explicitly deferred; DepMap CRISPR remains validation-only. Neither is replaced
  with fabricated values.
- Full methods, limitations and outputs are documented in `docs/STAGE5_DATASET.md` and
  `manifests/stage5_dataset_report.json`.

## Stage 6 baseline snapshot

- Fixed split: 11/3/3 known positives for training/validation/test, with 11,967/3,989/3,989
  unlabelled comparators.
- Feature baseline: linear non-negative positive-unlabelled logistic model trained from scratch.
- Graph baseline: weighted STRING random walk with restart and explicit self-loops.
- Corrected exploratory feature result: observed-label AUROC 0.9653 and average precision 0.22465;
  two of three known positives were recovered in the test subset's top 100.
- Graph diffusion: AUROC 0.9658, average precision 0.08844. It reduces average precision despite
  a slightly higher AUROC. These results supersede the defective v1 experiment.
- These are provisional known-label retrieval metrics with only three test positives—not clinical
  sensitivity, specificity or proof that top-ranked unlabelled genes are drivers.
- Full design and interpretation are in `docs/STAGE6_BASELINES.md` and
  `manifests/stage6_baseline_report.json`.
- All six audit findings and their corrections are recorded in `docs/AUDIT_CORRECTIONS.md`.
- `data/processed/stage6/model.npz` contains a complete reloadable model, including its intercept.

## Initial data plan

| Purpose | Preferred source |
|---|---|
| Coherent PANC-1/HPNE epigenome and Hi-C | GEO GSE149103 |
| Patient-derived ATAC-seq | GEO GSE124229 |
| PAAD mutation, CNV, RNA and methylation | TCGA-PAAD |
| Label-support cohorts | ICGC/PCAWG and CPTAC through the pinned IntOGen snapshot |
| Independent Hi-C validation | GEO GSE185069 |
| Orthogonal functional validation | DepMap pancreas lineage |
| Held-out curated evidence | OncoKB or CIViC pancreatic evidence acquired after model lock |

The canonical coordinate system is **GRCh38**. GSE149103 and GSE124229 processed tracks are
reported as hg19 and must be harmonized before feature extraction.

## Frozen Stage 5 feature groups

- PANC-1 versus HPNE differential ATAC, CTCF, H3K4me3 and H3K27ac signals.
- Patient-level PDAC ATAC mean and variability.
- Consequence-filtered somatic-mutation prevalence.
- Separate amplification and deletion prevalence.
- Distance-normalized local and distal promoter contacts (explicitly deferred from `stage5-core-v1`).
- Tumor RNA expression.
- DepMap pancreatic dependency as validation evidence, not automatic driver labels.

The implemented definitions are frozen in `configs/features.json`; deferred features cannot enter a
later dataset version until their source files and genome builds pass the same QC gates.

## Running the foundation checks

From this directory with Python 3.11 or newer:

```powershell
$env:PYTHONPATH = "src"
python -m genopred_v2 status
python -m genopred_v2 validate-config
python -m genopred_v2 acquire-stage3-metadata
python -m genopred_v2 stage3-qc
python -m genopred_v2 complete-stage3
python -m genopred_v2 complete-stage4
python -m genopred_v2 complete-stage5
python -m genopred_v2 complete-stage6
python -m pytest -q
```

For editable installation during development:

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -e ".[dev]"
genopred-v2 validate-config
```

## Repository layout

```text
GenoPred V2/
├── configs/             Source and feature manifests
├── data/                Raw, interim and processed data (not committed)
├── docs/                Scientific and engineering decisions
├── manifests/           Versioned inventories, checksums and acquisition reports
├── scripts/             Reproducible acquisition entry points
├── src/genopred_v2/     Reusable pipeline code
├── tests/               Automated safeguards
├── pyproject.toml
└── README.md
```

## Stage log

### 2026-08-21 — Stage 1 started

- Created the clean V2 project rather than modifying the original notebooks.
- Selected GRCh38 as the canonical build.
- Registered the first eight data sources and thirteen proposed features.
- Added an explicit guard for the hg19/GRCh38 mismatch found in V1.
- Defined the three-state label policy and untouched-test evaluation rule.

### 2026-08-21 — Stage 2 completed

- Snapshotted the official GEO SOFT records and the complete open TCGA-PAAD GDC file inventory.
- Classified every discovered file into core, supporting, validation or not-selected tiers.
- Downloaded and checksum-verified 561 core assay files (6.94 GiB total).
- Selected harmonized TCGA-PAAD STAR counts, masked mutation calls and ASCAT3 gene-level CNV.
- Kept the 50.9 GiB GSE124229 bigWig tier reproducibly addressable but deferred pending QC.
- Added resumable downloads, transient-failure recovery, disk safety checks and per-source reports.

### 2026-08-26 — Stage 3 started

- Added `docs/POTENTIAL_IDEAS.md` and recorded GSE114704, GSE51372 and EMT scoring as
  controlled candidate additions inspired by Mohit Kumar Jolly's work.
- Inspected all 561 core files and confirmed a single consistent schema for each TCGA modality.
- Snapshotted exact GDC file, case, sample-type and associated-aliquot metadata for all 542 files.
- Mapped the TCGA cohort to 185 cases and 375 sample barcodes without modifying raw inputs.
- Confirmed the 80-column GSE124229 ATAC matrix structure and all GSE149103 file containers.
- Enforced the build gate: 19 hg19 GEO files cannot enter feature extraction before liftover.
- Expanded the automated safeguard suite to 19 passing tests.

### 2026-08-26 — Stage 3 completed

- Acquired five checksum-recorded UCSC/GENCODE reference assets, including the repository-verified
  hg19-to-hg38 chain.
- Validated 16 bigWig chromosome trees and fully inventoried two Hi-C archives without modifying
  raw data.
- Materialized the GSE124229 atlas in GRCh38 with 99.860% interval and 98.442% signal retention.
- Built a 60,660-record GENCODE v36 Ensembl/HGNC identity map with explicit pseudoautosomal and
  duplicate-symbol collision policy.
- Verified 100% gene-ID coverage for representative TCGA RNA and CNV inputs.
- Selected 178 primary-tumour RNA cases and a final 165-case three-modality modelling cohort with
  no duplicate primary-tumour RNA assignments.
- Passed all 22 automated tests and every Stage 3 completion gate.
- Added `docs/GenoPred_V2_Mentor_Brief.docx` as the reviewed project-audit handoff document.

### 2026-08-27 — Project glossary added

- Added `docs/GenoPred_V2_Glossary_and_Notation_Guide.docx` as a 21-page study and mentor-reference
  guide covering project-specific biological, data-engineering and machine-learning terminology.
- Documented the original 15-feature order, graph tensor shapes, three-state label policy,
  evaluation metrics, transfer-learning notation and common terminology confusions.

### 2026-08-27 — Stage 4 completed

- Acquired and checksum-pinned the CC0 IntOGen driver and cohort archives advertised as release
  2024.09.20.
- Corrected the Stage 3 GENCODE parser so `hgnc_id` contains actual `HGNC:` identifiers rather than
  a duplicate gene-symbol field; regenerated the identity map and coverage report.
- Built a 19,962-gene protein-coding label universe with 17 positives supported by at least two
  external families among CPTAC, ICGC PACA-AU and ICGC PACA-CA.
- Retained PCAWG as overlapping support and TCGA-PAAD as audit-only evidence; neither can inflate
  the count of independent families used to assign labels.
- Left 19,945 genes unlabelled and created zero reliable negatives, establishing positive-unlabelled
  learning as the approved downstream regime.
- Added `docs/STAGE4_LABELS.md`, long-form evidence provenance, deterministic compressed outputs,
  Stage 4 CLI commands and six new tests.

### 2026-09-03 — Stage 5 completed

- Froze `stage5-core-v1` with 10 biological PAAD features and 10 missingness indicators across all
  19,962 Stage 4 protein-coding genes.
- Aggregated a fixed 165-case TCGA-PAAD cohort for expression, damaging mutation and absolute CNV
  prevalence without reusing TCGA evidence as a label source.
- Extracted promoter features from 80 patient PDAC ATAC profiles and differential PANC-1/HPNE
  ATAC, CTCF, H3K4me3 and H3K27ac tracks with explicit hg19/GRCh38 handling.
- Acquired and checksum-pinned the STRING v12.0 physical network under CC BY 4.0, creating 170,342
  directed high-confidence edges.
- Preserved every missing value through a companion flag, deferred scaling to Stage 6 folds and
  retained DepMap as validation-only evidence.
- Documented the isolated-node limitation, deferred Hi-C rather than synthesizing it, and passed
  all 31 automated checks.

### 2026-09-17 — Stage 6 completed

- Added an exact deterministic gene split that preserves 11/3/3 known positives across
  training/validation/test and never consults features while splitting.
- Implemented a dependency-light linear nnPU baseline so unlabelled genes are not silently treated
  as verified negatives.
- Refit both missing-value medians and scaling only on allowed training partitions, then froze
  class-prior and regularization choices before opening the test labels.
- Added weighted STRING random-walk diffusion with self-loops and an explicit transductive scope.
- Found that the feature-only baseline outperformed graph diffusion on the locked test; retained the
  weaker graph result and documented the `LRP1B` isolated-node sensitivity.
- Added deterministic score/model artifacts, a complete selection audit, CLI integration and four
  Stage 6 safeguards, bringing the suite to 35 passing tests.

### 2026-09-17 — Audit corrections and v2 rebuild

- Archived the defective v1 code, documentation, manifests and Stage 5/6 outputs under
  `audit_archive/pre_correction_v1`; earlier stage-log completion claims are historical.
- Fixed tab-separated MAF parsing, required-column checks and tumour/case validation. Mutation
  prevalence is now nonzero for 3,525 genes; KRAS is 94/165 and TP53 is 90/165.
- Added promoter-length, orientation and interior liftover checks, and collapsed 80 ATAC libraries
  to 40 patient profiles before computing cohort summaries.
- Replaced fixed CNV thresholds with per-case autosomal modal baselines; deferred sex-chromosome
  CNV. Chromosome X now contributes 5 of the feature model's top 100 genes, versus 68 previously.
- Corrected tie handling and persisted the intercept and preprocessing in a complete model bundle.
- Passed 41 tests and reproduced feature and graph predictions exactly from saved parameters.
- Rebuilt Stages 5 and 6 as v2; retained the same split and explicitly labelled results exploratory.

### 2026-09-18 — Plain-language project guide

- Added `docs/GenoPred_V2_Stages_1_to_6_Explained_Simply.docx`, a 16-page visual guide for readers
  without a biology background. It explains basic biology, all six stages, a worked gene example,
  corrected v2 results, audit repairs, limitations and a short glossary. Fourteen diagrams, charts
  and compact tables now provide a visual explanation alongside the plain-language text.
- Checked the rendered document page by page. This documentation update does not change models
  or advance the project beyond Stage 6.

### 2026-11-04 — Stage 6 PU tree experiments completed

- Added shallow PU tree bagging and PU Random Forest bagging with validation-only selection
  across three fixed seeds, preserving the corrected data and original gene split.
- Observed-label test average precision: 0.41727 for tree bagging, 0.37870 for forest bagging,
  and 0.22465 for existing linear nnPU. Both tree methods recover 2 of 3 positives in the top
  500, versus 3 of 3 for nnPU; no universal winner is claimed.
- Saved independent model bundles, gene scores and a complete comparison report. Both models
  reproduce saved predictions exactly. All 44 automated tests pass.
- See `docs/STAGE6_TREES.md` and `manifests/stage6_tree_report.json`. These are exploratory
  results on the previously examined test split. Stage 7 remains pending.

### 2026-11-04 — Broad Stage 6 model benchmark completed

- Completed all 55 registered experiments, a validation-selected rank ensemble and four earlier
  references: 60 comparison rows, with no failed attempts. The screen covers 20 PU-bagged learner
  families, alternative learning assumptions, graph methods, all ten individual features, seven
  feature ablations, four preprocessing variants and three controls.
- Compared 95 registered candidate settings using three fixed seeds. Choices and ensemble members
  were saved before test evaluation. Models, preprocessing, scores, source hashes, environment
  versions, warnings and seed sensitivity are preserved under `data/processed/stage6_benchmark_v1`.
- The shallow PU decision-tree screen has the highest observed test AP (0.66721), followed by
  gradient boosting (0.55641) and LightGBM (0.55627). Each recovers 2 of 3 positives in the top 500;
  the original nnPU recovers all 3. These results do not establish a validated winning model.
- All saved-model round trips have zero prediction difference. All 50 tests pass and dependency
  checks pass. QDA produced recorded covariance-rank warnings despite regularization.
- Added a complete practical catalog and explicit unrun research backlog in
  `docs/STAGE6_EXPERIMENT_CATALOG.md`, the 60-row comparison in
  `docs/STAGE6_BENCHMARK_RESULTS.md`, and `docs/figures/stage6_model_comparison.png`.
- Run with `.venv\Scripts\python.exe -m genopred_v2.stage6_benchmark .`. This is a broad,
  bounded exploratory screen on the consumed split, not an exhaustive search of every possible
  algorithm. Stage 7 and independent external validation remain pending.

### Next milestone after correction

Stage 7: pre-register controlled transfer-learning and negative-transfer experiments. The consumed
Stage 6 test cannot be used for Stage 7 tuning; later performance claims need a newly locked or
external validation design.

### 2026-11-04 — Initial external functional assessment completed

- Froze nnPU, PU decision tree, PU gradient boosting, PU LightGBM and a mutation-prevalence
  control before evaluating external functional measurements. Original predictions remain intact.
- Downloaded and checksum-pinned the official Sanger Project Score combined Sanger v2/Broad
  21Q2 archive (June 2025), following HTTP 403 responses from the proposed Broad archive.
- Evaluated 15,767 eligible unlabelled, non-common-essential genes across 33 pancreatic
  adenocarcinoma donors, with 856 other-cancer donors for a specificity comparison. Excluded
  PANC-1 and its donor-related KP-1N; collapsed models by donor and checked cohort disjointness.
- Primary rank correlations: nnPU 0.0229, decision tree 0.0670, gradient boosting 0.0504,
  LightGBM 0.0695. These are small effects; specificity is close to zero. No shortlisted model
  has top-100 enrichment. Decision-tree top-500 includes 1,684 genes because of ties.
- This measures cancer-cell survival dependency, not all driver mechanisms or clinical accuracy.
  No winning driver-gene model is established and no model was tuned on external outcomes.
- Added `docs/EXTERNAL_FUNCTIONAL_VALIDATION.md`, a comparison chart, immutable protocol,
  schema rules, cohort audit, gene-level evidence, metrics, source hashes and output manifest.
- All 55 tests passed; saved-output checks confirm exact frozen scores and eligibility/exclusions.
- Stage 7 remains pending. This is an initial Stage 8 functional assessment only. These outcomes
  are now consumed and must not be called an untouched final test for future model changes.

### 2026-11-04 — MSK driver-evidence audit and provisional benchmark

- Downloaded and pinned the MSK 2,336-patient supplementary tables, common IMPACT341 panel and HGNC symbol mapping.
- Identified 80 new recurrent sequence-evidence genes among 326 eligible panel genes after excluding original positives. These are provisional evidence labels, not 80 independently confirmed PDAC drivers.
- Evaluated frozen models without retraining: observed-label AP 0.4179 nnPU, 0.3490 decision tree, 0.4365 gradient boosting and 0.4056 LightGBM. nnPU top-100 recall is 47.5%.
- Saved the full audit, exclusions, source hashes, frozen protocol and results in `docs/MSK_DRIVER_EVIDENCE_AUDIT.md` and `data/processed/msk_driver_audit_v1/`.
- Cross-study donor overlap and variant-level annotation remain unresolved. Mouse screen preparation is blocked by unavailable supplement downloads. No fully independent gold-standard claim is made.
- Original labels and models are unchanged; Stage 7 remains pending.

### 2026-11-04 — Final exploratory evaluation snapshot

- Completed 5 repeats x 3 outer folds with 3-fold inner model selection for nnPU, shallow PU trees, gradient boosting and LightGBM; added mutation baseline, feature ablations and shuffled-label controls.
- Mean observed-label AP: nnPU 0.5517, mutation baseline 0.5025, LightGBM 0.4490, gradient boosting 0.3160, shallow tree 0.1365. Mutation features dominate; the historical single-split tree result is not the final headline.
- Audited all 80 provisional MSK evidence genes: 917/924 supporting variant rows match public records; explicit per-variant driver fields and cross-cohort donor verification remain unavailable.
- Added a separate 30-donor autopsy evidence assessment: author-selected sequence records yield 15 new genes, of which 9 overlap MSK. Documented source-count discrepancies and quarantined the overly broad LiFD extraction.
- All 58 tests pass, 60 nested model round trips are exact, and all repeated test partitions cover each gene once.
- Final results, figures and report outline: `docs/FINAL_RESULTS.md` and `docs/REPORT_WRITING_PLAN.md`. Original labels/models are unchanged. Stage 7 is not implemented; strict independent and clinical validation remain unresolved; mouse evaluation remains unavailable.

### 2026-11-04 — Separately tuned features and final sensitivity experiments

- Completed independently tuned mutation-only, mutation+CNV, mutation+expression and full-feature comparisons for nnPU, gradient boosting and LightGBM.
- nnPU mutation-only AP 0.5025; full nnPU AP 0.5517. Highest mean AP in this comparison: nnPU / Full features (0.5517).
- Completed 17 leave-one-positive-out tests per four procedures, held-out candidate-stability analysis and fixed-prior sensitivity.
- All 64 tests pass; 248 saved model round trips are exact. Original labels/models remain unchanged.
- See `docs/FINAL_EXTENSION_RESULTS.md`; this is exploratory sensitivity analysis, not a new independent external test.

## Report and presentation

- [B.Tech project report](report/GenoPred_BTP_Report.pdf)
- [Editable LaTeX source](report/latex/main.tex)
- [10-minute presentation](presentation/GenoPred_10_Minute_Presentation.pptx)

The Python package and CLI retain their existing `genopred_v2` and `genopred-v2` names for compatibility. Raw and processed datasets, trained models, local environments and historical snapshots are excluded.
