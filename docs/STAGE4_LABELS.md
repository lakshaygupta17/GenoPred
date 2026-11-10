# Stage 4 label design and QC

## Outcome

Stage 4 constructs a versioned gene-level evidence table for pancreatic adenocarcinoma without
using the TCGA-PAAD feature cohort to define its own answers. The completed table contains 19,962
GENCODE v36 protein-coding genes:

- 17 `positive` genes supported by at least two qualifying external cohort families.
- 19,945 `unlabelled` genes.
- 0 `reliable_negative` genes.

The zero-negative result is deliberate. Absence from a driver catalog is not proof that a gene is a
true non-driver. Downstream modelling must therefore use a positive-unlabelled learning or ranking
formulation unless a separately reviewed negative-evidence rule is added.

## Evidence source

The reproducible source is the IntOGen driver catalog release labelled 2024.09.20 on the official
download page:

- <https://www.intogen.org/download>
- Driver and cohort archives are stored as immutable raw inputs.
- Both archives are pinned by SHA-256 and verified as valid ZIP files.
- The bundled license is CC0 1.0 Universal.

The archive's internal folder is dated 2024-06-18 and its README still says 20230531. The project
therefore identifies the snapshot by the official download filename, advertised release date, exact
SHA-256 and archive contents rather than trusting the stale internal README label.

OncoKB is not used in the automatic Stage 4 build because its current API requires a registered
token and license. It remains an optional future evidence adapter, not an undocumented manual
substitute.

## Positive-label rule

Only IntOGen rows with `CANCER_TYPE = PAAD` and `IS_DRIVER = True` are considered. Pancreatic
neuroendocrine tumours and non-specific metastatic `PANCREAS` cohorts are excluded.

A gene is positive when it appears in at least two of these independent external families:

1. CPTAC pancreatic adenocarcinoma.
2. ICGC PACA-AU.
3. ICGC PACA-CA.

WXS and WGS analyses from the same program/country count as one family. The PCAWG pancreatic
adenocarcinoma cohort is retained as supporting evidence but is not counted as independent because
it overlaps ICGC donors. TCGA_WXS_PAAD is also retained in the evidence ledger but is never allowed
to qualify a label because TCGA-PAAD mutation prevalence is a planned model feature.

The 17 positive genes are:

`ACVR1B`, `ACVR2A`, `ARID1A`, `ATM`, `CDKN2A`, `FAT1`, `KDM6A`, `KMT2C`, `KMT2D`, `KRAS`,
`LRP1B`, `RBM10`, `RNF43`, `SMAD4`, `SMARCA4`, `TGFBR2`, and `TP53`.

These are high-confidence computational labels under the stated rule. They are not a claim that the
list is complete, and they are not clinical assertions.

## Identifier policy

IntOGen HGNC symbols are joined only when they map to exactly one GENCODE v36 protein-coding gene.
All 146 PAAD evidence rows resolved uniquely. The canonical output key is the unversioned Ensembl
gene ID, while the approved symbol and HGNC identifier are retained for interpretation.

During this work, the Stage 3 GENCODE parser was corrected: the source metadata contains transcript,
gene-symbol and HGNC-ID columns, and the original implementation had read the symbol column into the
`hgnc_id` field. The identity map and coverage report were regenerated after the fix.

## Outputs

- `configs/labels.json` — frozen source hashes, cohort-family mapping and label rules.
- `data/processed/labels/paad_gene_labels.tsv.gz` — one row per protein-coding gene.
- `data/processed/labels/paad_label_evidence.tsv.gz` — long-form evidence ledger.
- `manifests/stage4_source_report.json` — source archive provenance and checksums.
- `manifests/stage4_label_report.json` — label counts, positive list and QC decisions.
- `manifests/stage4_completion_report.json` — Stage 4 completion gate.

Compressed outputs use a fixed gzip timestamp so identical inputs and rules produce identical file
hashes.

## Leakage policy for Stage 5 and later

- IntOGen TCGA_WXS_PAAD evidence is audit-only and cannot define a positive label.
- Unlabelled genes must not be converted wholesale into negatives.
- Any future label source must declare whether its patients overlap TCGA, CPTAC, ICGC or PCAWG.
- CPTAC and ICGC are consumed by this label version through IntOGen and therefore cannot later be
  described as untouched external validation for the same task. Stage 8 must use orthogonal or
  genuinely held-out evidence.
- A label-derived feature or feature-derived label is prohibited unless evaluated with an
  independent cohort or an explicit ablation.
- Final test genes and any external validation source remain unavailable to model selection.

## Completion gates

- Source archives match pinned SHA-256 values and pass ZIP CRC checks.
- All PAAD evidence symbols map uniquely to the GENCODE gene universe.
- Every positive has at least two qualifying external cohort families.
- Canonical controls `KRAS`, `TP53`, `CDKN2A` and `SMAD4` are recovered.
- No TCGA-only evidence can create a positive label.
- No reliable negatives are manufactured while the negative rule is disabled.
- All automated tests pass.
