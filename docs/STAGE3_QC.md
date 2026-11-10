# Stage 3 QC and Harmonization Decisions

Stage 3 converts verified repository files into explicit, testable data contracts. Raw files remain
immutable. No signal is aggregated into gene features until its identifiers and coordinate system
pass the gates below.

Last updated: **2026-08-26**

## Current audit result

| Input | Files | Schema result | Build gate |
|---|---:|---|---|
| TCGA-PAAD STAR counts | 183 | One nine-column GENCODE v36 schema | GRCh38 compatible |
| TCGA-PAAD ASCAT3 gene CNV | 179 | One gene-level CNV schema | GRCh38 compatible |
| TCGA-PAAD masked MAF | 180 | One 140-column GDC MAF schema | GRCh38 compatible |
| GSE149103 bigWig/Hi-C | 18 | 16 chromosome trees and 2 archive inventories valid | hg19; controlled query/liftover policy |
| GSE124229 ATAC atlas | 1 | Three coordinate plus 80 assay columns | Materialized in GRCh38 |

The TCGA files map to 185 cases. After enforcing Primary Tumor RNA and exactly one selected CNV and
MAF file, there are 178 RNA cases, 179 CNV cases, 178 mutation cases and 165 complete cases. All
patient-level train/validation/test splitting must occur on the case identifier before any
preprocessing is fitted.

## TCGA RNA contract

- Source schema declares GENCODE v36.
- Remove STAR summary rows beginning with `N_` before constructing gene matrices.
- Preserve the original versioned Ensembl ID in provenance; use an explicit unversioned ID only as
  a separate join key.
- Keep raw unstranded counts for count-aware QC and use a declared transformation such as
  `log2(TPM + 1)` for continuous model features.
- Do not combine primary, normal and metastatic aliquots without a declared sample-selection rule.
- Duplicate files for one case must be resolved before patient-level aggregation.

## TCGA copy-number contract

- Use the selected ASCAT3 gene-level workflow only; supporting CNV workflows cannot be silently
  mixed into the primary matrix.
- Retain amplification and deletion evidence as separate features.
- Treat blank copy-number values as missing, not as neutral copy number.
- Use Ensembl gene ID as the primary join key and retain chromosome/start/end for build auditing.

## TCGA mutation contract

- MAF records declare `NCBI_Build=GRCh38`.
- Aggregate mutations using `Tumor_Sample_Barcode` from each record. File-level GDC associations
  include matched-normal material and therefore cannot determine the tumour observation by
  themselves.
- Define the accepted consequence set before calculating prevalence.
- Keep synonymous, noncoding and low-confidence calls available for QC, but do not silently include
  them in damaging-mutation prevalence.
- Calculate recurrence and prevalence inside training data when they influence feature selection.

## GEO coordinate contract

- GSE149103 and GSE124229 are hg19/GRCh37 and cannot join the GRCh38 matrix directly.
- Record the exact chain-file URL and SHA-256 before liftover.
- For intervals, lift both boundaries and reject split, inverted, cross-chromosome or unmapped
  results unless an assay-specific policy explicitly handles them.
- Report the fraction of intervals and signal mass retained per file.
- Never replace a failed liftover with the original hg19 coordinate.

## Completed gates

1. Parsed all 16 bigWig chromosome trees with a portable read-only parser; shared chromosome sizes
   are consistent with hg19 and tracks may omit chromosomes with no signal.
2. Stream-inventoried both Hi-C tar archives: each contains 24 chromosome `.hic` files plus 24
   harmless macOS metadata entries.
3. Pinned the UCSC hg19-to-hg38 chain; MD5 `35887f73fe5e2231656504d1f6430900`
   and SHA-256 `5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d`.
4. Lifted 119,299/119,466 GSE124229 intervals (99.860%) while retaining 98.442% of signal mass;
   rejected 167 ambiguous, split or strand-discordant intervals.
5. Built the 60,660-record GENCODE v36 identity map. Versionless Ensembl IDs retain `_PAR_Y`;
   pseudoautosomal pairs and reused symbols are never silently collapsed.
6. Verified 100% identity-map coverage for representative GDC RNA and CNV inputs.
7. Selected Primary Tumor RNA only and produced a 165-case complete three-modality cohort without
   duplicate tumour-RNA assignments.

## Generated evidence

- `manifests/stage3_metadata_report.json`
- `manifests/stage3_schema_report.json`
- `manifests/stage3_file_inventory.csv`
- `manifests/stage3_sample_inventory.csv`
- `manifests/stage3_reference_report.json`
- `manifests/GSE149103_container_qc.json`
- `manifests/GSE124229_liftover_report.json`
- `manifests/gene_identity_report.json`
- `manifests/TCGA-PAAD_gene_identity_coverage.json`
- `manifests/TCGA-PAAD_canonical_cohort.csv`
- `manifests/stage3_completion_report.json`

These files are regenerated with `genopred-v2 acquire-stage3-metadata`, `genopred-v2 stage3-qc`
and `genopred-v2 complete-stage3`.
