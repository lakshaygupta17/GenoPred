# Stage 5: model-ready dataset

Stage 5 freezes the first leakage-controlled GenoPred V2 input dataset. It contains gene-level
features and a physical interaction graph, but it does not train or tune a model.

## Frozen scope

- Dataset version: `stage5-core-v2` (supersedes the defective v1 dataset)
- Universe: 19,962 GENCODE v36 protein-coding genes from the Stage 4 label table
- Labels: 17 positive and 19,945 unlabelled genes; unlabelled does not mean negative
- Cohort: 165 TCGA-PAAD primary-tumour cases having RNA, CNV and mutation files
- Matrix: 19,962 rows by 20 columns (10 biological features plus 10 explicit missing flags)
- Graph: STRING v12.0 human physical links at combined score >=700

Rows use the same deterministic `node_index` in the feature table, label tensor and graph. No
breast-cancer data or V1 labels enter this dataset.

## Biological features

| Feature | Calculation |
|---|---|
| PANC-1 vs HPNE ATAC | Difference of replicate-median `log1p` promoter signals |
| Patient ATAC mean | 10% trimmed mean across 40 patient profiles, each the median of two library `log1p(CPM)` profiles |
| Patient ATAC variability | Median absolute deviation across those 40 patient profiles |
| PANC-1 vs HPNE CTCF | Difference of replicate-median signed promoter signals |
| PANC-1 vs HPNE H3K4me3 | Difference of replicate-median signed promoter signals |
| PANC-1 vs HPNE H3K27ac | Difference of replicate-median signed promoter signals |
| Damaging SNV prevalence | Fraction of 165 cases with at least one selected damaging consequence |
| Amplification prevalence | Autosomes: fraction of cases with copy number >= twice the case's modal baseline |
| Deletion prevalence | Autosomes: fraction of cases below the case's modal baseline, including zero (loss) |
| Tumour expression | Median `log2(TPM+1)` across the 165 cases |

Promoter features use TSS +/-2 kb. GSE149103 hg19 windows are reached by unique reverse liftover
from the canonical GRCh38 gene coordinates. Ambiguous or unavailable mappings are missing, not
treated as zero.

Mapped windows must preserve orientation, have a length ratio between 0.8 and 1.25, and pass
unique, monotonic interior mappings every 100 bp. A local displacement error over 100 bp rejects
the window. This is a conservative filter, not proof of complete base-level alignment.

CNV baselines use the mode across observed autosomal protein-coding genes (larger mode on ties).
Sex-chromosome CNV values are marked missing pending validated sex-aware baselines. The definition
follows the relative-copy principle in the [GDC CNV category methodology](https://gdc.cancer.gov/how-are-five-categories-copy-number-changes-determined).

## Missing data and leakage controls

Every base feature has a companion `__missing` column. Non-finite or structurally unavailable
values are median-imputed across genes solely to create a finite matrix; the flag preserves the
fact that the value was imputed. Current missing counts are listed in the dataset report; they now
include unavailable sex-chromosome CNV and rejected promoter mappings.

Stage 5 performs no z-scoring, supervised feature selection, resampling or label-dependent
imputation. Stage 6 must fit any scaler or selector on training nodes only. TCGA mutation evidence
was excluded from Stage 4 label assignment because mutation prevalence is now an input feature.
DepMap remains validation-only.

The mutation prevalence denominator is the 165-case available cohort. Gene-specific exome
callability masks were unavailable, so this feature requires sensitivity analysis.

## Physical graph

The graph uses the immutable STRING v12.0 human physical-links and protein-info files under CC BY
4.0. Preferred symbols are mapped only when they identify exactly one Stage 4 gene. After filtering
at score 700, the graph has 85,171 undirected interactions, represented as 170,342 directed edges,
and 10,599 connected nodes.

There are 9,363 isolated genes. One positive, `LRP1B`, is isolated at this confidence threshold.
The threshold was not weakened to force positive connectivity. Stage 6 graph models must add
self-loops, preserve feature-only predictions for isolated nodes and report an isolated-positive
sensitivity analysis.

## Deferred inputs

- Hi-C local and distal promoter features: deferred until a validated Windows-compatible `.hic`
  extraction route and the independent GSE185069 genome-build declaration are available.
- DepMap pancreatic CRISPR dependency: held out for orthogonal validation rather than used as a
  training feature.

These deferrals are explicit; no zero-filled or synthetic substitute columns were created.

## Outputs

- `data/processed/stage5/gene_features.tsv.gz`: auditable gene table, labels and denominators
- `data/processed/stage5/x.npy`: finite float32 feature matrix
- `data/processed/stage5/y.npy`: positive/unlabelled label tensor (`1`/`-1`)
- `data/processed/stage5/edge_index.npy`: directed graph indices
- `data/processed/stage5/edge_attr.npy`: STRING scores divided by 1,000
- `data/processed/stage5/string_physical_edges.tsv.gz`: human-readable mapped edges
- `manifests/stage5_dataset_report.json`: dimensions, missingness, graph QC and SHA-256 hashes

The correction suite includes real-format MAF regression tests and explicit patient, CNV and
mapping tests. The matrix must contain no NaN or infinity after explicit imputation, and an
entirely zero mutation feature blocks completion. Rebuild with `python -m genopred_v2 complete-stage5`.

The pre-correction implementation, documentation and outputs are preserved under
`audit_archive/pre_correction_v1`. Its mutation values and resulting baseline claims are superseded.
