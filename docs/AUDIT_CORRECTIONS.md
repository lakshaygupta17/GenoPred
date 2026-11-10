# Audit corrections: Stages 5 and 6

The Stage 1-6 audit found six substantive defects. The code and data were corrected and rebuilt
as `stage5-core-v2` and `stage6-baselines-v2`. The earlier v1 completion and performance claims are
superseded. Historical code, reports, documentation and outputs remain recoverable under
`audit_archive/pre_correction_v1`. Raw data and Stage 4 labels were preserved.

| Finding | Correction | Verified evidence |
|---|---|---|
| MAF parsed as comma-separated data | Explicit tab delimiter, mandatory columns, tumour/case validation and per-case deduplication | 3,525 genes now have nonzero mutation prevalence; KRAS 94/165, TP53 90/165, SMAD4 34/165, CDKN2A 26/165 |
| Promoter liftover accepted enormous spans | Length ratio 0.8-1.25, orientation checks, unique monotonic interior checkpoints every 100 bp, local distortion limit 100 bp | ANKRD20A1 and TBC1D3E are now marked missing instead of querying megabase spans |
| 80 ATAC libraries treated as 80 patients | Two library log1p CPM profiles combined by median within each patient, followed by 10% trimmed cohort mean and patient MAD | Exactly 40 patient profiles; synthetic within-patient noise no longer inflates between-patient variability |
| Absolute CNV cutoffs ignored baseline | Per-case autosomal protein-coding mode, larger mode on ties; amplification >=2x baseline, loss below baseline | Case baselines: 133 at 2, 9 at 3, 22 at 4, 1 at 6; sex-chromosome CNV remains unavailable |
| Tied scores depended on gene order | Threshold-grouped average precision, fractional recall at boundary ties, midrank percentile | Constant-score AP equals prevalence under both tested gene orders |
| Saved model omitted intercept | Complete float64 NPZ bundle and feature/graph prediction loaders | Maximum reload prediction error is 0.0 for both models |

## Verification

The full suite passes 41 tests. Six new tests exercise MAF delimiters and case deduplication,
per-case CNV and sex-chromosome exclusion, patient aggregation, invalid promoter mappings, tied
metrics, and saved-model bias/feature order. Stage 5's finite matrix has 19,962 rows and 20 columns.
Stage 6 checks the Stage 5 version and every recorded input artifact hash before fitting; it records
input hashes and the experiment configuration hash in its report.

Copy-number handling follows the relative-copy principle in the
[GDC category documentation](https://gdc.cancer.gov/how-are-five-categories-copy-number-changes-determined).
The modal baseline is an integer approximation, not ASCAT's precise tumour-ploidy estimate.
Sex-chromosome CNV is explicitly missing pending sex-aware processing. This increases missing CNV
to 953 genes. The top 100 feature-model genes now include 5 on chromosome X versus 68 previously;
this diagnostic does not prove that all chromosome or annotation bias has been removed.

## Corrected exploratory results

| Model | Observed-label AUROC | Average precision | Recall at 100 | Recall at 500 |
|---|---:|---:|---:|---:|
| Feature-only nnPU | 0.965321 | 0.224649 | 2/3 | 3/3 |
| STRING diffusion | 0.965823 | 0.088443 | 2/3 | 3/3 |

These use the unchanged 11/3/3 positive train/validation/test allocation and the same candidate
parameter grid. The test contains 3 positives and 3,989 unlabelled genes. It was already examined
before correction, so these are exploratory rerun results and cannot establish independent
generalization. Restoring mutation prevalence explains why the results can differ substantially
from v1; this is not evidence from a new validation cohort.

Average precision is reported under the legacy `observed_label_auprc` key. Its ratio to observed
positive prevalence is a descriptive reference, not an exact finite-sample random-AP expectation.
Top-k metrics are computed within the test subset. LRP1B is a training positive: its feature score
is unchanged by diffusion, while its global rank moves from 24 to 29 as other scores change.

## Remaining boundaries

No verified-negative gene set exists; unlabelled genes remain unknown. Gene-specific mutation
callability masks are still unavailable. Conservative liftover checkpoints are not full basewise
alignment validation. Hi-C remains deferred and DepMap remains reserved for orthogonal validation.
The graph baseline is diffusion, not a trained graph neural network. A later scientific claim
requires genuinely held-out or orthogonal evidence; reshuffling the same examined genes is not
sufficient. None of these limits is hidden by a software-completion flag.
