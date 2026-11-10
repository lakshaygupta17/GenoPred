# Stage 6: leakage-controlled baselines

Current experiment: `stage6-baselines-v2`, using `stage5-core-v2`. Version 1 is superseded;
its code, results and documentation are preserved in `audit_archive/pre_correction_v1`.

## Outcome

Stage 6 establishes the first from-scratch PAAD ranking benchmarks. It does not use breast-cancer
pretraining. The primary baseline is a linear non-negative positive-unlabelled (nnPU) model; the
graph baseline diffuses its scores over the frozen STRING physical network.

The feature-only model has higher average precision than graph diffusion on the previously examined
test split; graph AUROC is slightly higher. These correction results are exploratory, not a fresh
test of generalization.

## Preserved split (previously examined)

Genes were assigned by a seeded SHA-256 ordering, stratified only by Stage 4 label state. Features
were never consulted when splitting.

| Split | Known positives | Unlabelled genes |
|---|---:|---:|
| Training | 11 | 11,967 |
| Validation | 3 | 3,989 |
| Test | 3 | 3,989 |

The test labels were not used for scaling, fitting, class-prior selection, regularization selection
or graph restart selection. Scaling was fitted first on training genes, then refitted on
training-plus-validation genes only after all hyperparameters were frozen.

Stage 5's global imputed values are only finite storage placeholders. Stage 6 uses the companion
missing flags to replace them with medians fitted on the permitted training partition before it
fits that partition's scaler, preventing test-distribution leakage through imputation.

## Models

### Feature-only nnPU logistic baseline

The model minimizes non-negative PU risk rather than declaring every unlabelled gene a true
negative. Class priors of 0.001, 0.003 and 0.01 and L2 penalties of 0.0001 and 0.01 were compared on
validation observed-label AUPRC. The selected prior was 0.001 and the selected L2 penalty was 0.01.

Coefficients are for standardized inputs and are associations, not causal effects. The full ordered
coefficient list is stored in the report.

### STRING random-walk baseline

The graph baseline applies a weighted random walk with restart to the feature-only scores. STRING
combined scores are the edge weights and every node receives a self-loop, so isolated genes retain
a prediction. Restart values 0.3, 0.5, 0.7 and 0.9 were compared using validation labels; 0.5 was
selected for the corrected dataset.

The graph experiment is transductive: topology and feature values for all genes are visible.
Validation labels select parameters and enter the final training-plus-validation refit. Test labels
do not enter fitting; however, they were already examined in the v1 run and audit.

## Exploratory corrected results

| Metric | Feature-only nnPU | STRING random walk |
|---|---:|---:|
| Observed-label AUROC | 0.9653 | 0.9658 |
| Observed-label average precision | 0.22465 | 0.08844 |
| AP / observed positive prevalence | 298.93x | 117.69x |
| Median positive percentile | 0.9987 | 0.9947 |
| Known-positive recall at 100 | 2/3 | 2/3 |
| Known-positive recall at 500 | 3/3 | 3/3 |

The positive prevalence reference is 0.0007515. This reference is not the exact finite-sample
expectation of AP under random rankings. Three positives (`KMT2C`, `ACVR2A`, `RBM10`) are too few
for a stable estimate, and this split is not an untouched holdout. Top-k metrics refer to the
3,992-gene test subset, not the full 19,962-gene ranking.

The isolated positive `LRP1B` has global ranks 24 (features) and 29 (graph). Its score is unchanged
by diffusion; its rank changes because scores of other genes move. It is a training positive, so
this is a diagnostic rather than held-out sensitivity evidence.

Average precision groups tied scores at common thresholds. Recall at a cutoff allocates a
fractional share to tied genes; percentile uses midranks. Thus metrics do not depend on gene order
inside ties. The legacy `observed_label_auprc` field stores non-interpolated average precision.

## Interpretation limits

- Unlabelled genes may contain undiscovered drivers. AUROC and AUPRC therefore measure recovery of
  the 17 observed labels, not discrimination against verified biological negatives.
- Three test positives are too few for a stable performance claim or clinical interpretation.
- This is one fixed split. The result is a software and modelling baseline, not a publication-grade
  estimate of generalization.
- Feature coefficients and high-ranking unlabelled genes are hypotheses, not validated drivers.
- Graph topology can encode study and annotation bias toward well-characterized genes.
- The Stage 6 test has now been evaluated and is consumed. Stage 7 must not tune transfer learning
  against it or repeatedly reuse it as if untouched; later claims require a newly locked or external
  validation design.

## Outputs

- `configs/stage6.json`: frozen split and hyperparameter search specification
- `data/processed/stage6/gene_scores.tsv.gz`: gene-level scores, ranks, labels and split membership
- `data/processed/stage6/*.npy`: fitted weights, scaler, scores and split tensor
- `data/processed/stage6/model.npz`: authoritative float64 weights, intercept, imputer, scaler,
  feature order and graph diffusion settings; loaders verify feature order and reproduce both models
- `manifests/stage6_baseline_report.json`: complete selection trail, coefficients, metrics and hashes
- `manifests/stage6_completion_report.json`: Stage 6 completion gate

The implementation uses only NumPy. Input hashes are verified before modelling; output manifests
record the source dataset and configuration hash. Both saved-model round trips have zero maximum
absolute prediction error. The correction suite passes 41 tests.
