# Stage 6 tree benchmark

Completed 4 November 2026 using corrected Stage 5 v2 and the preserved Stage 6 split.

## Results

| Method | Observed label average precision | AUROC | Positives in top 100 | Positives in top 500 |
|---|---:|---:|---:|---:|
| Linear nnPU | 0.22465 | 0.9653 | 2 of 3 | 3 of 3 |
| PU tree bagging | 0.41727 | 0.8608 | 2 of 3 | 2 of 3 |
| PU Random Forest bagging | 0.37870 | 0.9084 | 2 of 3 | 2 of 3 |

Trees improved concentration near the top as measured by average precision, but lost retrieval
of one positive within the top 500. Keep all three benchmarks: these results do not justify
declaring a universally superior replacement. Metrics refer to the same 3,992 test genes.

## Learning and selection

Every bag contains all permitted positive examples and a randomly sampled set of unlabelled
comparators. The comparators receive temporary zero labels for that learner; they are not
verified biological negatives. Balanced class weights are used. Scores average the positive
class outputs of the learners and are not calibrated driver probabilities.

For each family, depths 2 and 3 and unlabelled-to-positive sampling ratios 5 and 20 were compared.
All candidates used seeds 41, 137 and 907, 40 bags per seed and minimum leaf size 2.
Each forest bag contains 10 trees with sklearn's default square-root feature subsampling.
The standalone decision trees consider all features. The selection criterion was mean validation
average precision across the three seeds; ties prefer shallower trees and the smaller ratio.
All choices were written to selection.json before test evaluation.

Imputation uses only training examples during selection, and training plus validation during
the final refit. Trees do not require feature scaling. No test labels enter fitting or selection.
Final scores average all three seeds, rather than choosing the best test seed.
Training uses 11 positives during selection and 14 during final refitting.

## Stability and limitations

Individual-seed test average precision ranged from 0.36061 to 0.55620 for tree bagging and
0.37118 to 0.39496 for forest bagging. These ranges describe random-seed sensitivity on the
same split; they are not confidence intervals or evidence of independent validation.
Only three known positives are in the test group, which has already been examined in earlier
experiments. All results are exploratory. Unlabelled genes may contain undiscovered drivers.
Scores on training genes in the full gene table are in-sample and must not be used as held-out
performance evidence. Trees were evaluated without graph diffusion in this experiment.

## Reproduction and saved artifacts

Run from the project directory:

```powershell
.\.venv\Scripts\python.exe -m pip install scikit-learn==1.7.2
.\.venv\Scripts\python.exe -m genopred_v2.stage6_trees .
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

Implementation: src/genopred_v2/stage6_trees.py. The module's CONFIG is authoritative;
configs/stage6_trees.json is its exported record, not an editable input configuration.
The run verifies Stage 5 hashes and equality with the original split. The report records versions,
input hashes, configuration and code hashes, all validation candidates, selected settings,
seed results, and final metrics.

Outputs are in data/processed/stage6_trees. The two joblib bundles contain models, imputation
medians, feature order and settings. gene_scores.tsv.gz contains scores in original gene order;
sort scores descending to obtain a ranking. selection.json records frozen choices.
Detailed comparison: manifests/stage6_tree_report.json.

Both saved bundles reproduced predictions with zero maximum absolute error. All 44 tests passed,
including deterministic ensemble fitting, feature-order rejection, PU group checks, and training-only
imputation. Only load trusted local joblib files; serialization is Python pickle based.

This extends Stage 6. It does not complete Stage 7 or establish external biological validity.
