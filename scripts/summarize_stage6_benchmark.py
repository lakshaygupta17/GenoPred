"""Produce the human-readable catalog and scientific comparison figure from saved results."""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


DESCRIPTIONS = {
 'decision_tree':'A shallow sequence of feature-based decisions.',
 'random_forest':'Averages trees trained with bootstrap samples and feature subsampling.',
 'extra_trees':'Averages trees using more randomized splitting thresholds.',
 'gradient_boosting':'Adds small trees to correct errors left by earlier trees.',
 'hist_gradient_boosting':'Boosting with numerical features grouped into bins.',
 'adaboost':'Repeatedly emphasizes harder training examples using decision stumps.',
 'xgboost':'Regularized gradient-boosted trees.',
 'lightgbm':'Histogram-based gradient boosting with limited leaves.',
 'catboost':'Symmetric boosted trees; all current inputs are numeric.',
 'logistic_l2':'A weighted linear score with shrinkage.',
 'logistic_l1':'A sparse linear score that can set some weights to zero.',
 'ridge':'A regularized linear least-squares classifier.',
 'linear_svm':'A linear boundary fitted with a margin-based objective.',
 'rbf_svm':'A curved boundary based on similarity between gene profiles.',
 'knn':'Scores a gene from nearby training examples in feature space.',
 'gaussian_nb':'Combines feature-wise Gaussian distributions under conditional independence.',
 'lda':'Uses a shared, regularized covariance model for the two temporary classes.',
 'qda':'Uses separate heavily regularized covariance models for the temporary classes.',
 'mlp':'A neural network with one hidden layer containing eight units.',
 'gaussian_process':'A fixed-kernel probabilistic classifier; no kernel optimization.',
}


def main(root):
 report=json.loads((root/'manifests/stage6_benchmark_report.json').read_text())
 folder=root/'data/processed/stage6_benchmark_v1'
 plan=json.loads((folder/'plan.json').read_text())
 names=json.loads((root/'data/processed/stage5/feature_names.json').read_text())
 results=report['results']
 lines=['# Comprehensive Stage 6 experiment catalog','',
 'Completed 4 November 2026. This is a bounded catalog of practical experiments for the current',
 'gene-level dataset, plus a separately identified research backlog. There is no finite exhaustive',
 'list of every algorithm, architecture and parameter combination. The frozen registry makes the',
 'scope of this run explicit and reproducible.','',
 f"The run attempted {len(plan['experiments'])} registered experiments, added one validation-selected",
 f"rank ensemble and compared four earlier references. It produced {len(results)} result rows.",
 f"Recorded failed attempts: {len(report['failures'])}.",'',
 '## How to interpret the experiment','',
 'The task remains ranking genes for pancreatic-cancer research. A high score is a research',
 'hypothesis, not a patient diagnosis or a calibrated probability of causing cancer. The data',
 'contain 17 observed positive genes and 19,945 unlabelled genes, with no verified negative set.',
 'Learning from unlabelled genes requires assumptions. None of these experiments establishes',
 'that a low-ranking gene is harmless.','',
 'The preserved split has 11/3/3 positives in training/validation/test. Selection uses training',
 'genes and the three validation positives; final refitting uses training plus validation. The',
 'test subset contains 3,992 genes. It was examined in earlier work, so all comparisons remain',
 'exploratory. Searching many models increases the chance of finding a flattering result by luck.',
 'The largest test score must not be promoted as a validated winner.','',
 '## All tree and other supervised learners tried with PU bagging','',
 'Each learner below was trained in PU bags: all permitted positive genes plus sampled unlabelled',
 'comparators. The latter receive temporary negative labels for that learner. All candidates use',
 'six bags per seed, three fixed seeds (41, 137, 907), and ten unlabelled examples per positive.',
 'This is a modest screening budget shared by these families. Earlier tree references used 40',
 'bags per seed and a different sampling-ratio search, so they are retained as separate runs.','',
 '| Experiment | What it tries | Status |','|---|---|---|']
 for family,desc in DESCRIPTIONS.items():
  name='pu_'+family
  lines.append(f"| {name} | {desc} | {'Completed' if name in results else 'Failed; see failure log'} |")
 lines+=['','Most models use balanced sample weights. KNN and Gaussian processes do not support',
 'those weights in the chosen interfaces and retain the bag class ratio. LDA and QDA use equal',
 'class priors. SVM and ridge margins are sigmoid-transformed only to put bag outputs on a',
 'bounded scale; that operation is not probability calibration.','',
 'Two settings per family were tried. Tree depth is 2 or 3; boosted models use 60 iterations.',
 'Logistic/SVM strength is 0.1 or 1; ridge penalty is 10 or 1; KNN uses 3 or 7 neighbours;',
 'Gaussian NB smoothing is 0.001 or 0.1; LDA shrinkage and QDA regularization are 0.5 or 0.9;',
 'the eight-unit MLP penalty is 1 or 10; Gaussian-process RBF length scale is 1 or 3.',
 'AdaBoost uses depth-one stumps with learning rates 0.1 or 0.5. Full constructors are in',
 'src/genopred_v2/stage6_benchmark.py and dependency versions are recorded with the outputs.','',
 '## Alternative learning assumptions and graph experiments','',
 '| Experiment | What changes | Interpretation |','|---|---|---|',
 '| biased_logistic_l2 | Fits all unlabelled examples with temporary zero labels | A deliberately biased comparison; not a verified-negative model |',
 '| biased_random_forest | Same comparison using a forest | Tests the effect of sampling versus using the full unlabelled pool |',
 '| two_step_provisional_negatives | Selects the lowest 10 or 30 percent of training-unlabelled pilot scores, then fits a forest | These are provisional comparators, never reliable biological negatives |',
 '| centroid | Distance or cosine similarity to the positive-gene mean | Learns a simple positive profile |',
 '| nearest_positive | Distance to the nearest 1 or 3 positive genes | Tests local similarity without assigning negatives |',
 '| positive_kde | Positive-only kernel density with bandwidth 1 or 3 | High density means similarity to known positives |',
 '| positive_ocsvm | One-class SVM fitted to positives with nu 0.1 or 0.3 | Small positive sample limits stability |',
 '| positive_isolation_forest | Isolation Forest fitted to positives | Uses inlier similarity, not generic anomaly detection on all genes |',
 '| rbf_features_nnpu | 64 random Fourier features followed by the existing nnPU learner | Nonlinear PU risk model; gamma 0.05 or 0.2, prior 0.001, L2 0.01, 300 epochs |',
 '| positive_seed_graph_diffusion | Propagates permitted positive labels through STRING | All graph topology is visible; held-out positive labels never seed propagation |',
 '| graph_degree | Ranks by summed edge confidence | A study/connectivity bias baseline, not proof of cancer relevance |','',
 'Positive-seeded diffusion tries restart values 0.3 and 0.7 with 100 iterations and self-loops.',
 'Graph results use the entire frozen gene network and are transductive. Their saved bundles',
 'require the same gene order and pinned graph. The positive-only Isolation Forest has no',
 'effective level-dependent parameter; its two registered settings are identical, and this',
 'redundancy is recorded rather than counted as additional scientific evidence.','',
 '## Every biological feature tested alone','',
 'For each feature, selection tries high-value-first and low-value-first ranking. Missing values',
 'are placed last in either direction. The chosen direction depends only on validation AP.','',
 '| Experiment | Feature |','|---|---|']
 for i,feature in enumerate(names[:10]):lines.append(f'| feature_{i:02d} | {feature} |')
 lines+=['','## Feature removal and preprocessing experiments','',
 'These use the same PU logistic L2 learner at C=0.1, preserving all other registered settings.','',
 '| Experiment | Inputs or preparation |','|---|---|',
 '| ablation_mutation_only | Mutation prevalence and its missingness flag |',
 '| ablation_cnv_only | Copy gain/loss prevalence and their flags |',
 '| ablation_expression_only | Tumour expression and its flag |',
 '| ablation_regulation_only | Six ATAC/chromatin measurements and their flags |',
 '| ablation_no_mutation | Every input except mutation prevalence and its flag |',
 '| ablation_biology_without_flags | Ten biological columns, excluding missingness flags |',
 '| ablation_missing_flags_only | Ten availability flags, excluding biological values |',
 '| preprocess_robust | Median/interquartile scaling |',
 '| preprocess_quantile | Training-fitted quantile transformation to a normal distribution |',
 '| preprocess_pca | Standard scaling then five principal components |',
 '| preprocess_polynomial | Standard scaling then pairwise interaction features |','',
 'All imputation and transformations are fitted only on the permitted partition. Standard',
 'scaling is used elsewhere, including for tree families for a consistent preparation path.',
 'It is not necessary for trees. Ablations test one reference learner, so they do not establish',
 'that the same feature group is useless for every possible model.','',
 '## Controls and ensemble','',
 '| Experiment | Purpose |','|---|---|',
 '| random_rank_control | Random rankings averaged across fixed seeds |',
 '| constant_score_control | All genes tied, exercising the tie-aware metric policy |',
 '| shuffled_training_labels_control | PU logistic model after shuffling training labels only |',
 '| validation_selected_rank_ensemble | Equal percentile-rank average of the three PU families with highest mean validation AP |','',
 'The rank ensemble members were fixed before computing test metrics: '+', '.join(report['selected'][n]['spec']['name'] for n in results['validation_selected_rank_ensemble']['members'])+'.',
 'Percentile ranks are computed over the frozen full gene universe. The ensemble has no',
 'test-selected weights. Validation AP displayed for the ensemble is computed from the',
 'combined scores; other new-model validation AP values are means across seed-specific APs.','',
 'Four references were loaded from existing fitted artifacts: linear nnPU, nnPU plus STRING',
 'diffusion, the earlier PU tree ensemble and the earlier PU Random Forest ensemble.','',
 '## Checks and evidence limits','',
 'The run verified the corrected Stage 5 hashes and exact equality with the original gene split.',
 'Model settings and ensemble members were saved before test evaluation. Every new fitted model',
 'was saved and reloaded; predictions were checked to absolute tolerance 1e-12. Test metrics',
 'use threshold-grouped average precision and fractional credit for ties at top-k cutoffs.',
 'The project suite passes 50 tests, including held-out-label exclusion, training-only',
 'preprocessing, graph-seed exclusion, missing-value ranking, and feature-order checks.','',
 'Per-seed test metrics and rank correlations measure random-seed sensitivity on this one split.',
 'They are not confidence intervals. Three additional stratified folds on the development genes',
 'were evaluated for the three validation-selected families using fixed settings and seed 41.',
 'That check reuses development data involved in model selection; it is a sensitivity diagnostic,',
 'not nested cross-validation or independent confirmation. No test genes enter those folds.','',
 'The gene score exports also contain in-sample scores for training genes. Use the split column',
 'when interpreting them; a training-positive rank is not held-out evidence. Warning messages',
 'are retained in warnings.json, including covariance or convergence issues.','',
 '## Research backlog beyond this completed screen','',
 'These items are listed explicitly as not run. They require new data, a new experimental design,',
 'or a separately specified implementation; they are not implied successes of this benchmark.','',
 '| Direction | What is needed before a meaningful experiment |','|---|---|',
 '| GCN, GraphSAGE and graph attention models with PU objectives | Define graph supervision, graph leakage controls and a training protocol; a learned graph model is not implemented by diffusion |',
 '| Deep tabular transformers and large neural networks | Much more labelled evidence or a justified pretraining source; an eight-unit MLP was the small-network screen here |',
 '| Breast-to-pancreas and pan-cancer transfer | Compatible source features, source labels, overlap audit and negative-transfer controls; belongs to Stage 7 |',
 '| Hi-C, DepMap, pathways, embeddings and extra omics | Acquire, harmonize and independently audit the new measurements before changing the frozen matrix |',
 '| Elkan-Noto propensity correction and calibrated PU probabilities | Assess its positive-labelling assumptions and obtain enough independent positives for a stable calibration holdout |',
 '| Selection-biased PU learning and class-prior estimation | Evidence about how genes enter the known-positive list and uncertainty in the hidden-positive prevalence |',
 '| Positive Bayesian networks and discrete PU Naive Bayes | Predeclare discretization and sparse-count handling; continuous Gaussian NB bagging was tested here |',
 '| Gaussian-mixture and manifold representations | Predeclare a train-only representation protocol and validate it without extending this reused-test search |',
 '| Multitask, survival and patient-level prediction | Appropriate patient outcomes and a patient-level dataset; current rows represent genes |',
 '| Nested repeated cross-validation and label-sensitivity analysis | Predeclare a new assessment protocol; the existing test is consumed and cannot become unseen again |',
 '| Temporal, chromosome-grouped or gene-family-grouped validation | Define scientifically justified grouping or timestamp evidence and adequate positive coverage |',
 '| Independent functional validation | Obtain evidence not already used to construct labels, including suitable perturbation or held-out study data |',
 '| SHAP, permutation importance and biological follow-up | Choose a candidate based on a defensible protocol first; explanations of this screen are exploratory associations |',
 '| Broad hyperparameter optimization and stacked meta-models | A larger development/validation evidence base; further searches on three positives would amplify selection noise |','',
 '## Reproduce and inspect','',
 '```powershell',
 '.\\.venv\\Scripts\\python.exe -m pip install -r data/processed/stage6_benchmark_v1/requirements-lock.txt',
 '.\\.venv\\Scripts\\python.exe -m genopred_v2.stage6_benchmark .',
 '.\\.venv\\Scripts\\python.exe scripts/summarize_stage6_benchmark.py .',
 '.\\.venv\\Scripts\\python.exe -m unittest discover -s tests -q',
 '```','',
 'Run commands from the project root. The registry in the module is authoritative. The run writes',
 'plan.json before fitting and selection_frozen.json before test evaluation. A repeated run',
 'replaces this benchmark output folder contents by filename; preserve a copy before changing',
 'settings. It does not overwrite earlier Stage 6 or tree benchmark model files.','',
 'The comparison is in STAGE6_BENCHMARK_RESULTS.md and manifests/stage6_benchmark_report.json.',
 'All scores, model bundles, settings, dependency versions, warning records, and test-positive',
 'ranks are in data/processed/stage6_benchmark_v1. Joblib uses Python pickle: load only trusted',
 'local bundles. General model bundles store preprocessing and feature order. Graph and',
 'percentile-rank ensemble reproduction also requires the pinned gene universe.','',
 '## Method references','',
 '- [Scikit-learn model families](https://scikit-learn.org/stable/supervised_learning.html)',
 '- [PU bagging implementation and references](https://pulearn.github.io/pulearn/doc/pulearn/bagging.html)',
 '- [Mordelet and Vert on PU bagging](https://arxiv.org/abs/1010.0772)',
 '- [XGBoost API](https://xgboost.readthedocs.io/en/stable/python/python_api.html)',
 '- [LightGBM classifier](https://lightgbm.readthedocs.io/en/stable/pythonapi/lightgbm.LGBMClassifier.html)',
 '- [CatBoost classifier](https://catboost.ai/docs/en/concepts/python-reference_catboostclassifier)','']
 findings=['## What the completed screen showed','',
 'The largest observed test AP was 0.66721 for shallow PU decision-tree bagging, followed by',
 '0.55641 for gradient boosting and 0.55627 for LightGBM. These are descriptive results after',
 'a large search. They do not establish independent superiority over nnPU.','',
 'The shallow tree ensemble ranked KMT2C first and RBM10 second among the test genes, but ACVR2A',
 'was tied across positions 1,017 to 1,855. The original nnPU placed those genes at 3, 6 and 412.',
 'This explains why tree AP can rise while recovery within the top 500 falls.','',
 'For the fixed PU logistic ablation, removing mutation data gave AP 0.00163, compared with',
 '0.11940 for mutation alone. This supports mutation dependence in that reference learner;',
 'it does not prove other measurements are biologically unimportant.','',
 'The three validation-selected families were KNN, the small MLP and AdaBoost. Their rank',
 'ensemble reached test AP 0.20221, below original nnPU at 0.22465. Validation leadership did',
 'not transfer consistently to the test subset. Development-fold AP ranged from 0.20467 to',
 '0.67857 for KNN, 0.30851 to 0.83692 for MLP, and 0.24894 to 0.50069 for AdaBoost. This',
 'variation reinforces the need for a stronger assessment before choosing a final model.','',
 'QDA completed but emitted 54 positive-class covariance-rank warnings across selection and',
 'refitting. Its results should be treated cautiously. All other recorded warning groups were',
 'empty in this run. The exhaustive warning log is saved alongside the models.','']
 insertion=lines.index('## How to interpret the experiment')
 lines[insertion:insertion]=findings
 (root/'docs/STAGE6_EXPERIMENT_CATALOG.md').write_text('\n'.join(lines),encoding='utf-8')
 # Plot every primary PU learner against original references; no selection by test results.
 plotted=[n for n in results if n.startswith('pu_')]+[n for n in results if n.startswith('reference_')]
 plotted=sorted(plotted,key=lambda n:results[n].get('validation_mean_ap',-1),reverse=True)
 labels=[n.replace('reference_','Earlier ').replace('pu_','PU ').replace('_',' ') for n in plotted]
 fig,ax=plt.subplots(figsize=(11,10))
 ypos=np.arange(len(plotted)); vals=[results[n]['test']['observed_label_auprc'] for n in plotted]
 ax.barh(ypos,vals,color=['#748398' if n.startswith('reference_') else '#257D89' for n in plotted])
 ax.set_yticks(ypos,labels);ax.invert_yaxis();ax.set_xlim(0,1)
 ax.set_xlabel('Observed-label test average precision')
 ax.set_title('Stage 6 exploratory model comparison\nPrimary PU models ordered by validation score, with earlier references',loc='left',pad=15)
 for i,v in enumerate(vals):ax.text(v+.008,i,f'{v:.3f}',va='center',fontsize=8)
 ax.spines[['top','right']].set_visible(False)
 fig.text(.04,.015,'Only three known test positives. Previously examined split. Scores are not clinical accuracy.',fontsize=10)
 fig.tight_layout(rect=(0,.04,1,1))
 figures=root/'docs/figures';figures.mkdir(exist_ok=True)
 fig.savefig(figures/'stage6_model_comparison.png',dpi=170);plt.close(fig)
 print('Catalog and comparison figure created')


if __name__=='__main__':main(Path(sys.argv[1]).resolve())
