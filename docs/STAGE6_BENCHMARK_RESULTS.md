# Broad Stage 6 benchmark results

All results are exploratory on the reused 3,992-gene test subset with three known positives.
Rows are sorted by test average precision for inspection, not used to select a deployment model.

| Experiment | Validation AP | Test AP | Test AUROC | Recall at 100 | Recall at 500 |
|---|---:|---:|---:|---:|---:|
| pu_decision_tree | 0.35446 | 0.66721 | 0.88025 | 0.667 | 0.667 |
| pu_gradient_boosting | 0.39642 | 0.55641 | 0.90215 | 0.667 | 0.667 |
| pu_lightgbm | 0.41052 | 0.55627 | 0.88293 | 0.667 | 0.667 |
| pu_catboost | 0.43463 | 0.44540 | 0.91201 | 0.667 | 0.667 |
| pu_random_forest | 0.04230 | 0.42973 | 0.92772 | 0.667 | 0.667 |
| reference_pu_tree | prior run | 0.41727 | 0.86078 | 0.667 | 0.667 |
| pu_xgboost | 0.63028 | 0.40109 | 0.92304 | 0.667 | 0.667 |
| reference_pu_random_forest | prior run | 0.37870 | 0.90841 | 0.667 | 0.667 |
| pu_adaboost | 0.66750 | 0.37552 | 0.83956 | 0.667 | 0.667 |
| biased_random_forest | 0.01524 | 0.36537 | 0.95028 | 0.667 | 0.667 |
| pu_hist_gradient_boosting | 0.57477 | 0.35026 | 0.87223 | 0.667 | 0.667 |
| preprocess_polynomial | 0.63038 | 0.33440 | 0.68472 | 0.333 | 0.333 |
| positive_kde | 0.00947 | 0.33389 | 0.49394 | 0.333 | 0.333 |
| pu_qda | 0.03195 | 0.27843 | 0.87106 | 0.667 | 0.667 |
| ablation_biology_without_flags | 0.66721 | 0.26293 | 0.91802 | 0.667 | 0.667 |
| pu_logistic_l1 | 0.66733 | 0.26263 | 0.88510 | 0.667 | 0.667 |
| reference_nnpu | prior run | 0.22465 | 0.96532 | 0.667 | 1.000 |
| pu_lda | 0.66723 | 0.22440 | 0.96131 | 0.667 | 1.000 |
| pu_gaussian_nb | 0.66739 | 0.22297 | 0.88752 | 0.667 | 0.667 |
| pu_mlp | 0.66757 | 0.21837 | 0.80254 | 0.667 | 0.667 |
| pu_logistic_l2 | 0.66723 | 0.21479 | 0.83296 | 0.667 | 0.667 |
| validation_selected_rank_ensemble | 0.66880 | 0.20221 | 0.81507 | 0.667 | 0.667 |
| biased_logistic_l2 | 0.66740 | 0.18571 | 0.83805 | 0.667 | 0.667 |
| centroid | 0.66714 | 0.17700 | 0.95187 | 0.667 | 1.000 |
| pu_knn | 0.67105 | 0.17161 | 0.78658 | 0.333 | 0.667 |
| positive_ocsvm | 0.33469 | 0.16799 | 0.79067 | 0.333 | 0.333 |
| pu_linear_svm | 0.66723 | 0.12547 | 0.82226 | 0.667 | 0.667 |
| pu_gaussian_process | 0.66694 | 0.12297 | 0.94276 | 0.667 | 0.667 |
| feature_06 | 0.55581 | 0.11940 | 0.96549 | 0.667 | 0.889 |
| ablation_mutation_only | 0.55581 | 0.11940 | 0.96549 | 0.667 | 0.889 |
| reference_nnpu_graph | prior run | 0.08844 | 0.96582 | 0.667 | 1.000 |
| feature_00 | 0.00158 | 0.06297 | 0.81324 | 0.667 | 0.667 |
| two_step_provisional_negatives | 0.06116 | 0.04704 | 0.97209 | 0.333 | 1.000 |
| positive_seed_graph_diffusion | 0.01419 | 0.02606 | 0.89571 | 0.667 | 0.667 |
| pu_rbf_svm | 0.66701 | 0.02416 | 0.95596 | 0.667 | 1.000 |
| preprocess_quantile | 0.38215 | 0.01454 | 0.95137 | 0.333 | 1.000 |
| rbf_features_nnpu | 0.06423 | 0.00638 | 0.88702 | 0.333 | 0.667 |
| graph_degree | 0.00573 | 0.00468 | 0.87524 | 0.000 | 0.667 |
| pu_extra_trees | 0.00595 | 0.00458 | 0.83463 | 0.000 | 0.333 |
| feature_05 | 0.00175 | 0.00431 | 0.69867 | 0.000 | 0.333 |
| feature_01 | 0.00300 | 0.00293 | 0.78971 | 0.000 | 0.333 |
| ablation_missing_flags_only | 0.00076 | 0.00254 | 0.65050 | 0.203 | 0.393 |
| ablation_regulation_only | 0.00120 | 0.00246 | 0.79226 | 0.000 | 0.000 |
| preprocess_pca | 0.00128 | 0.00197 | 0.72591 | 0.000 | 0.333 |
| positive_isolation_forest | 0.00125 | 0.00185 | 0.52887 | 0.000 | 0.333 |
| feature_09 | 0.00217 | 0.00185 | 0.69048 | 0.000 | 0.333 |
| ablation_expression_only | 0.00217 | 0.00185 | 0.69048 | 0.000 | 0.333 |
| nearest_positive | 0.00413 | 0.00172 | 0.44531 | 0.000 | 0.333 |
| pu_ridge | 0.11639 | 0.00170 | 0.71538 | 0.000 | 0.000 |
| ablation_no_mutation | 0.00156 | 0.00163 | 0.70452 | 0.000 | 0.000 |
| feature_02 | 0.00205 | 0.00161 | 0.66788 | 0.000 | 0.000 |
| preprocess_robust | 0.00149 | 0.00152 | 0.68346 | 0.000 | 0.000 |
| feature_04 | 0.00120 | 0.00097 | 0.47514 | 0.000 | 0.000 |
| ablation_cnv_only | 0.00216 | 0.00097 | 0.43282 | 0.000 | 0.000 |
| feature_07 | 0.00383 | 0.00080 | 0.32468 | 0.000 | 0.000 |
| feature_03 | 0.00124 | 0.00076 | 0.37144 | 0.000 | 0.000 |
| constant_score_control | 0.00075 | 0.00075 | 0.50000 | 0.025 | 0.125 |
| random_rank_control | 0.00110 | 0.00073 | 0.32999 | 0.000 | 0.000 |
| feature_08 | 0.00239 | 0.00059 | 0.20139 | 0.000 | 0.000 |
| shuffled_training_labels_control | 0.11343 | 0.00055 | 0.12100 | 0.000 | 0.000 |

## Failed attempts

[]

See STAGE6_EXPERIMENT_CATALOG.md for scope, methods, limitations and reproduction.
