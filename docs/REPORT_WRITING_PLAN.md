# Report writing plan

Use FINAL_RESULTS.md as the frozen numerical source. Do not describe this project as cancer diagnosis or completed breast-to-pancreas transfer learning.

1. **Introduction:** pancreatic cancer driver prioritization, incomplete labels and why unlabelled genes cannot be assumed negative.
2. **Objectives:** reproducible molecular data integration, conservative evidence labels, model comparison and evaluation of limitations.
3. **Data and preprocessing:** source cohorts, gene identity, genome-coordinate QC, patient-level aggregation, 10 biological features plus missingness, provenance and corrected V1 defects.
4. **Methods:** linear nnPU, PU tree bagging, gradient boosting, LightGBM, mutation baseline, original exploratory screen, repeated nested gene-level evaluation and fixed-setting ablations.
5. **Results:** lead with nested evaluation; then ablations, MSK provisional evidence, secondary source-reconciled autopsy assessment and weak CRISPR support. Each table must name its gene universe, positive count and endpoint.
6. **Discussion:** why the old shallow-tree advantage was split-sensitive, mutation dominance, model disagreement and the difference between driver evidence and survival dependency.
7. **Limitations:** 17 positives, no true negatives, retrospective screening, external annotation and donor-overlap uncertainty, computational rather than experimentally verified labels, incomplete mouse assessment.
8. **Conclusion:** reproducible exploratory prioritization, not validated novel drivers or clinical accuracy.
9. **Future work:** licensed variant annotation, independent donor provenance, expert review, untouched functional/driver evidence and separately preregistered transfer-learning studies.
10. **Appendices:** full 60-experiment historical screen, all fold metrics, per-gene/variant audit, source discrepancy log, test and artifact verification.

Suggested figures: pipeline diagram; final_model_results.png; final_feature_ablations.png; existing external_functional_validation.png. Keep original single-split numbers in the exploratory-screen section or appendix.

Required wording: use observed-label AP rather than accuracy; provisional MSK evidence rather than confirmed new drivers; descriptive repeat range rather than confidence interval; external computational evidence rather than certified independent ground truth.
