import tempfile
import unittest
from pathlib import Path

import numpy as np

from genopred_v2.stage6 import (
    Stage6Error,
    deterministic_split,
    diffuse_scores,
    impute_from_partition,
    observed_label_metrics,
    predict_saved_model,
)


class Stage6BaselineTests(unittest.TestCase):
    def test_metrics_are_invariant_to_tied_gene_order(self):
        for labels in ([1,-1,-1,-1], [-1,-1,-1,1]):
            metrics = observed_label_metrics(np.array(labels), np.zeros(4), (1,2))
            self.assertEqual(.25, metrics['observed_label_auprc'])
            self.assertEqual(.5, metrics['observed_label_auroc'])
            self.assertEqual(.25, metrics['known_positive_recall_at_1'])
            self.assertEqual(.5, metrics['known_positive_recall_at_2'])

    def test_saved_model_restores_bias_and_rejects_wrong_feature_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'model.npz'
            np.savez(path, weights=np.array([2.,0.]), bias=np.array(-3.), medians=np.array([2.]), mean=np.array([1.,0.]), scale=np.ones(2), feature_names=np.array(['a','a__missing']))
            scores = predict_saved_model(path, np.array([[2.,0.],[999.,1.]]), ['a','a__missing'])
            np.testing.assert_allclose(scores, np.full(2,1/(1+np.exp(1))))
            with self.assertRaises(Stage6Error):
                predict_saved_model(path, np.zeros((1,2)), ['a__missing','a'])

    def test_split_is_deterministic_and_stratified(self):
        genes = [f"G{i}" for i in range(20)]
        labels = np.asarray([1] * 5 + [-1] * 15)
        arguments = (genes, labels, 7, {"train": 3, "validation": 1, "test": 1},
                     {"train": 0.6, "validation": 0.2, "test": 0.2})
        first = deterministic_split(*arguments)
        second = deterministic_split(*arguments)
        np.testing.assert_array_equal(first, second)
        self.assertEqual([3, 1, 1], [int(np.sum((first == i) & (labels == 1))) for i in range(3)])
        self.assertEqual([9, 3, 3], [int(np.sum((first == i) & (labels == -1))) for i in range(3)])

    def test_observed_metrics_are_perfect_for_separated_scores(self):
        result = observed_label_metrics(np.asarray([1, 1, -1, -1]), np.asarray([0.9, 0.8, 0.2, 0.1]))
        self.assertEqual(1.0, result["observed_label_auroc"])
        self.assertEqual(1.0, result["observed_label_auprc"])

    def test_diffusion_preserves_an_isolated_self_loop_node(self):
        scores = np.asarray([0.8, 0.2, 0.6])
        edges = np.asarray([[0, 1], [1, 0]])
        result, _ = diffuse_scores(scores, edges, np.ones(2), 0.5, 100, 1e-12, 1.0)
        self.assertAlmostEqual(0.6, result[2])

    def test_imputation_uses_only_the_permitted_partition(self):
        # Columns are two base features followed by their two missing flags.
        values = np.asarray([[1.0, 10.0, 0, 0], [3.0, 99.0, 0, 1], [100.0, 30.0, 1, 0]])
        result, medians = impute_from_partition(values, np.asarray([True, True, False]), 2)
        np.testing.assert_array_equal(medians, np.asarray([2.0, 10.0]))
        self.assertEqual(2.0, result[2, 0])
        self.assertEqual(10.0, result[1, 1])


if __name__ == "__main__":
    unittest.main()
