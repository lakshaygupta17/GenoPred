import unittest
import numpy as np
from genopred_v2.stage6 import Stage6Error, impute_from_partition
from genopred_v2.stage6_trees import fit_bags, score, predict_bundle


class TreeTests(unittest.TestCase):
    def test_deterministic_both_families_and_positive_scores(self):
        x = np.arange(60, dtype=float).reshape(30, 2)
        y = np.r_[np.ones(4), -np.ones(26)]
        for family in ['pu_tree', 'pu_random_forest']:
            a = fit_bags(x,y,family,2,5,41,bags=3)
            b = fit_bags(x,y,family,2,5,41,bags=3)
            np.testing.assert_array_equal(score(a,x), score(b,x))
            self.assertGreater(score(a,x[:4]).mean(), score(a,x[-4:]).mean())

    def test_partition_imputation_and_feature_order(self):
        x=np.array([[2.,0.],[4.,0.],[900.,0.],[999.,1.]])
        clean,med=impute_from_partition(x,np.array([True,True,False,False]),1)
        self.assertEqual(clean[-1,0],3.)
        bundle=dict(feature_names=['a','missing'],medians=med,models=[])
        with self.assertRaises(Stage6Error):
            predict_bundle(bundle,x,['missing','a'])

    def test_requires_both_pu_groups(self):
        with self.assertRaises(Stage6Error):
            fit_bags(np.zeros((5,2)),np.ones(5),'pu_tree',2,5,41)


if __name__=='__main__':
    unittest.main()
