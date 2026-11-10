import unittest
import numpy as np
from genopred_v2.stage6 import Stage6Error
from genopred_v2.stage6_benchmark import (
    registry, columns_for, fit_entry, predict, prepare, transform_saved,
)


class BroadBenchmarkTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(8)
        self.x=np.c_[rng.normal(size=(60,10)),np.zeros((60,10))]
        self.y=np.full(60,-1);self.y[[0,1,2,3,40,41]]=1
        self.mask=np.arange(60)<40
        self.names=[str(i) for i in range(20)]

    def test_registry_names_are_unique_and_all_subsets_valid(self):
        entries=registry()
        self.assertEqual(len(entries),len({e['name'] for e in entries}))
        for e in entries:
            cols=columns_for(e.get('columns'),20)
            self.assertTrue(cols and min(cols)>=0 and max(cols)<20)
        self.assertNotIn(6,columns_for('no_mutation',20))
        self.assertNotIn(16,columns_for('no_mutation',20))

    def test_heldout_labels_do_not_change_fitting(self):
        spec=dict(kind='bag',family='logistic_l2',level=0)
        a=fit_entry(self.x,self.y,self.mask,spec,41,self.names)
        altered=self.y.copy();altered[~self.mask]*=-1
        b=fit_entry(self.x,altered,self.mask,spec,41,self.names)
        np.testing.assert_array_equal(predict(a,self.x,self.names),predict(b,self.x,self.names))

    def test_heldout_values_do_not_change_preprocessor(self):
        spec=dict(kind='bag',family='logistic_l2',level=0)
        _,a=prepare(self.x,self.mask,spec,41)
        altered=self.x.copy();altered[~self.mask,:10]=100000
        _,b=prepare(altered,self.mask,spec,41)
        np.testing.assert_array_equal(a['medians'],b['medians'])
        np.testing.assert_array_equal(a['scaler'].mean_,b['scaler'].mean_)
        damaged=self.x.copy();damaged[50,0]=999;damaged[50,10]=1
        expected=(a['medians'][0]-a['scaler'].mean_[0])/a['scaler'].scale_[0]
        self.assertAlmostEqual(transform_saved(damaged,a)[50,0],expected)

    def test_graph_seeds_use_only_permitted_positives(self):
        state=fit_entry(self.x,self.y,self.mask,dict(kind='graph',level=0),41,self.names)
        self.assertEqual(state['initial'].sum(),4)
        self.assertFalse(state['initial'][~self.mask].any())

    def test_rejects_wrong_feature_order(self):
        state=fit_entry(self.x,self.y,self.mask,dict(kind='constant',level=0),41,self.names)
        with self.assertRaises(Stage6Error):predict(state,self.x,list(reversed(self.names)))

    def test_missing_feature_ranks_last_in_either_direction(self):
        for level in (0,1):
            state=fit_entry(self.x,self.y,self.mask,dict(kind='feature',level=level,index=0),41,self.names)
            altered=self.x.copy();altered[50,0]=1e12;altered[50,10]=1
            scores=predict(state,altered,self.names)
            self.assertEqual(scores.argmin(),50)


if __name__=='__main__': unittest.main()
